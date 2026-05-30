#include "types.h"
#include "parser.h"
#include "floorplan.h"
#include "channel.h"
#include "router.h"
#include "sa_optimizer.h"
#include "output.h"
#include <iostream>
#include <string>
#include <chrono>
#include <future>
#include <thread>
#include <vector>
#include <cstdlib>
#include <climits>

struct SolverOptions {
    double time_limit = 7000.0;
    bool enable_connectivity = false;
    double conn_weight_max = 0.0;
    int conn_top_k = 30;
    bool enable_congestion = true;
    double cong_weight_max = 8000.0;
    int cong_bins = 10;
    double edge_penalty_coeff = 2e6;
    bool edge_use_best_location = true;
    bool enable_repair = true;
    double repair_sa_time = 12.0;
    double hard_center_weight_max = 3e5;
    int repair_trigger_open = 1;
};

struct WorkerResult {
    double cost = 1e18;
    int routing_open = INT_MAX;
    Design d;
};

static bool better_result(int open_a, double cost_a, int open_b, double cost_b) {
    if (open_a != open_b) return open_a < open_b;
    return cost_a < cost_b;
}

static unsigned mix_seed(unsigned base, unsigned idx) {
    uint64_t x = (uint64_t)base + 0x9e3779b97f4a7c15ULL * (uint64_t)(idx + 1);
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return (unsigned)(x & 0xffffffffu);
}

static int count_routing_open(const std::vector<int>& failed_conn) {
    return (int)failed_conn.size();
}

static void log_failed_pairs(const Design& d, const std::vector<int>& failed_conn,
                             const char* tag) {
    std::cerr << tag << " routing_open=" << failed_conn.size()
              << "/" << d.connections.size();
    int show = std::min(8, (int)failed_conn.size());
    for (int k = 0; k < show; k++) {
        int ci = failed_conn[k];
        const auto& c = d.connections[ci];
        std::cerr << " " << d.blocks[c.from].name << "-" << d.blocks[c.to].name;
    }
    if ((int)failed_conn.size() > show) std::cerr << " ...";
    std::cerr << "\n";
}

static void reroute_design(Design& d, std::vector<int>* failed_conn = nullptr) {
    double ow = d.outline.cur_width;
    double oh = d.outline.cur_height;
    d.channels = ChannelCalculator::compute(d.blocks, ow, oh);
    GlobalRouter gr;
    gr.init(d.blocks, d.channels, ow, oh);
    gr.route_all(d, 10, failed_conn);
}

// Short SA pass after routing failures: ramp hard-at-center penalty via sequence-pair moves.
static void run_repair_sa(Floorplan& fp, const SolverOptions& opt, unsigned seed) {
    if (!opt.enable_repair || opt.repair_sa_time <= 0.0) return;
    if (opt.hard_center_weight_max <= 0.0) return;

    std::cerr << "[Repair-SA] start budget=" << opt.repair_sa_time
              << "s hard_center_weight_max=" << opt.hard_center_weight_max << "\n";

    SAOptimizer sa(fp, seed);
    sa.time_limit_sec = opt.repair_sa_time;
    sa.T_init = 5e7;
    sa.T_final = 1.0;
    sa.cool_rate = 0.99;
    sa.moves_per_temp = 150;
    sa.conn_weight_max = 0.0;
    sa.cong_weight_max = opt.enable_congestion ? opt.cong_weight_max * 0.5 : 0.0;
    sa.cong_bins = opt.cong_bins;
    sa.hard_center_weight_max = opt.hard_center_weight_max;
    sa.run();
    fp.pack();
    std::cerr << "[Repair-SA] done\n";
}

// Derive a representative point for HPWL based on edge adjacency.
static std::pair<double, double> get_guiding_point(
    const std::string& r1_name, int e1,
    const std::string& r2_name, int e2,
    const Design& d) 
{
    double lx1=0, ly1=0, w1=0, h1=0;
    double lx2=0, ly2=0, w2=0, h2=0;
    
    auto lookup = [&](const std::string& name, double& lx, double& ly, double& w, double& h) {
        for(const auto& b : d.blocks) if(b.name == name) { lx=b.lx; ly=b.ly; w=b.width; h=b.height; return; }
        for(const auto& c : d.channels) if(c.name == name) { lx=c.lx; ly=c.ly; w=c.width; h=c.height; return; }
    };
    lookup(r1_name, lx1, ly1, w1, h1);
    lookup(r2_name, lx2, ly2, w2, h2);

    if (e1 == 3 && e2 == 1) { 
        return {lx1 + w1, (std::max(ly1, ly2) + std::min(ly1 + h1, ly2 + h2)) / 2.0};
    } else if (e1 == 1 && e2 == 3) {
        return {lx1, (std::max(ly1, ly2) + std::min(ly1 + h1, ly2 + h2)) / 2.0};
    } else if (e1 == 2 && e2 == 4) {
        return {(std::max(lx1, lx2) + std::min(lx1 + w1, lx2 + w2)) / 2.0, ly1 + h1};
    } else if (e1 == 4 && e2 == 2) {
        return {(std::max(lx1, lx2) + std::min(lx1 + w1, lx2 + w2)) / 2.0, ly1};
    }
    return {0, 0};
}

// Final cost uses routed paths when available, otherwise HPWL on block centers.
static double compute_final_cost(const Design& d) {
    double area = d.outline.cur_width * d.outline.cur_height;
    double hpwl = 0;
    
    if (!d.paths.empty()) {
        for (auto& path : d.paths) {
            if (path.segments.empty()) continue;
            double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
            bool has_pts = false;
            for (size_t i = 0; i + 1 < path.segments.size(); i++) {
                auto& s1 = path.segments[i];
                auto& s2 = path.segments[i+1];
                if (s1.rect_name == s2.rect_name) continue;
                auto pt = get_guiding_point(s1.rect_name, s1.edge_out, s2.rect_name, s2.edge_in, d);
                min_x = std::min(min_x, pt.first); max_x = std::max(max_x, pt.first);
                min_y = std::min(min_y, pt.second); max_y = std::max(max_y, pt.second);
                has_pts = true;
            }
            if (has_pts) hpwl += path.nets * ((max_x - min_x) + (max_y - min_y));
        }
    } else {
        // Fallback for SA process
        for (auto& conn : d.connections) {
            auto& ba = d.blocks[conn.from];
            auto& bb = d.blocks[conn.to];
            hpwl += conn.nets * (std::abs((bb.lx + bb.width/2) - (ba.lx + ba.width/2)) + 
                                 std::abs((bb.ly + bb.height/2) - (ba.ly + ba.height/2)));
        }
    }
    return area + d.alpha * hpwl;
}

// Run one complete solve cycle. Returns final cost and routing_open count.
// Modifies d with the result (block positions, paths, channels).
static std::pair<double, int> run_once(Design& d_in, double sa1_time, double sa2_time,
                       unsigned seed1, unsigned seed2, const SolverOptions& opt) {
    Design d = d_in; 

    Floorplan fp(d);
    fp.edge_penalty_coeff = opt.edge_penalty_coeff;
    fp.edge_use_best_location = opt.edge_use_best_location;
    // Phase 1: coarse SA for topology and outline feasibility.
    SAOptimizer sa(fp, seed1);
    sa.time_limit_sec = sa1_time;
    sa.T_init = 1e9;
    sa.T_final = 1e3;
    sa.cool_rate = 0.995;
    sa.moves_per_temp = 200;
    sa.conn_weight_max = opt.enable_connectivity ? opt.conn_weight_max : 0.0;
    sa.conn_top_k = opt.conn_top_k;
    sa.cong_weight_max = opt.enable_congestion ? opt.cong_weight_max : 0.0;
    sa.cong_bins = opt.cong_bins;
    sa.run();

    auto [tw, th] = fp.pack();
    d.outline.cur_width = tw;
    d.outline.cur_height = th;
    d.channels = ChannelCalculator::compute(d.blocks, tw, th);

        // Phase 2: FT-aware SA with updated soft-block areas.
        SAOptimizer sa2(fp, seed2);
    sa2.time_limit_sec = sa2_time;
    sa2.T_init = 1e6;
    sa2.T_final = 1.0;
    sa2.cool_rate = 0.99;
    sa2.moves_per_temp = 100;
    sa2.conn_weight_max = opt.enable_connectivity ? opt.conn_weight_max : 0.0;
    sa2.conn_top_k = opt.conn_top_k;
    sa2.cong_weight_max = opt.enable_congestion ? opt.cong_weight_max : 0.0;
    sa2.cong_bins = opt.cong_bins;
    sa2.update_ft_areas();
    fp.apply_ft_areas();
    sa2.run();

    auto [tw2, th2] = fp.pack();
    d.outline.cur_width = tw2;
    d.outline.cur_height = th2;
    d.channels = ChannelCalculator::compute(d.blocks, tw2, th2);

    // Phase 3: global routing with negotiated congestion.
    GlobalRouter gr;
    gr.init(d.blocks, d.channels, tw2, th2);
    gr.route_all(d, 10);

    // FT update: only extend blocks with actual FT nets (don't reset the rest).
    std::fill(fp.ft_nets.begin(), fp.ft_nets.end(), 0);
    for (auto& path : d.paths) {
        for (int si = 1; si < (int)path.segments.size() - 1; si++) {
            auto& seg = path.segments[si];
            if (seg.rect_name.substr(0,2) != "CH") {
                int bi = d.block_idx(seg.rect_name);
                if (bi >= 0 && d.blocks[bi].type == BlockType::SOFT)
                    fp.ft_nets[bi] += path.nets;
            }
        }
    }
    fp.apply_ft_areas();
    auto [tw3, th3] = fp.pack();
    d.outline.cur_width = tw3;
    d.outline.cur_height = th3;

    d.channels = ChannelCalculator::compute(d.blocks, tw3, th3);
    // Final reroute after FT-based resizing.
    GlobalRouter gr2;
    gr2.init(d.blocks, d.channels, tw3, th3);
    std::vector<int> failed_conn;
    gr2.route_all(d, 10, &failed_conn);
    int open_after_route = count_routing_open(failed_conn);
    log_failed_pairs(d, failed_conn, "[Route] after main flow");

    int final_open = open_after_route;
    if (open_after_route < opt.repair_trigger_open) {
        std::cerr << "[Repair-SA] skip (routing_open=" << open_after_route
                  << " < trigger=" << opt.repair_trigger_open << ")\n";
    } else if (!opt.enable_repair) {
        std::cerr << "[Repair-SA] skip (disabled, routing_open=" << open_after_route << ")\n";
    } else {
        int open_before = open_after_route;
        log_failed_pairs(d, failed_conn, "[Repair-SA] before");
        run_repair_sa(fp, opt, mix_seed(seed2, 77));
        auto [tw_r, th_r] = fp.pack();
        d.outline.cur_width = tw_r;
        d.outline.cur_height = th_r;
        failed_conn.clear();
        reroute_design(d, &failed_conn);
        final_open = count_routing_open(failed_conn);
        std::cerr << "[Repair-SA] routing_open " << open_before
                  << " -> " << final_open << "\n";
        log_failed_pairs(d, failed_conn, "[Repair-SA] after");
    }

    double cost = compute_final_cost(d);
    for (auto& ch : d.channels)
        if (ch.overflowed()) cost += 1e9;

    d_in = d;
    return {cost, final_open};
}

static WorkerResult run_search(const Design& d,
                                            double time_limit,
                                            unsigned seed_base,
                                            const SolverOptions& opt) {
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    };

    double total_budget = time_limit * 0.95;
    const double min_restart = 10.0;

    Design best_d = d;
    double best_cost = 1e18;
    int best_open = INT_MAX;
    int restart = 0;

    while (elapsed() < total_budget) {
        double remaining = total_budget - elapsed();
        if (remaining < min_restart) break;

        double sa1_t = remaining * 0.75;
        double sa2_t = remaining * 0.20;

        unsigned seed1 = mix_seed(seed_base, (unsigned)(restart * 2));
        unsigned seed2 = mix_seed(seed_base, (unsigned)(restart * 2 + 1));

        Design trial = d;
        auto [cost, open] = run_once(trial, sa1_t, sa2_t, seed1, seed2, opt);

        if (better_result(open, cost, best_open, best_cost)) {
            best_open = open;
            best_cost = cost;
            best_d = trial;
        }
        restart++;
    }

    return {best_cost, best_open, best_d};
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.csv> <output.cfg> [time_limit_sec]\n"
                  << "       " << argv[0] << " <input.csv> <output.cfg> [--time SEC] [options]\n"
                  << "Options:\n"
                  << "  --enable-congestion | --disable-congestion\n"
                  << "  --cong-weight <value> --cong-bins <int>\n"
                  << "  --enable-connectivity | --disable-connectivity\n"
                  << "  --conn-weight <value> --conn-topk <int>\n"
                  << "  --edge-penalty <value>\n"
                  << "  --edge-best-location | --edge-active-location\n"
                  << "  --enable-repair | --disable-repair  (Repair-SA on routing failure)\n"
                  << "  --repair-sa-time <sec> --hard-center-weight <value>\n"
                  << "  --repair-trigger-open <int>\n";
        return 1;
    }
    std::string in_path = argv[1];
    std::string out_path = argv[2];
    SolverOptions opt;

    int idx = 3;
    if (idx < argc && argv[idx][0] != '-') {
        opt.time_limit = std::stod(argv[idx]);
        idx++;
    }
    while (idx < argc) {
        std::string arg = argv[idx++];
        auto need_value = [&](const std::string& name) {
            if (idx >= argc) {
                std::cerr << "ERROR: Missing value for " << name << "\n";
                std::exit(1);
            }
            return std::string(argv[idx++]);
        };
        if (arg == "--time") opt.time_limit = std::stod(need_value(arg));
        else if (arg == "--enable-congestion") opt.enable_congestion = true;
        else if (arg == "--disable-congestion") opt.enable_congestion = false;
        else if (arg == "--cong-weight") opt.cong_weight_max = std::stod(need_value(arg));
        else if (arg == "--cong-bins") opt.cong_bins = std::max(2, std::stoi(need_value(arg)));
        else if (arg == "--enable-connectivity") opt.enable_connectivity = true;
        else if (arg == "--disable-connectivity") opt.enable_connectivity = false;
        else if (arg == "--conn-weight") opt.conn_weight_max = std::stod(need_value(arg));
        else if (arg == "--conn-topk") opt.conn_top_k = std::max(1, std::stoi(need_value(arg)));
        else if (arg == "--edge-penalty") opt.edge_penalty_coeff = std::stod(need_value(arg));
        else if (arg == "--edge-best-location") opt.edge_use_best_location = true;
        else if (arg == "--edge-active-location") opt.edge_use_best_location = false;
        else if (arg == "--enable-repair" || arg == "--enable-repair-sa") opt.enable_repair = true;
        else if (arg == "--disable-repair" || arg == "--disable-repair-sa") opt.enable_repair = false;
        else if (arg == "--repair-sa-time") opt.repair_sa_time = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--hard-center-weight") opt.hard_center_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--repair-trigger-open") opt.repair_trigger_open = std::max(0, std::stoi(need_value(arg)));
        else {
            std::cerr << "ERROR: Unknown option " << arg << "\n";
            return 1;
        }
    }
    double time_limit = opt.time_limit;

    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    };

    std::cerr << "[Phase 1] Parsing " << in_path << "\n";
    Design d = Parser::load_csv(in_path);
    if (d.blocks.empty()) {
        std::cerr << "ERROR: No blocks loaded.\n";
        return 1;
    }
    std::cerr << "  Loaded " << d.blocks.size() << " blocks, "
              << d.connections.size() << " connections\n"
              << "  Outline max: " << d.outline.max_width << " x " << d.outline.max_height << "\n"
              << "  alpha = " << d.alpha << "\n";

    unsigned hc = std::thread::hardware_concurrency();
//    int workers = (hc > 2) ? (int)hc - 2 : 1;
    int workers = hc / 2;
    if (workers < 1) workers = 1;

    std::vector<std::future<WorkerResult>> futures;
    futures.reserve((size_t)workers);

    for (int w = 0; w < workers; w++) {
        unsigned seed_base = 12345u + (unsigned)w * 101u;
        futures.push_back(std::async(std::launch::async, [=]() {
            return run_search(d, time_limit, seed_base, opt);
        }));
    }

    Design best_d = d;
    double best_cost = 1e18;
    int best_open = INT_MAX;
    for (auto& fut : futures) {
        WorkerResult wr = fut.get();
        if (better_result(wr.routing_open, wr.cost, best_open, best_cost)) {
            best_open = wr.routing_open;
            best_cost = wr.cost;
            best_d = wr.d;
        }
    }

    std::cerr << "[Final] Selected worker result: routing_open=" << best_open
              << " cost=" << best_cost << "\n";
    std::cerr << "[Final] Writing output\n";
    OutputWriter::print_summary(best_d);
    OutputWriter::write(best_d, out_path);
    std::cerr << "Total time: " << elapsed() << "s\n";
    return 0;
}
