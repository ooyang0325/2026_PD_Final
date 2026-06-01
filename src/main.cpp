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
    bool enable_congestion = true;
    double cong_weight_max = 8000.0;
    int cong_bins = 10;
    double edge_penalty_coeff = 2e6;
    bool edge_use_best_location = true;
    bool enable_repair = true;
    double repair_sa_time = 12.0;
    double hard_center_weight_max = 1e5;
    double repair_hard_center_scale = 0.35;
    double sa1_slack_weight_max = 2.5e6;
    double sa2_slack_weight_max = 1.2e6;
    double sa2_slack_ramp_start = 0.3;
    double sa2_slack_min_frac = 0.0;
    double repair_slack_weight_max = 2.0e6;
    double repair_slack_min_frac = 0.8;
    double repair_cong_scale = 0.70;
    double repair_conn_weight_max = 8e4;
    int repair_conn_top_k = 12;
    int repair_trigger_open = 1;
    double sa2_ce_weight_max = 3e4;
    double repair_cost_improve_eps = 0.002;
    double sa3_time = 10.0;
    bool enable_sa3 = true;
    bool has_seed_base = false;
    unsigned seed_base = 12345u;
};

struct RepairScore {
    int total_fail_proxy = INT_MAX;
    int routing_open = INT_MAX;
    int edge_violation_proxy = INT_MAX;
    int outline_violation_proxy = INT_MAX;
    int overlap_violation_proxy = INT_MAX;
    int overflow_count = INT_MAX;
    double final_cost = 1e18;
};

struct WorkerResult {
    RepairScore score;
    Design d;
};

static bool better_repair_score(const RepairScore& a, const RepairScore& b) {
    if (a.total_fail_proxy != b.total_fail_proxy) return a.total_fail_proxy < b.total_fail_proxy;
    if (a.routing_open != b.routing_open) return a.routing_open < b.routing_open;
    if (a.edge_violation_proxy != b.edge_violation_proxy) return a.edge_violation_proxy < b.edge_violation_proxy;
    if (a.outline_violation_proxy != b.outline_violation_proxy) return a.outline_violation_proxy < b.outline_violation_proxy;
    if (a.overlap_violation_proxy != b.overlap_violation_proxy) return a.overlap_violation_proxy < b.overlap_violation_proxy;
    if (a.overflow_count != b.overflow_count) return a.overflow_count < b.overflow_count;
    return a.final_cost < b.final_cost;
}

static bool accept_repair_result(const RepairScore& before, const RepairScore& after,
                                 double cost_improve_eps) {
    if (better_repair_score(after, before)) return true;
    if (after.total_fail_proxy != before.total_fail_proxy) return false;
    if (after.routing_open != before.routing_open) return after.routing_open < before.routing_open;
    if (after.overflow_count != before.overflow_count) return after.overflow_count < before.overflow_count;
    return after.final_cost < before.final_cost * (1.0 - cost_improve_eps);
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

static int count_channel_overflow(const Design& d) {
    int c = 0;
    for (const auto& ch : d.channels) {
        if (ch.overflowed()) c++;
    }
    return c;
}

static int outline_violation_proxy(const Design& d) {
    int v = 0;
    if (d.outline.cur_width > d.outline.max_width + 1e-6) v++;
    if (d.outline.cur_height > d.outline.max_height + 1e-6) v++;
    for (const auto& b : d.blocks) {
        if (b.lx < -1e-6 || b.ly < -1e-6) v++;
        if (b.lx + b.width > d.outline.cur_width + 1e-6) v++;
        if (b.ly + b.height > d.outline.cur_height + 1e-6) v++;
    }
    return v;
}

static int overlap_violation_proxy(const Design& d) {
    int v = 0;
    for (int i = 0; i < (int)d.blocks.size(); i++) {
        const auto& a = d.blocks[i];
        for (int j = i + 1; j < (int)d.blocks.size(); j++) {
            const auto& b = d.blocks[j];
            double ox = std::max(0.0, std::min(a.lx + a.width, b.lx + b.width) - std::max(a.lx, b.lx));
            double oy = std::max(0.0, std::min(a.ly + a.height, b.ly + b.height) - std::max(a.ly, b.ly));
            if (ox > 1e-6 && oy > 1e-6) v++;
        }
    }
    return v;
}

static bool edge_location_satisfied(const Block& b, const std::string& loc, double ow, double oh) {
    const double tol = 1e-3;
    if (loc.find('L') != std::string::npos && std::abs(b.lx) > tol) return false;
    if (loc.find('R') != std::string::npos && std::abs((b.lx + b.width) - ow) > tol) return false;
    if (loc.find('B') != std::string::npos && std::abs(b.ly) > tol) return false;
    if (loc.find('T') != std::string::npos && std::abs((b.ly + b.height) - oh) > tol) return false;
    return true;
}

static int edge_violation_proxy(const Design& d) {
    int v = 0;
    double ow = d.outline.cur_width;
    double oh = d.outline.cur_height;
    for (const auto& b : d.blocks) {
        if (b.type != BlockType::EDGE || b.locations.empty()) continue;
        bool ok = false;
        for (const auto& loc : b.locations) {
            if (edge_location_satisfied(b, loc, ow, oh)) {
                ok = true;
                break;
            }
        }
        if (!ok) v++;
    }
    return v;
}

static double compute_final_cost(const Design& d);

static RepairScore score_design(const Design& d, int routing_open) {
    RepairScore s;
    s.routing_open = routing_open;
    s.edge_violation_proxy = edge_violation_proxy(d);
    s.outline_violation_proxy = outline_violation_proxy(d);
    s.overlap_violation_proxy = overlap_violation_proxy(d);
    s.overflow_count = count_channel_overflow(d);
    s.final_cost = compute_final_cost(d);
    s.total_fail_proxy = s.routing_open + s.edge_violation_proxy
                       + s.outline_violation_proxy + s.overlap_violation_proxy;
    return s;
}

static void log_repair_score(const char* tag, const RepairScore& s) {
    std::cerr << tag
              << " fail_proxy=" << s.total_fail_proxy
              << " open=" << s.routing_open
              << " edge=" << s.edge_violation_proxy
              << " outline=" << s.outline_violation_proxy
              << " overlap=" << s.overlap_violation_proxy
              << " overflow=" << s.overflow_count
              << " cost=" << s.final_cost << "\n";
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

// Short SA pass after routing failures: slack + optional hard-center / failed-pair conn.
static void run_repair_sa(Floorplan& fp, const SolverOptions& opt, unsigned seed,
                          const std::vector<int>& failed_conn, int routing_open,
                          double time_budget) {
    if (!opt.enable_repair || time_budget <= 0.0) return;
    const bool use_hard_center = routing_open > 0 && opt.hard_center_weight_max > 0.0;
    const bool use_slack = opt.repair_slack_weight_max > 0.0;
    const bool use_conn = opt.repair_conn_weight_max > 0.0 && !failed_conn.empty();
    if (!use_hard_center && !use_slack && !use_conn) return;

    std::cerr << "[Repair-SA] start budget=" << time_budget
              << "s routing_open=" << routing_open
              << " hard_center=" << (use_hard_center ? "on" : "off")
              << " slack=" << (use_slack ? "on" : "off") << "\n";

    SAOptimizer sa(fp, seed);
    sa.time_limit_sec = time_budget;
    sa.T_init = 5e7;
    sa.T_final = 1.0;
    sa.moves_per_temp = 150;
    sa.cong_weight_max = opt.enable_congestion ? opt.cong_weight_max * opt.repair_cong_scale : 0.0;
    sa.cong_bins = opt.cong_bins;
    sa.hard_center_weight_max = use_hard_center
        ? opt.hard_center_weight_max * opt.repair_hard_center_scale : 0.0;
    sa.slack_weight_max = use_slack ? opt.repair_slack_weight_max : 0.0;
    sa.slack_ramp_start = 0.0;
    sa.slack_min_frac = opt.repair_slack_min_frac;
    sa.repair_failed_conn_weight_max = use_conn ? opt.repair_conn_weight_max : 0.0;
    sa.repair_failed_conn_top_k = opt.repair_conn_top_k;
    fp.set_repair_failed_connections(failed_conn);
    sa.run();
    fp.pack();
    fp.set_repair_failed_connections({});
    std::cerr << "[Repair-SA] done\n";
}

// Post-route cost-only SA: swap moves, minimize area + HPWL proxy.
static void run_cost_sa3(Floorplan& fp, double time_budget, unsigned seed) {
    if (time_budget <= 0.0) return;
    std::cerr << "[SA3] cost-only start budget=" << time_budget << "s\n";
    SAOptimizer sa(fp, seed);
    sa.time_limit_sec = time_budget;
    sa.T_init = 2e6;
    sa.T_final = 1.0;
    sa.moves_per_temp = 120;
    sa.cost_only_mode = true;
    sa.swap_moves_only = true;
    sa.run();
    fp.pack();
    std::cerr << "[SA3] done\n";
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

static bool outline_exceeds_max(const Design& d) {
    return d.outline.cur_width > d.outline.max_width + 1e-6
        || d.outline.cur_height > d.outline.max_height + 1e-6;
}

// Run one complete solve cycle. Returns evaluator-aligned score.
static RepairScore run_once(Design& d_in, double sa1_time, double sa2_time,
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
    sa.moves_per_temp = 200;
    sa.cong_weight_max = opt.enable_congestion ? opt.cong_weight_max : 0.0;
    sa.cong_bins = opt.cong_bins;
    sa.slack_weight_max = opt.sa1_slack_weight_max;
    sa.slack_ramp_start = 0.0;
    sa.slack_min_frac = 0.7;
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
    sa2.moves_per_temp = 100;
    sa2.cong_weight_max = opt.enable_congestion ? opt.cong_weight_max : 0.0;
    sa2.cong_bins = opt.cong_bins;
    sa2.slack_weight_max = opt.sa2_slack_weight_max;
    sa2.slack_ramp_start = opt.sa2_slack_ramp_start;
    sa2.slack_min_frac = opt.sa2_slack_min_frac;
    sa2.slack_decay_late = true;
    sa2.slack_decay_start = 0.4;
    sa2.ce_weight_max = opt.sa2_ce_weight_max;
    sa2.ce_ramp_start = 0.4;
    sa2.ce_ramp_end = 0.7;
    sa2.ce_gate_feasible = true;
    sa2.cost_only_start = 0.7;
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
    RepairScore base_score = score_design(d, open_after_route);
    log_repair_score("[Repair-SA] baseline", base_score);
    const bool need_repair = open_after_route >= opt.repair_trigger_open
                          || outline_exceeds_max(d)
                          || base_score.overlap_violation_proxy > 0;
    if (!need_repair) {
        std::cerr << "[Repair-SA] skip (routing_open=" << open_after_route
                  << " outline_proxy=" << base_score.outline_violation_proxy << ")\n";
    } else if (!opt.enable_repair) {
        std::cerr << "[Repair-SA] skip (disabled, routing_open=" << open_after_route
                  << " outline_proxy=" << base_score.outline_violation_proxy << ")\n";
    } else {
        Design before_repair = d;
        RepairScore before_score = base_score;
        int open_before = open_after_route;
        double repair_budget = opt.repair_sa_time;
        if (open_after_route < opt.repair_trigger_open) repair_budget *= 0.5;
        log_failed_pairs(d, failed_conn, "[Repair-SA] before");
        run_repair_sa(fp, opt, mix_seed(seed2, 77), failed_conn, open_after_route, repair_budget);
        auto [tw_r, th_r] = fp.pack();
        d.outline.cur_width = tw_r;
        d.outline.cur_height = th_r;
        failed_conn.clear();
        reroute_design(d, &failed_conn);
        final_open = count_routing_open(failed_conn);
        RepairScore after_score = score_design(d, final_open);
        bool accept = accept_repair_result(before_score, after_score, opt.repair_cost_improve_eps);
        std::cerr << "[Repair-SA] routing_open " << open_before
                  << " -> " << final_open
                  << " (" << (accept ? "accepted" : "rollback") << ")\n";
        log_failed_pairs(d, failed_conn, "[Repair-SA] after");
        log_repair_score("[Repair-SA] after score", after_score);
        if (!accept) {
            d = before_repair;
            final_open = open_before;
            failed_conn.clear();
            reroute_design(d, &failed_conn);
            std::cerr << "[Repair-SA] reverted to baseline state\n";
            log_failed_pairs(d, failed_conn, "[Repair-SA] reverted");
            log_repair_score("[Repair-SA] reverted score", before_score);
        }
    }

    RepairScore final_score = score_design(d, final_open);
    if (opt.enable_sa3 && opt.sa3_time > 0.0 && final_score.total_fail_proxy == 0) {
        Design before_sa3 = d;
        RepairScore before_sa3_score = final_score;
        int open_before_sa3 = final_open;
        run_cost_sa3(fp, opt.sa3_time, mix_seed(seed2, 99));
        auto [tw_sa3, th_sa3] = fp.pack();
        d.outline.cur_width = tw_sa3;
        d.outline.cur_height = th_sa3;
        failed_conn.clear();
        reroute_design(d, &failed_conn);
        final_open = count_routing_open(failed_conn);
        RepairScore after_sa3_score = score_design(d, final_open);
        bool accept_sa3 = accept_repair_result(before_sa3_score, after_sa3_score, opt.repair_cost_improve_eps);
        std::cerr << "[SA3] open " << open_before_sa3 << " -> " << final_open
                  << " cost " << before_sa3_score.final_cost << " -> " << after_sa3_score.final_cost
                  << " (" << (accept_sa3 ? "accepted" : "rollback") << ")\n";
        if (!accept_sa3) {
            d = before_sa3;
            final_open = open_before_sa3;
            failed_conn.clear();
            reroute_design(d, &failed_conn);
        } else {
            final_score = after_sa3_score;
        }
    }

    d_in = d;
    return final_score;
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
    RepairScore best_score;
    int restart = 0;

    while (elapsed() < total_budget) {
        double remaining = total_budget - elapsed();
        if (remaining < min_restart) break;

        double sa3_reserve = (opt.enable_sa3 && opt.sa3_time > 0.0) ? opt.sa3_time : 0.0;
        double sa1_t = remaining * 0.75;
        double sa2_t = std::max(0.0, remaining * 0.20 - sa3_reserve);

        unsigned seed1 = mix_seed(seed_base, (unsigned)(restart * 2));
        unsigned seed2 = mix_seed(seed_base, (unsigned)(restart * 2 + 1));

        Design trial = d;
        RepairScore trial_score = run_once(trial, sa1_t, sa2_t, seed1, seed2, opt);

        if (better_repair_score(trial_score, best_score)) {
            best_score = trial_score;
            best_d = trial;
        }
        restart++;
    }

    return {best_score, best_d};
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.csv> <output.cfg> [time_limit_sec]\n"
                  << "       " << argv[0] << " <input.csv> <output.cfg> [--time SEC] [options]\n"
                  << "Options:\n"
                  << "  --enable-congestion | --disable-congestion\n"
                  << "  --cong-weight <value> --cong-bins <int>\n"
                  << "  --edge-penalty <value>\n"
                  << "  --edge-best-location | --edge-active-location\n"
                  << "  --enable-repair | --disable-repair  (Repair-SA on routing failure)\n"
                  << "  --repair-sa-time <sec> --hard-center-weight <value>\n"
                  << "  --sa1-slack-weight <value> --sa2-slack-weight <value>\n"
                  << "  --sa2-slack-ramp-start <frac> --repair-slack-weight <value>\n"
                  << "  --sa2-ce-weight <value> --sa3-time <sec> --disable-sa3\n"
                  << "  --repair-hard-center-scale <value> --repair-cong-scale <value>\n"
                  << "  --repair-conn-weight <value> --repair-conn-topk <int>\n"
                  << "  --seed-base <uint>\n"
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
        else if (arg == "--edge-penalty") opt.edge_penalty_coeff = std::stod(need_value(arg));
        else if (arg == "--edge-best-location") opt.edge_use_best_location = true;
        else if (arg == "--edge-active-location") opt.edge_use_best_location = false;
        else if (arg == "--enable-repair" || arg == "--enable-repair-sa") opt.enable_repair = true;
        else if (arg == "--disable-repair" || arg == "--disable-repair-sa") opt.enable_repair = false;
        else if (arg == "--repair-sa-time") opt.repair_sa_time = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--hard-center-weight") opt.hard_center_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--sa1-slack-weight") opt.sa1_slack_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--sa2-slack-weight") opt.sa2_slack_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--sa2-slack-ramp-start") opt.sa2_slack_ramp_start = std::max(0.0, std::min(0.95, std::stod(need_value(arg))));
        else if (arg == "--repair-slack-weight") opt.repair_slack_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--sa2-ce-weight") opt.sa2_ce_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--sa3-time") opt.sa3_time = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--disable-sa3") opt.enable_sa3 = false;
        else if (arg == "--repair-hard-center-scale") opt.repair_hard_center_scale = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--repair-cong-scale") opt.repair_cong_scale = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--repair-conn-weight") opt.repair_conn_weight_max = std::max(0.0, std::stod(need_value(arg)));
        else if (arg == "--repair-conn-topk") opt.repair_conn_top_k = std::max(1, std::stoi(need_value(arg)));
        else if (arg == "--seed-base") { opt.seed_base = (unsigned)std::stoul(need_value(arg)); opt.has_seed_base = true; }
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
        unsigned seed_base = (opt.has_seed_base ? opt.seed_base : 12345u) + (unsigned)w * 101u;
        futures.push_back(std::async(std::launch::async, [=]() {
            return run_search(d, time_limit, seed_base, opt);
        }));
    }

    Design best_d = d;
    RepairScore best_score;
    for (auto& fut : futures) {
        WorkerResult wr = fut.get();
        if (better_repair_score(wr.score, best_score)) {
            best_score = wr.score;
            best_d = wr.d;
        }
    }

    std::cerr << "[Final] Selected worker: fail_proxy=" << best_score.total_fail_proxy
              << " open=" << best_score.routing_open
              << " cost=" << best_score.final_cost << "\n";
    std::cerr << "[Final] Writing output\n";
    OutputWriter::print_summary(best_d);
    OutputWriter::write(best_d, out_path);
    std::cerr << "Total time: " << elapsed() << "s\n";
    return 0;
}
