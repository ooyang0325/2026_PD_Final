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
#include <cmath>

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

// Run one complete solve cycle. Returns final cost (without overflow penalty).
// Modifies d with the result (block positions, paths, channels).
// One full solve: coarse SA -> channelize -> fine SA -> route -> FT update -> reroute.
static double run_once(Design& d_in, double sa1_time, double sa2_time, unsigned seed1, unsigned seed2) {
    Design d = d_in; 

    Floorplan fp(d);
    SAOptimizer sa(fp, seed1);
    sa.time_limit_sec = sa1_time;
    sa.run();

    auto [_, _] = fp.pack();
    // 關鍵修改：強行將全局版面設定為大會給定的 MAX_OUTLINE
    d.outline.cur_width = d.outline.max_width;
    d.outline.cur_height = d.outline.max_height;
    
    // 將 Edge Blocks 精準歸位至 MAX_OUTLINE 的邊界，因為 SA 已經預留了完美空洞
    fp.finalize_edge_blocks(d.outline.max_width, d.outline.max_height);
    
    // 生成 Channel 時以 max_width/height 為畫布，這將吸收所有剩餘空間作為強大通道
    d.channels = ChannelCalculator::compute(d.blocks, d.outline.max_width, d.outline.max_height);

    // Phase 2: 細緻調整與預佈線
    SAOptimizer sa2(fp, seed2);
    sa2.time_limit_sec = sa2_time;
    sa2.T_init = 1e6;
    sa2.T_final = 1.0;
    sa2.update_ft_areas();
    fp.apply_ft_areas();
    sa2.run();

    fp.pack();
    fp.finalize_edge_blocks(d.outline.max_width, d.outline.max_height);
    d.channels = ChannelCalculator::compute(d.blocks, d.outline.max_width, d.outline.max_height);

    GlobalRouter gr;
    gr.init(d.blocks, d.channels, d.outline.max_width, d.outline.max_height);
    gr.route_all(d, 20); // 增加 Reroute 次數，保證 50 blocks 能繞通

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
    
    fp.apply_ft_areas(false);
    fp.pack();
    fp.finalize_edge_blocks(d.outline.max_width, d.outline.max_height);
    d.channels = ChannelCalculator::compute(d.blocks, d.outline.max_width, d.outline.max_height);

    GlobalRouter gr2;
    gr2.init(d.blocks, d.channels, d.outline.max_width, d.outline.max_height);
    gr2.route_all(d, 20);

    double cost = compute_final_cost(d);
    for (auto& ch : d.channels) if (ch.overflowed()) cost += 1e9;

    d_in = d; 
    return cost;
}

static unsigned mix_seed(unsigned base, unsigned idx) {
    // Simple mix to derive distinct seeds per worker/restart.
    uint64_t x = (uint64_t)base + 0x9e3779b97f4a7c15ULL * (uint64_t)(idx + 1);
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return (unsigned)(x & 0xffffffffu);
}

static std::pair<double, Design> run_search(const Design& d,
                                            double time_limit,
                                            unsigned seed_base) {
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    };

    double total_budget = time_limit * 0.95;
    const double min_restart = 10.0;

    Design best_d = d;
    double best_cost = 1e18;
    int restart = 0;

    while (elapsed() < total_budget) {
        double remaining = total_budget - elapsed();
        if (remaining < min_restart) break;

        double sa1_t = remaining * 0.75;
        double sa2_t = remaining * 0.20;

        unsigned seed1 = mix_seed(seed_base, (unsigned)(restart * 2));
        unsigned seed2 = mix_seed(seed_base, (unsigned)(restart * 2 + 1));

        Design trial = d;
        double cost = run_once(trial, sa1_t, sa2_t, seed1, seed2);

        if (cost < best_cost) {
            best_cost = cost;
            best_d = trial;
        }
        restart++;
    }

    return {best_cost, best_d};
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.csv> <output.cfg> [time_limit_sec]\n";
        return 1;
    }
    std::string in_path = argv[1];
    std::string out_path = argv[2];

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

    double time_limit = std::max(30.0, std::min(20 * std::pow((double)1.124, (double)d.blocks.size()), 7100.0)); // Scale time limit with block count

    std::cerr << "[Phase 2] Running search with time limit " << time_limit << "s\n";

    unsigned hc = std::thread::hardware_concurrency();
    int workers = (hc > 2) ? (int)hc - 2 : 1;
    if (workers < 1) workers = 1;

    std::vector<std::future<std::pair<double, Design>>> futures;
    futures.reserve((size_t)workers);

    for (int w = 0; w < workers; w++) {
        unsigned seed_base = 12345u + (unsigned)w * 101u;
        futures.push_back(std::async(std::launch::async, [=]() {
            return run_search(d, time_limit, seed_base);
        }));
    }

    Design best_d = d;
    double best_cost = 1e18;
    for (auto& fut : futures) {
        auto [cost, cand] = fut.get();
        if (cost < best_cost) {
            best_cost = cost;
            best_d = cand;
        }
    }

    std::cerr << "[Final] Writing output\n";
    OutputWriter::print_summary(best_d);
    OutputWriter::write(best_d, out_path);
    std::cerr << "Total time: " << elapsed() << "s\n";
    return 0;
}
