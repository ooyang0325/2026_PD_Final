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
// Flow: coarse SA -> fine SA -> 4x (route -> FT expand via full repack -> finalize).
static double run_once(Design& d_in, double sa1_time, double sa2_time, unsigned seed1, unsigned seed2) {
    Design d = d_in;
    Floorplan fp(d);

    // Full finalize: pack + snap edge blocks to the compact boundary + recompute channels.
    // Only call this for the final output step to avoid snap-induced overlaps during
    // intermediate FT repacks (the sequence pair does not guarantee edge blocks stay
    // at the extremes after soft-block size changes).
    auto finalize_output = [&]() {
        auto [tw, th] = fp.pack();
        double aw = std::min(tw, d.outline.max_width);
        double ah = std::min(th, d.outline.max_height);
        d.outline.cur_width  = aw;
        d.outline.cur_height = ah;
        fp.finalize_edge_blocks(aw, ah);
        d.channels = ChannelCalculator::compute(d.blocks, aw, ah);
    };

    // Collect FT net loads from current routing paths into fp.ft_nets.
    auto collect_ft = [&]() {
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
    };

    // ── SA phase 1 ──────────────────────────────────────────────────────────
    SAOptimizer sa(fp, seed1);
    sa.time_limit_sec = sa1_time;
    sa.run();
    finalize_output();

    // ── SA phase 2 (fine) ───────────────────────────────────────────────────
    // No FT pre-expansion: pack at base area to satisfy the fixed outline; the
    // convergence loop below sizes soft blocks to the real feedthrough.
    SAOptimizer sa2(fp, seed2);
    sa2.time_limit_sec = sa2_time;
    sa2.run();
    finalize_output();

    // Block overlap test on the committed layout.
    auto has_overlap = [&]() {
        int nb = (int)d.blocks.size();
        for (int i = 0; i < nb; i++) {
            for (int j = i+1; j < nb; j++) {
                double ox = std::min(d.blocks[i].lx + d.blocks[i].width,
                                     d.blocks[j].lx + d.blocks[j].width)
                          - std::max(d.blocks[i].lx, d.blocks[j].lx);
                double oy = std::min(d.blocks[i].ly + d.blocks[i].height,
                                     d.blocks[j].ly + d.blocks[j].height)
                          - std::max(d.blocks[i].ly, d.blocks[j].ly);
                if (ox > 0.01 && oy > 0.01) return true;
            }
        }
        return false;
    };

    // Full validity: every block inside the (compact) outline AND non-overlapping.
    // Any violation here would be an evaluator FAIL, so such layouts are never kept.
    auto is_valid = [&]() {
        double W = d.outline.cur_width, H = d.outline.cur_height;
        if (W > d.outline.max_width + 1e-3 || H > d.outline.max_height + 1e-3) return false;
        for (auto& b : d.blocks) {
            if (b.lx < -1e-3 || b.ly < -1e-3 ||
                b.lx + b.width  > W + 1e-3 ||
                b.ly + b.height > H + 1e-3) return false;
        }
        return !has_overlap();
    };

    // Penalty-aware score for ranking valid candidates: real contest cost plus a
    // heavy term per channel-overflow / feedthrough-overflow so the search drives
    // those penalties to zero before optimizing area+HPWL.
    const double PEN = 1e7;
    auto route_and_score = [&]() {
        GlobalRouter gr;
        gr.init(d.blocks, d.channels, d.outline.max_width, d.outline.max_height);
        gr.route_all(d, 20);
        collect_ft();

        int pen = 0;
        for (auto& ch : d.channels) {
            if (ch.nets_x > ch.cap_x() + 1e-3) pen++;
            if (ch.nets_y > ch.cap_y() + 1e-3) pen++;
        }
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT || fp.ft_nets[i] <= 0) continue;
            double required = d.blocks[i].get_target_area(fp.ft_nets[i]);
            if (fp.W[i] * fp.H[i] < required - 1.0) pen++;
        }
        return compute_final_cost(d) + PEN * pen;
    };

    // ── Robust FT convergence ───────────────────────────────────────────────
    // Each round: (re)snap + route the current layout; if it is fully valid, keep
    // it as a candidate (best score wins).  Then, if any soft block is undersized
    // for the observed feedthrough, expand and repack for another round.  Because
    // only valid layouts are ever recorded, the returned solution can never be an
    // evaluator FAIL — at worst it carries some channel/FT penalties.
    Design best_d; std::vector<double> best_W, best_H;
    double best_score = 1e18; bool have_best = false;

    for (int ft_iter = 0; ft_iter <= 8; ft_iter++) {
        finalize_output();           // pack + snap edge blocks + channels
        if (!is_valid()) break;      // expansion broke the fit; keep best-so-far

        double s = route_and_score();
        if (s < best_score) {
            best_score = s; best_d = d; best_W = fp.W; best_H = fp.H; have_best = true;
        }

        bool needs_expand = false;
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT || fp.ft_nets[i] <= 0) continue;
            double required = d.blocks[i].get_target_area(fp.ft_nets[i]);
            if (fp.W[i] * fp.H[i] < required - 1.0) { needs_expand = true; break; }
        }
        if (!needs_expand) break;

        fp.apply_ft_areas(false);    // grow undersized soft blocks for next round
    }

    if (have_best) {
        d = best_d; fp.W = best_W; fp.H = best_H;
    } else {
        // Emergency fallback: no expanded layout was valid.  Shrink soft blocks
        // back to base area (most compact) and emit a clean snapped+routed layout.
        std::fill(fp.ft_nets.begin(), fp.ft_nets.end(), 0);
        fp.apply_ft_areas(true);
        finalize_output();
        GlobalRouter gr;
        gr.init(d.blocks, d.channels, d.outline.max_width, d.outline.max_height);
        gr.route_all(d, 20);
        best_score = route_and_score();
    }

    d_in = d;

    // Validity gate: an out-of-bounds/overlapping layout is an evaluator FAIL.
    // Return a prohibitive-but-discriminating cost (base 1e15 plus the total
    // boundary/overlap violation) so that (a) any valid layout always wins, and
    // (b) among invalid layouts the *least-bad* one is kept — never the trivial
    // all-at-origin input design.
    if (!is_valid()) {
        double Wc = d.outline.cur_width, Hc = d.outline.cur_height;
        double viol = 0.0;
        int nb = (int)d.blocks.size();
        for (auto& b : d.blocks) {
            if (b.lx < 0) viol += -b.lx;
            if (b.ly < 0) viol += -b.ly;
            if (b.lx + b.width  > Wc) viol += b.lx + b.width  - Wc;
            if (b.ly + b.height > Hc) viol += b.ly + b.height - Hc;
        }
        for (int i = 0; i < nb; i++)
            for (int j = i + 1; j < nb; j++) {
                double ox = std::min(d.blocks[i].lx + d.blocks[i].width,  d.blocks[j].lx + d.blocks[j].width)
                          - std::max(d.blocks[i].lx, d.blocks[j].lx);
                double oy = std::min(d.blocks[i].ly + d.blocks[i].height, d.blocks[j].ly + d.blocks[j].height)
                          - std::max(d.blocks[i].ly, d.blocks[j].ly);
                if (ox > 0 && oy > 0) viol += ox + oy;
            }
        return 1e15 + viol;
    }
    return best_score;
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
    if(argv[3] != nullptr) time_limit = std::stod(argv[3]);
    
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
