#include "types.h"
#include "parser.h"
#include "floorplan.h"
#include "channel.h"
#include "router.h"
#include "sa_optimizer.h"
#include "legalize_loop.h"
#include "analytical_legalizer.h"
#include "output.h"
#include "config.h"
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

    auto route_now = [&]() {
        GlobalRouter gr;
        gr.init(d.blocks, d.channels, d.outline.max_width, d.outline.max_height);
        gr.route_all(d, 20);
    };

    // ── SA phase 1 ──────────────────────────────────────────────────────────
    SAOptimizer sa(fp, seed1);
    sa.time_limit_sec = sa1_time;
    sa.run();
    finalize_output();

    // ── SA phase 2 (fine) ───────────────────────────────────────────────────
    // Pack at base area to satisfy the fixed outline; the convergence loop below
    // sizes soft blocks to the real feedthrough.
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

    // Count edge blocks not sitting at their required boundary (each is an
    // evaluator FAIL).  The overlap-safe snap may leave a block unsnapped when
    // snapping it would collide; and FT expansion can shift the packing so an
    // edge block no longer reaches its corner.  Mirrors evaluator.check_edge_location.
    auto count_edge_fails = [&]() {
        int fails = 0;
        double W = d.outline.cur_width, H = d.outline.cur_height;
        for (int i : fp.edge_block_idx) {
            const auto& b = d.blocks[i];
            bool ok = false;
            for (const auto& loc : b.locations) {
                bool m = true;
                if (loc.find('T') != std::string::npos && std::abs(b.ly + b.height - H) > 1e-3) m = false;
                if (loc.find('B') != std::string::npos && std::abs(b.ly - 0.0)            > 1e-3) m = false;
                if (loc.find('L') != std::string::npos && std::abs(b.lx - 0.0)            > 1e-3) m = false;
                if (loc.find('R') != std::string::npos && std::abs(b.lx + b.width - W)    > 1e-3) m = false;
                if (m) { ok = true; break; }
            }
            if (!ok) fails++;
        }
        return fails;
    };

    // Penalty-aware score for ranking valid candidates.  The provided evaluator
    // counts each overflowing channel-direction / undersized soft block, so the
    // score counts overflows (heavily) then breaks ties by the proportional
    // overflow magnitude and finally the real contest cost (area + alpha*HPWL).
    // Edge-constraint FAILs are weighted above everything (a FAIL disqualifies).
    // Lexicographic priority: edge-FAIL > overflow penalty > cost.  FAILW must
    // dominate any sum of penalties; PEN must dominate any cost (area + α·HPWL)
    // change so a transition that trades cost for a penalty reduction always
    // wins.  Real cases: cost ~1e9, ~100 overflows max → PEN=1e10 gives a
    // single penalty (1e10) > any plausible cost swing; ~1000 max penalties
    // (1e13) << FAILW=1e15, so 1 FAIL still beats every overflow combined.
    const double FAILW = 1e15;
    const double PEN   = 1e10;
    auto route_and_score = [&]() {
        GlobalRouter gr;
        gr.init(d.blocks, d.channels, d.outline.max_width, d.outline.max_height);
        gr.route_all(d, 20);
        collect_ft();

        int pen = 0;
        double mag = 0.0;
        for (auto& ch : d.channels) {
            if (ch.nets_x > ch.cap_x() + 1e-3) { pen++; mag += ch.nets_x - ch.cap_x(); }
            if (ch.nets_y > ch.cap_y() + 1e-3) { pen++; mag += ch.nets_y - ch.cap_y(); }
        }
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT || fp.ft_nets[i] <= 0) continue;
            double required = d.blocks[i].get_target_area(fp.ft_nets[i]);
            double actual   = fp.W[i] * fp.H[i];
            if (actual < required - 1.0) { pen++; mag += (required - actual); }
        }
        return compute_final_cost(d) + FAILW * count_edge_fails() + PEN * pen + mag;
    };

    // ── Route-driven legalize loop ──────────────────────────────────────────
    // Replaces the legacy global-expand ft_iter loop with a hotspot-targeted
    // action cascade (rotate → displace → expand) under strict rollback.  The
    // loop never returns a worse layout than it was given, so the worst case
    // here is identical to "just take the SA output as-is."
    double best_score = 1e18;
    bool have_best = false;

    if (cfg::LEG_ENABLE) {
        printf("[LegalizeLoop] Starting with score %.3f\n", route_and_score());
        finalize_output();
        if (is_valid()) {
            LegalizeLoop loop(fp, d,
                              /*finalize=*/ finalize_output,
                              /*score   =*/ route_and_score,
                              /*is_valid=*/ is_valid,
                              seed1 ^ 0xACE5EEDu);
            loop.max_iters      = cfg::LEG_ITERS;
            loop.max_seconds    = cfg::LEG_TIME;
            loop.displace_tries = cfg::LEG_TRIES;
            best_score = loop.run();
            have_best  = is_valid();
        }

        // ── Final-pass analytical legalizer ──────────────────────────────────
        // Force-directed redistribution within the current compact outline.
        // Breaks the B*-tree by design (user-approved): operates directly on
        // (lx, ly) coordinates.  Strict rollback if it doesn't improve.
        if (have_best && cfg::ANA_ENABLE) {
            printf("[Analytical] Starting with score %.3f\n", best_score);
            // Snapshot the full pre-analytical state.
            auto snap_bst    = fp.bst.save();
            auto snap_W      = fp.W;
            auto snap_H      = fp.H;
            auto snap_loc    = fp.active_loc;
            auto snap_ft     = fp.ft_nets;
            auto snap_blocks = d.blocks;
            auto snap_chs    = d.channels;
            auto snap_paths  = d.paths;
            auto snap_outl   = d.outline;
            double snap_score = best_score;

            // Run the analytical pass.  It only mutates d.blocks[i].lx/ly for
            // non-pinned blocks (and ONLY clips within d.outline.cur_*).
            AnalyticalLegalizer ana(fp, d);
            ana.iterations = cfg::ANA_ITERS;
            ana.step_frac  = cfg::ANA_STEP;
            ana.repel_w    = cfg::ANA_REPEL;
            bool moved = ana.run();

            if (moved) {
                // Recompute channels from the new coords WITHOUT B*-tree pack
                // (the tree is now stale; coords are authoritative).  Outline
                // stays at the pre-analytical compact value — the legalizer
                // clipped within it so the max extent fits.
                d.channels = ChannelCalculator::compute(
                    d.blocks, d.outline.cur_width, d.outline.cur_height);

                if (is_valid()) {
                    double new_score = route_and_score();
                    if (new_score < snap_score - 1.0) {
                        best_score = new_score; // accept
                    } else {
                        // Rollback — analytical didn't help.
                        fp.bst.restore(snap_bst);
                        fp.W = snap_W; fp.H = snap_H;
                        fp.active_loc = snap_loc;
                        fp.ft_nets = snap_ft;
                        d.blocks = snap_blocks;
                        d.channels = snap_chs;
                        d.paths = snap_paths;
                        d.outline = snap_outl;
                    }
                } else {
                    // Invalid layout (overlap left after MTV resolution, or
                    // out-of-outline) — rollback.
                    fp.bst.restore(snap_bst);
                    fp.W = snap_W; fp.H = snap_H;
                    fp.active_loc = snap_loc;
                    fp.ft_nets = snap_ft;
                    d.blocks = snap_blocks;
                    d.channels = snap_chs;
                    d.paths = snap_paths;
                    d.outline = snap_outl;
                }
            }
        }
    } else {
        // Legacy global-expand convergence (kept for A/B testing via FP_LEG_ENABLE=0).
        // Records best valid state across rounds so an expansion that overshoots
        // doesn't leave d/fp in a worse state than an earlier round.
        Design best_d;
        std::vector<double> best_W, best_H;
        for (int ft_iter = 0; ft_iter <= 8; ft_iter++) {
            finalize_output();
            if (!is_valid()) break;
            double s = route_and_score();
            if (s < best_score) {
                best_score = s; best_d = d; best_W = fp.W; best_H = fp.H;
                have_best = true;
            }
            bool needs_expand = false;
            for (int i = 0; i < (int)d.blocks.size(); i++) {
                if (d.blocks[i].type != BlockType::SOFT || fp.ft_nets[i] <= 0) continue;
                double required = d.blocks[i].get_target_area(fp.ft_nets[i]);
                if (fp.W[i] * fp.H[i] < required - 1.0) { needs_expand = true; break; }
            }
            if (!needs_expand) break;
            fp.apply_ft_areas(false);
        }
        if (have_best) { d = best_d; fp.W = best_W; fp.H = best_H; }
    }

    if (!have_best) {
        // Emergency fallback: nothing valid — base area, clean snapped + routed.
        std::fill(fp.ft_nets.begin(), fp.ft_nets.end(), 0);
        fp.apply_ft_areas(true);
        finalize_output();
        route_now();
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

    cfg::load_from_env();

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
