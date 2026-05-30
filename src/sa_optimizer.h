#pragma once
#include "floorplan.h"
#include "channel.h"
#include "router.h"
#include <random>
#include <cmath>
#include <iostream>
#include <chrono>

class SAOptimizer {
public:
    Floorplan& fp;
    GlobalRouter router;
    std::mt19937 rng;

    double T_init = 1e8;
    double T_final = 10.0;
    double cool_rate = 0.99;
    int moves_per_temp = 200;
    double time_limit_sec = 6000.0;
    // Connectivity penalty is currently disabled for stability testing.
    double conn_weight_max = 0.0;
    int conn_top_k = 30;
    double cong_weight_max = 3500.0;
    int cong_bins = 10;
    double hard_center_weight_max = 0.0;

    SAOptimizer(Floorplan& fp_, unsigned seed = 42)
        : fp(fp_), rng(seed) {}

    // Time-budgeted SA with coarse->fine transition.
    // Hot loop uses O(1) incremental undo (no per-iteration vector copies).
    void run() {
        auto t0 = std::chrono::steady_clock::now();
        auto elapsed = [&]() {
            return std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
        };

        int n = fp.sp.n;
        double alpha = fp.d.alpha;
        double max_w = fp.d.outline.max_width;
        double max_h = fp.d.outline.max_height;

        fp.apply_ft_areas();
        auto [tw0, th0] = fp.eval();
        double cur_cost = fp.compute_cost(tw0, th0, alpha, max_w, max_h, 0.0, conn_top_k, 0.0, cong_bins, 0.0);

        // Best-state checkpoint (saved only on improvement, not per iter)
        SequencePair::State best_state = fp.sp.save();
        std::vector<double> best_W = fp.W, best_H = fp.H;
        std::vector<int>    best_loc = fp.active_loc;
        double best_cost = cur_cost;

        int iter = 0;
        int phase = 0; // 0=coarse, 1=fine
        double T = T_init;

        std::uniform_int_distribution<int> move_type(1, 6);
        std::uniform_int_distribution<int> bidx(0, n-1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        while (elapsed() < time_limit_sec) {
            double t_frac = std::min(1.0, elapsed() / time_limit_sec);
            double conn_weight = conn_weight_max * std::max(0.0, (t_frac - 0.2) / 0.8);
            double cong_weight = cong_weight_max * std::max(0.0, (t_frac - 0.3) / 0.7);
            double hard_center_weight = hard_center_weight_max * t_frac;
            if (phase == 0 && T < T_init * 1e-4) {
                phase = 1;
                update_ft_areas();
                fp.apply_ft_areas();
                auto [tw, th] = fp.eval();
                cur_cost = fp.compute_cost(tw, th, alpha, max_w, max_h, conn_weight, conn_top_k, cong_weight, cong_bins, hard_center_weight);
                if (cur_cost < best_cost) {
                    best_cost = cur_cost;
                    best_state = fp.sp.save();
                    best_W = fp.W; best_H = fp.H;
                    best_loc = fp.active_loc;
                }
            }

            for (int m = 0; m < moves_per_temp; m++) {
                if (elapsed() >= time_limit_sec) break;

                int mv = move_type(rng);
                int bi = bidx(rng);
                int bj = (mv <= 3) ? bidx(rng) : 0;
                if (mv <= 3) {
                    while (bj == bi) bj = bidx(rng);
                }

                // Incremental undo state (stack-only, no heap alloc)
                double saved_W = 0, saved_H = 0;
                int saved_loc = 0;
                bool changed = false;

                switch (mv) {
                    case 1: fp.sp.swap_gp(bi, bj); changed = true; break;
                    case 2: fp.sp.swap_gm(bi, bj); changed = true; break;
                    case 3: fp.sp.swap_both(bi, bj); changed = true; break;
                    case 4: {
                        // Rotate: swap W and H when the aspect ratio is within allowed limits
                        if (fp.d.blocks[bi].type == BlockType::SOFT) {
                            double target_area = fp.d.blocks[bi].get_target_area(fp.ft_nets[bi]);
                            double new_w = fp.H[bi], new_h = fp.W[bi];
                            double ar = (new_h > 0) ? new_w / new_h : 1.0;
                            double mn = fp.d.blocks[bi].min_ar;
                            double mx = fp.d.blocks[bi].max_ar;
                            if (ar >= mn - 1e-6 && ar <= mx + 1e-6) {
                                saved_W = fp.W[bi]; saved_H = fp.H[bi];
                                // Round up to two decimal places
                                fp.W[bi] = std::ceil(new_w * 100.0) / 100.0;
                                fp.H[bi] = std::ceil((target_area / fp.W[bi]) * 100.0) / 100.0;
                                changed = true;
                            }
                        }
                        break;
                    }
                    case 5: {
                        // Resize
                        if (fp.d.blocks[bi].type == BlockType::SOFT) {
                            double target_area = fp.d.blocks[bi].get_target_area(fp.ft_nets[bi]);
                            double mn = fp.d.blocks[bi].min_ar;
                            double mx = fp.d.blocks[bi].max_ar;
                            std::uniform_real_distribution<double> ar_dist(mn, mx);
                            double new_ar = ar_dist(rng);
                            
                            saved_W = fp.W[bi]; saved_H = fp.H[bi];
                            double raw_w = std::sqrt(target_area * new_ar);
                            // Round up to two decimal places
                            fp.W[bi] = std::ceil(raw_w * 100.0) / 100.0;
                            fp.H[bi] = std::ceil((target_area / fp.W[bi]) * 100.0) / 100.0;
                            changed = true;
                        }
                        break;
                    }
                    case 6: {
                        // Flip active location for an edge block
                        if (fp.d.blocks[bi].type == BlockType::EDGE &&
                            fp.d.blocks[bi].locations.size() > 1) {
                            saved_loc = fp.active_loc[bi];
                            fp.flip_location(bi);
                            changed = true;
                        }
                        break;
                    }
                }

                if (!changed) continue;

                auto [ntw, nth] = fp.eval();
                double new_cost = fp.compute_cost(ntw, nth, alpha, max_w, max_h, conn_weight, conn_top_k, cong_weight, cong_bins, hard_center_weight);
                double delta = new_cost - cur_cost;

                bool accept = (delta <= 0) || (unit(rng) < std::exp(-delta / T));
                if (accept) {
                    cur_cost = new_cost;
                    if (new_cost < best_cost) {
                        best_cost = new_cost;
                        best_state = fp.sp.save();
                        best_W = fp.W;
                        best_H = fp.H;
                        best_loc = fp.active_loc;
                    }
                } else {
                    // O(1) incremental undo
                    switch (mv) {
                        case 1: fp.sp.swap_gp(bi, bj); break;
                        case 2: fp.sp.swap_gm(bi, bj); break;
                        case 3: fp.sp.swap_both(bi, bj); break;
                        case 4:
                        case 5:
                            fp.W[bi] = saved_W; fp.H[bi] = saved_H;
                            break;
                        case 6:
                            fp.active_loc[bi] = saved_loc;
                            break;
                    }
                }
            }

            // Time-proportional cooling: T(t) = T_init * (T_final/T_init)^(t/total)
            t_frac = std::min(1.0, elapsed() / time_limit_sec);
            T = T_init * std::pow(T_final / T_init, t_frac);
            iter++;
            if (iter % 5000 == 0) {
                std::cerr << "[SA] iter=" << iter << " T=" << T
                          << " best=" << best_cost
                          << " time=" << elapsed() << "s\n";
            }
        }

        // Restore best state and finalize d.blocks via pack() (which calls commit()).
        fp.sp.restore(best_state);
        fp.W = best_W;
        fp.H = best_H;
        fp.active_loc = best_loc;
        fp.pack();

        std::cerr << "[SA] Done. Best cost=" << best_cost
                  << " in " << elapsed() << "s\n";
    }

    // Bounding-box heuristic for FT estimation. Reads from sp.x/sp.y/W/H so
    // it works even if d.blocks is stale during SA.
    void update_ft_areas() {
        std::fill(fp.ft_nets.begin(), fp.ft_nets.end(), 0);
        const auto& spx = fp.sp.x;
        const auto& spy = fp.sp.y;
        const auto& W = fp.W;
        const auto& H = fp.H;

        // Pre-collect soft block indices to avoid type check in inner loop
        std::vector<int> soft_idx;
        soft_idx.reserve(fp.d.blocks.size());
        for (int i = 0; i < (int)fp.d.blocks.size(); i++)
            if (fp.d.blocks[i].type == BlockType::SOFT) soft_idx.push_back(i);

        for (auto& conn : fp.d.connections) {
            int a = conn.from, b = conn.to;
            double cax = spx[a] + W[a] * 0.5, cay = spy[a] + H[a] * 0.5;
            double cbx = spx[b] + W[b] * 0.5, cby = spy[b] + H[b] * 0.5;
            double min_x = std::min(cax, cbx);
            double max_x = std::max(cax, cbx);
            double min_y = std::min(cay, cby);
            double max_y = std::max(cay, cby);

            for (int i : soft_idx) {
                if (i == a || i == b) continue;
                double bx0 = spx[i], bx1 = spx[i] + W[i];
                double by0 = spy[i], by1 = spy[i] + H[i];
                if (bx0 < max_x && bx1 > min_x && by0 < max_y && by1 > min_y) {
                    fp.ft_nets[i] += conn.nets;
                }
            }
        }
    }
};
