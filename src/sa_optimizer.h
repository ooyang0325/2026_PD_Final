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
    double time_limit_sec = 6000.0; // 100 min, leave headroom for routing

    SAOptimizer(Floorplan& fp_, unsigned seed = 42)
        : fp(fp_), rng(seed) {}

    // Time-budgeted SA with coarse->fine transition and FT-aware updates.
    void run() {
        auto t0 = std::chrono::steady_clock::now();
        auto elapsed = [&]() {
            auto now = std::chrono::steady_clock::now();
            return std::chrono::duration<double>(now - t0).count();
        };

        int n = fp.sp.n;
        double& alpha = fp.d.alpha;
        double max_w = fp.d.outline.max_width;
        double max_h = fp.d.outline.max_height;

        // Initial pack and baseline cost
        fp.apply_ft_areas();
        auto [tw, th] = fp.pack();
        double cur_cost = fp.compute_cost(tw, th, alpha, max_w, max_h);

        SequencePair::State best_state = fp.sp.save();
        std::vector<double> best_W = fp.W, best_H = fp.H;
        std::vector<int> best_loc = fp.active_loc;
        double best_cost = cur_cost;

        // Time-proportional cooling: T(t) = T_init * (T_final/T_init)^(t/time_limit)
        // This ensures we explore broadly early and refine near the end.
        int iter = 0;
        int phase = 0; // 0=coarse, 1=fine
        double T = T_init;

        std::uniform_int_distribution<int> move_type(1, 6);
        std::uniform_int_distribution<int> bidx(0, n-1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        while (elapsed() < time_limit_sec) {
            // Switch to fine phase when temperature drops sufficiently
            if (phase == 0 && T < T_init * 1e-4) {
                phase = 1;
                // Update FT areas (run quick routing)
                update_ft_areas();
                fp.apply_ft_areas();
            }

            for (int m = 0; m < moves_per_temp; m++) {
                if (elapsed() >= time_limit_sec) break;

                // Choose a move
                int mv = move_type(rng);
                int bi = bidx(rng), bj = bidx(rng);
                while (bj == bi) bj = bidx(rng);

                SequencePair::State prev_sp = fp.sp.save();
                std::vector<double> prev_W = fp.W, prev_H = fp.H;
                std::vector<int> prev_loc = fp.active_loc;

                bool changed = false;
                if (mv <= 3) {
                    // SP perturbation (topology change)
                    if (mv == 1) fp.sp.swap_gp(bi, bj);
                    else if (mv == 2) fp.sp.swap_gm(bi, bj);
                    else fp.sp.swap_both(bi, bj);
                    changed = true;
                } else if (mv == 4) {
                    // Rotate a soft block (swap W/H if within AR range)
                    if (fp.d.blocks[bi].type == BlockType::SOFT) {
                        double new_w = fp.H[bi], new_h = fp.W[bi];
                        double ar = new_w / new_h;
                        double mn = fp.d.blocks[bi].min_ar;
                        double mx = fp.d.blocks[bi].max_ar;
                        if (ar >= mn - 1e-6 && ar <= mx + 1e-6) {
                            fp.W[bi] = new_w; fp.H[bi] = new_h;
                            changed = true;
                        }
                    }
                } else if (mv == 5) {
                    // Resize: adjust AR of a soft block within bounds
                    if (fp.d.blocks[bi].type == BlockType::SOFT) {
                        double area_i = fp.W[bi] * fp.H[bi];
                        double mn = fp.d.blocks[bi].min_ar;
                        double mx = fp.d.blocks[bi].max_ar;
                        std::uniform_real_distribution<double> ar_dist(mn, mx);
                        double new_ar = ar_dist(rng);
                        fp.W[bi] = std::sqrt(area_i * new_ar);
                        fp.H[bi] = area_i / fp.W[bi];
                        changed = true;
                    }
                } else {
                    // Flip active location for an edge block
                    if (fp.d.blocks[bi].type == BlockType::EDGE &&
                        fp.d.blocks[bi].locations.size() > 1) {
                        fp.flip_location(bi);
                        changed = true;
                    }
                }

                if (!changed) continue;

                auto [ntw, nth] = fp.pack();
                double new_cost = fp.compute_cost(ntw, nth, alpha, max_w, max_h);
                double delta = new_cost - cur_cost;

                bool accept = (delta < 0) || (unit(rng) < std::exp(-delta / T));
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
                    // Revert
                    fp.sp.restore(prev_sp);
                    fp.W = prev_W;
                    fp.H = prev_H;
                    fp.active_loc = prev_loc;
                }
            }

            // Time-based temperature update
            double t_frac = std::min(1.0, elapsed() / time_limit_sec);
            T = T_init * std::pow(T_final / T_init, t_frac);
            iter++;
            if (iter % 5000 == 0) {
                std::cerr << "[SA] iter=" << iter << " T=" << T
                          << " best=" << best_cost
                          << " time=" << elapsed() << "s\n";
            }
        }

        // Restore best state found during SA
        fp.sp.restore(best_state);
        fp.W = best_W;
        fp.H = best_H;
        fp.active_loc = best_loc;
        fp.pack();

        std::cerr << "[SA] Done. Best cost=" << best_cost
                  << " in " << elapsed() << "s\n";
    }

    // Estimate FT nets per soft block using bounding-box overlap heuristic
    void update_ft_areas() {
        // Reset
        std::fill(fp.ft_nets.begin(), fp.ft_nets.end(), 0);

        // For each connection, check if any soft block lies between source and destination
        // Simple bounding-box check: if soft block overlaps HPWL bounding box of connection
        for (auto& conn : fp.d.connections) {
            auto& ba = fp.d.blocks[conn.from];
            auto& bb = fp.d.blocks[conn.to];
            double min_x = std::min(ba.lx + ba.width/2, bb.lx + bb.width/2);
            double max_x = std::max(ba.lx + ba.width/2, bb.lx + bb.width/2);
            double min_y = std::min(ba.ly + ba.height/2, bb.ly + bb.height/2);
            double max_y = std::max(ba.ly + ba.height/2, bb.ly + bb.height/2);

            for (int i = 0; i < (int)fp.d.blocks.size(); i++) {
                if (i == conn.from || i == conn.to) continue;
                if (fp.d.blocks[i].type != BlockType::SOFT) continue;
                auto& bm = fp.d.blocks[i];
                double bx0 = bm.lx, bx1 = bm.lx + bm.width;
                double by0 = bm.ly, by1 = bm.ly + bm.height;
                // Check if block intersects the HPWL bounding box
                if (bx0 < max_x && bx1 > min_x && by0 < max_y && by1 > min_y) {
                    fp.ft_nets[i] += conn.nets;
                }
            }
        }
    }
};
