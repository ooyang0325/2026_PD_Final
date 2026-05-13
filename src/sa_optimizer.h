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
    std::mt19937 rng;

    double T_init = 1e8;
    double T_final = 10.0;
    double cool_rate = 0.99;
    int moves_per_temp = 300; // 增加探索次數以應付 50 blocks
    double time_limit_sec = 6000.0;

    SAOptimizer(Floorplan& fp_, unsigned seed = 42) : fp(fp_), rng(seed) {}

    void run() {
        auto t0 = std::chrono::steady_clock::now();
        auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

        int n = fp.sp.n;
        double alpha = fp.d.alpha;
        double max_w = fp.d.outline.max_width;
        double max_h = fp.d.outline.max_height;

        fp.apply_ft_areas();
        auto [tw0, th0] = fp.eval();
        double cur_cost = fp.compute_cost(tw0, th0, alpha, max_w, max_h);

        SequencePair::State best_state = fp.sp.save();
        std::vector<double> best_W = fp.W, best_H = fp.H;
        std::vector<int>    best_loc = fp.active_loc;
        double best_cost = cur_cost;

        int iter = 0;
        int phase = 0;
        double T = T_init;

        std::uniform_int_distribution<int> move_type(1, 6);
        std::uniform_int_distribution<int> bidx(0, n-1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        while (elapsed() < time_limit_sec) {
            if (phase == 0 && T < T_init * 1e-4) {
                phase = 1;
            }

            // 動態 FT 面積反饋：每 500 次迭代，重新評估 FT 需求，提前擴大 Soft block
            if (phase == 1 && iter % 500 == 0) {
                update_ft_areas();
                fp.apply_ft_areas(false);
            }

            for (int m = 0; m < moves_per_temp; m++) {
                if (elapsed() >= time_limit_sec) break;

                int mv = move_type(rng);
                int bi = bidx(rng);
                int bj = (mv <= 3) ? bidx(rng) : 0;
                if (mv <= 3) { while (bj == bi) bj = bidx(rng); }

                double saved_W = 0, saved_H = 0;
                int saved_loc = 0;
                bool changed = false;

                switch (mv) {
                    case 1: fp.sp.swap_gp(bi, bj); changed = true; break;
                    case 2: fp.sp.swap_gm(bi, bj); changed = true; break;
                    case 3: fp.sp.swap_both(bi, bj); changed = true; break;
                    case 4: {
                        if (fp.d.blocks[bi].type == BlockType::SOFT) {
                            double target_area = fp.d.blocks[bi].get_target_area(fp.ft_nets[bi]);
                            double new_w = fp.H[bi], new_h = fp.W[bi];
                            double ar = (new_h > 0) ? new_w / new_h : 1.0;
                            if (ar >= fp.d.blocks[bi].min_ar - 1e-6 && ar <= fp.d.blocks[bi].max_ar + 1e-6) {
                                saved_W = fp.W[bi]; saved_H = fp.H[bi];
                                fp.W[bi] = std::ceil(new_w * 100.0) / 100.0;
                                fp.H[bi] = std::ceil((target_area / fp.W[bi]) * 100.0) / 100.0;
                                changed = true;
                            }
                        }
                        break;
                    }
                    case 5: {
                        if (fp.d.blocks[bi].type == BlockType::SOFT) {
                            double target_area = fp.d.blocks[bi].get_target_area(fp.ft_nets[bi]);
                            std::uniform_real_distribution<double> ar_dist(fp.d.blocks[bi].min_ar, fp.d.blocks[bi].max_ar);
                            double new_ar = ar_dist(rng);
                            saved_W = fp.W[bi]; saved_H = fp.H[bi];
                            double raw_w = std::sqrt(target_area * new_ar);
                            fp.W[bi] = std::ceil(raw_w * 100.0) / 100.0;
                            fp.H[bi] = std::ceil((target_area / fp.W[bi]) * 100.0) / 100.0;
                            changed = true;
                        }
                        break;
                    }
                    case 6: {
                        if (fp.d.blocks[bi].type == BlockType::EDGE && fp.d.blocks[bi].locations.size() > 1) {
                            saved_loc = fp.active_loc[bi];
                            fp.flip_location(bi);
                            changed = true;
                        }
                        break;
                    }
                }

                if (!changed) continue;

                auto [ntw, nth] = fp.eval();
                double new_cost = fp.compute_cost(ntw, nth, alpha, max_w, max_h);
                double delta = new_cost - cur_cost;

                if (delta <= 0 || unit(rng) < std::exp(-delta / T)) {
                    cur_cost = new_cost;
                    if (new_cost < best_cost) {
                        best_cost = new_cost;
                        best_state = fp.sp.save();
                        best_W = fp.W; best_H = fp.H;
                        best_loc = fp.active_loc;
                    }
                } else {
                    switch (mv) {
                        case 1: fp.sp.swap_gp(bi, bj); break;
                        case 2: fp.sp.swap_gm(bi, bj); break;
                        case 3: fp.sp.swap_both(bi, bj); break;
                        case 4: case 5: fp.W[bi] = saved_W; fp.H[bi] = saved_H; break;
                        case 6: fp.active_loc[bi] = saved_loc; break;
                    }
                }
            }

            double t_frac = std::min(1.0, elapsed() / time_limit_sec);
            T = T_init * std::pow(T_final / T_init, t_frac);
            iter++;
            
            if (iter % 50000 == 0) {
                std::cerr << "[SA] iter=" << iter << " T=" << T
                          << " best=" << best_cost
                          << " time=" << elapsed() << "s\n";
            }
        }

        fp.sp.restore(best_state);
        fp.W = best_W; fp.H = best_H; fp.active_loc = best_loc;
        fp.pack();
    }

    void update_ft_areas() {
        std::fill(fp.ft_nets.begin(), fp.ft_nets.end(), 0);
        const auto& spx = fp.sp.x; const auto& spy = fp.sp.y;
        const auto& W = fp.W; const auto& H = fp.H;

        for (auto& conn : fp.d.connections) {
            int a = conn.from, b = conn.to;
            double cax = spx[a] + W[a] * 0.5, cay = spy[a] + H[a] * 0.5;
            double cbx = spx[b] + W[b] * 0.5, cby = spy[b] + H[b] * 0.5;
            double min_x = std::min(cax, cbx), max_x = std::max(cax, cbx);
            double min_y = std::min(cay, cby), max_y = std::max(cay, cby);

            for (int i = 0; i < (int)fp.d.blocks.size(); i++) {
                if (i == a || i == b || fp.d.blocks[i].type != BlockType::SOFT) continue;
                if (spx[i] < max_x && spx[i] + W[i] > min_x && spy[i] < max_y && spy[i] + H[i] > min_y) {
                    fp.ft_nets[i] += conn.nets;
                }
            }
        }
    }
};