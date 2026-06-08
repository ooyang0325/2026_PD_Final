#pragma once
#include "floorplan.h"
#include "ft_estimator.h"
#include "config.h"
#include "channel.h"
#include <random>
#include <cmath>
#include <iostream>
#include <chrono>

// Simulated annealing over the B*-Tree representation.
// Moves:
//   1  swap two tree nodes' blocks
//   2  move a block to another tree position (delete-to-leaf + reinsert)
//   3  rotate a SOFT block (swap W/H within aspect-ratio limits)
//   4  change a SOFT block's aspect ratio
//   5  flip an EDGE block's allowed location
// Feedthrough area pre-sizing is driven by the FTAFP estimation model
// (ftest::estimate) instead of routing inside the loop.
class SAOptimizer {
public:
    Floorplan& fp;
    std::mt19937 rng;

    double T_init = 1e8;
    double T_final = 10.0;
    int moves_per_temp = 300;
    double time_limit_sec = 6000.0;

    SAOptimizer(Floorplan& fp_, unsigned seed = 42) : fp(fp_), rng(seed) {}

    void run() {
        auto t0 = std::chrono::steady_clock::now();
        auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

        int n = (int)fp.d.blocks.size();
        double alpha = fp.d.alpha;
        double max_w = fp.d.outline.max_width;
        double max_h = fp.d.outline.max_height;

        fp.apply_ft_areas();

        std::uniform_int_distribution<int> move_type(1, 5);
        std::uniform_int_distribution<int> bidx(0, n - 1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        // ── Normalization sampling (PA2-style) ───────────────────────────────
        // Sample random tree perturbations to estimate average area / HPWL so the
        // cost terms are O(1), and to derive an initial temperature from the
        // average uphill cost step.  Keep the most compact sample as the start.
        BStarTree::State move_snap;
        {
            auto [tw0, th0] = fp.eval();
            int samples = std::max(40, 3 * n);
            std::vector<double> As, Ws; As.reserve(samples); Ws.reserve(samples);
            std::vector<double> TWs, THs; TWs.reserve(samples); THs.reserve(samples);
            double a0 = tw0 * th0, w0 = fp.compute_hpwl();
            double f0 = ftest::cost(fp.d.blocks, fp.d.connections, fp.bst.x, fp.bst.y, fp.W, fp.H);
            As.push_back(a0); Ws.push_back(w0); TWs.push_back(tw0); THs.push_back(th0);
            double sumA = a0, sumW = w0, sumF = f0;

            BStarTree::State best_start = fp.bst.save();
            double best_metric = a0 + alpha * w0;

            for (int s = 1; s < samples; s++) {
                if (unit(rng) < 0.5) {
                    int i = bidx(rng), j = bidx(rng);
                    while (j == i) j = bidx(rng);
                    fp.bst.swap_nodes(i, j);
                } else {
                    fp.bst.move_random(rng);
                }
                auto [t1, t2] = fp.eval();
                double a = t1 * t2, w = fp.compute_hpwl();
                double f = ftest::cost(fp.d.blocks, fp.d.connections, fp.bst.x, fp.bst.y, fp.W, fp.H);
                As.push_back(a); Ws.push_back(w); TWs.push_back(t1); THs.push_back(t2);
                sumA += a; sumW += w; sumF += f;
                double metric = a + alpha * w;
                if (metric < best_metric) { best_metric = metric; best_start = fp.bst.save(); }
            }

            fp.Anorm = std::max(1.0, sumA / samples);
            fp.Wnorm = std::max(1.0, sumW / samples);
            fp.Fnorm = std::max(1.0, sumF / samples);
            fp.ftw   = cfg::FTW; // feedthrough-minimizing penalty: pack connected blocks adjacent

            // Average uphill step of the FULL normalized cost INCLUDING the outline
            // violation term (gamma).  Calibrating the temperature to area+HPWL only
            // (as a naive normalization would) makes T ~1000x too cold for the real
            // landscape, where the gamma violation dominates — the SA then freezes in
            // the initial aspect ratio and never reshapes to fit the outline.
            auto full_cost = [&](int i) {
                double base = As[i] / fp.Anorm + alpha * Ws[i] / fp.Wnorm;
                double wv = std::max(0.0, TWs[i] - max_w) / max_w;
                double hv = std::max(0.0, THs[i] - max_h) / max_h;
                return base + fp.gamma * (wv + hv);
            };
            double dsum = 0; int dc = 0;
            double cp = full_cost(0);
            for (int i = 1; i < samples; i++) {
                double cc = full_cost(i);
                if (cc > cp) { dsum += cc - cp; dc++; }
                cp = cc;
            }
            double delta_avg = dc ? dsum / dc : 1e-3;
            delta_avg = std::max(1e-6, delta_avg);

            // accept ~85% of uphill moves initially; anneal to ~1e-4 of that
            T_init  = -delta_avg / std::log(0.85);
            T_final = std::max(1e-9, T_init * 1e-4);

            fp.bst.restore(best_start); // start from the most compact sample
        }

        auto [tw0, th0] = fp.eval();
        double cur_cost = fp.compute_cost(tw0, th0, alpha, max_w, max_h);

        // Best-state checkpoint
        BStarTree::State best_tree = fp.bst.save();
        std::vector<double> best_W = fp.W, best_H = fp.H;
        std::vector<int>    best_loc = fp.active_loc;
        double best_cost = cur_cost;

        int iter = 0;
        double T = T_init;

        while (elapsed() < time_limit_sec) {
            // Soft blocks stay at BASE area during the SA so a fitting layout
            // always exists (the validity anchor).  Feedthrough is handled two
            // ways instead: (1) the FTAFP estimate is folded into the cost as a
            // penalty (fp.ftw) so connected blocks pack adjacent and real routing
            // feedthrough stays low; (2) the robust convergence loop in main.cpp
            // grows soft blocks to the real feedthrough afterwards, but only ever
            // keeps layouts that remain inside the outline.

            for (int m = 0; m < moves_per_temp; m++) {
                if (elapsed() >= time_limit_sec) break;

                int mv = move_type(rng);
                int bi = bidx(rng);

                int saved_a = 0, saved_b = 0;     // for swap undo
                double saved_W = 0, saved_H = 0;  // for resize undo
                int saved_loc = 0;                // for edge flip undo
                bool changed = false;

                switch (mv) {
                    case 1: { // swap two nodes' blocks
                        int bj = bidx(rng);
                        while (bj == bi) bj = bidx(rng);
                        saved_a = bi; saved_b = bj;
                        fp.bst.swap_nodes(bi, bj);
                        changed = true;
                        break;
                    }
                    case 2: { // relocate a block in the tree
                        fp.bst.save(move_snap);
                        fp.bst.move_random(rng);
                        changed = true;
                        break;
                    }
                    case 3: { // rotate a SOFT block
                        if (fp.d.blocks[bi].type == BlockType::SOFT) {
                            double new_w = fp.H[bi], new_h = fp.W[bi];
                            double ar = (new_h > 0) ? new_w / new_h : 1.0;
                            if (ar >= fp.d.blocks[bi].min_ar - 1e-6 && ar <= fp.d.blocks[bi].max_ar + 1e-6) {
                                saved_W = fp.W[bi]; saved_H = fp.H[bi];
                                fp.W[bi] = new_w; fp.H[bi] = new_h;
                                changed = true;
                            }
                        }
                        break;
                    }
                    case 4: { // change a SOFT block's aspect ratio
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
                    case 5: { // flip an EDGE block location
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
                        best_tree = fp.bst.save();
                        best_W = fp.W; best_H = fp.H;
                        best_loc = fp.active_loc;
                    }
                } else {
                    switch (mv) {
                        case 1: fp.bst.swap_nodes(saved_a, saved_b); break;
                        case 2: fp.bst.restore(move_snap); break;
                        case 3: case 4: fp.W[bi] = saved_W; fp.H[bi] = saved_H; break;
                        case 5: fp.active_loc[bi] = saved_loc; break;
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

        fp.bst.restore(best_tree);
        fp.W = best_W; fp.H = best_H; fp.active_loc = best_loc;
        fp.pack();
    }
};
