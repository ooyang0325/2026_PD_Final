#pragma once
#include "types.h"
#include "sequence_pair.h"
#include "channel.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

class Floorplan {
public:
    Design& d;
    SequencePair sp;

    std::vector<double> W, H;
    std::vector<bool> rotatable;
    std::vector<int> ft_nets;
    std::vector<int> active_loc;

    // Cached: indices of edge blocks (constant after construction)
    std::vector<int> edge_block_idx;

    Floorplan(Design& d_) : d(d_), sp((int)d_.blocks.size()),
        W(d_.blocks.size()), H(d_.blocks.size()),
        rotatable(d_.blocks.size(), false),
        ft_nets(d_.blocks.size(), 0),
        active_loc(d_.blocks.size(), 0) {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            W[i] = d.blocks[i].width;
            H[i] = d.blocks[i].height;
            rotatable[i] = (d.blocks[i].type == BlockType::SOFT);
            if (d.blocks[i].type == BlockType::EDGE && !d.blocks[i].locations.empty())
                edge_block_idx.push_back(i);
        }
    }

    void flip_location(int i) {
        if (d.blocks[i].type != BlockType::EDGE) return;
        int nlocs = (int)d.blocks[i].locations.size();
        if (nlocs <= 1) return;
        active_loc[i] = (active_loc[i] + 1) % nlocs;
    }

    // Update soft block dimensions based on estimated feedthrough nets.
    void apply_ft_areas(bool reset_zero = true) {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT) continue;
            double base = d.blocks[i].area;

            if (ft_nets[i] <= 0) {
                if (reset_zero) {
                    W[i] = std::ceil(std::sqrt(base) * 100.0) / 100.0;
                    H[i] = std::ceil((base / W[i]) * 100.0) / 100.0;
                }
                continue;
            }

            double rate = d.blocks[i].ft_conversion_rate(ft_nets[i]);
            double base_side = std::sqrt(base);
            double extend = ((double)ft_nets[i] / 25.0) * rate / 2.0;
            double target_area = (base_side + extend) * (base_side + extend);

            double cur_ar = (H[i] > 0) ? W[i] / H[i] : 1.0;
            double mn = d.blocks[i].min_ar, mx = d.blocks[i].max_ar;
            cur_ar = std::max(mn, std::min(mx, cur_ar));

            double raw_w = std::sqrt(target_area * cur_ar);
            W[i] = std::ceil(raw_w * 100.0) / 100.0;
            H[i] = std::ceil((target_area / W[i]) * 100.0) / 100.0;
        }
    }

    // Fast eval for SA inner loop: updates sp.x/sp.y, returns (tw, th).
    // Does NOT touch d.blocks (caller can read sp.x/sp.y/W/H directly).
    std::pair<double,double> eval() {
        return sp.evaluate(W, H);
    }

    // Eval + write back to d.blocks. Use outside SA hot loop (routing, output).
    std::pair<double,double> pack() {
        auto [tw, th] = sp.evaluate(W, H);
        commit();
        return {tw, th};
    }

    // Write sp.x/sp.y/W/H into d.blocks so external code (routing,
    // ChannelCalculator, output) sees current state.
    void commit() {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            d.blocks[i].lx = sp.x[i];
            d.blocks[i].ly = sp.y[i];
            d.blocks[i].width  = W[i];
            d.blocks[i].height = H[i];
        }
    }

    // HPWL on block centers, reading directly from sp.x/sp.y/W/H (not d.blocks).
    // This decouples cost eval from d.blocks staleness during SA.
    double compute_hpwl() const {
        double total = 0;
        const auto& spx = sp.x;
        const auto& spy = sp.y;
        for (auto& conn : d.connections) {
            int a = conn.from, b = conn.to;
            double cx_a = spx[a] + W[a] * 0.5;
            double cy_a = spy[a] + H[a] * 0.5;
            double cx_b = spx[b] + W[b] * 0.5;
            double cy_b = spy[b] + H[b] * 0.5;
            total += conn.nets * (std::abs(cx_b - cx_a) + std::abs(cy_b - cy_a));
        }
        return total;
    }

    double compute_cost(double chip_w, double chip_h, double alpha,
                        double max_w, double max_h) const {
        double area = chip_w * chip_h;
        double hpwl = compute_hpwl();
        double cost = area + alpha * hpwl;

        if (chip_w > max_w) cost += 1e8 * (chip_w - max_w);
        if (chip_h > max_h) cost += 1e8 * (chip_h - max_h);

        cost += edge_block_penalty(chip_w, chip_h);
        return cost;
    }

    // Edge block location penalty, reading from sp.x/sp.y/W/H.
    // Iterates only over cached edge block indices.
    double edge_block_penalty(double chip_w, double chip_h) const {
        double penalty = 0;
        const double W_PENALTY = 1e6;
        for (int i : edge_block_idx) {
            const auto& b = d.blocks[i];
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];
            char side = loc[0];
            double lx = sp.x[i], ly = sp.y[i];
            double w = W[i], h = H[i];
            if (side == 'T') {
                penalty += W_PENALTY * std::abs((ly + h) - chip_h);
            } else if (side == 'B') {
                penalty += W_PENALTY * std::abs(ly);
            } else if (side == 'L') {
                penalty += W_PENALTY * std::abs(lx);
            } else if (side == 'R') {
                penalty += W_PENALTY * std::abs((lx + w) - chip_w);
            }
        }
        return penalty;
    }
};
