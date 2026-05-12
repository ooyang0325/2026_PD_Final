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

    void apply_ft_areas() {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT) continue;
            
            double target_area = d.blocks[i].get_target_area(ft_nets[i]);
            double cur_ar = (H[i] > 0) ? W[i] / H[i] : 1.0;
            double mn = d.blocks[i].min_ar, mx = d.blocks[i].max_ar;
            cur_ar = std::max(mn, std::min(mx, cur_ar));

            double raw_w = std::sqrt(target_area * cur_ar);
            // Force rounding up to two decimal places to prevent area shrinkage after output
            W[i] = std::ceil(raw_w * 100.0) / 100.0;
            H[i] = std::ceil((target_area / W[i]) * 100.0) / 100.0;
        }
    }

    std::pair<double,double> eval() {
        return sp.evaluate(W, H);
    }

    std::pair<double,double> pack() {
        auto [tw, th] = sp.evaluate(W, H);
        commit();
        // Do not add any forced modifications to lx, ly here!
        return {tw, th};
    }

    void commit() {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            d.blocks[i].lx = sp.x[i];
            d.blocks[i].ly = sp.y[i];
            d.blocks[i].width  = W[i];
            d.blocks[i].height = H[i];
        }
    }

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

        // Use a large penalty so SA naturally places edge blocks at the boundary
        cost += edge_block_penalty(chip_w, chip_h);
        return cost;
    }

    double edge_block_penalty(double chip_w, double chip_h) const {
        double penalty = 0;
        const double W_PENALTY = 1e6; // Very high displacement penalty
        for (int i : edge_block_idx) {
            const auto& b = d.blocks[i];
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];
            
            // If the SP-generated coordinates are not flush with the target edge, add a large penalty
            if (loc.find('T') != std::string::npos) penalty += W_PENALTY * std::abs((sp.y[i] + H[i]) - chip_h);
            if (loc.find('B') != std::string::npos) penalty += W_PENALTY * std::abs(sp.y[i]);
            if (loc.find('L') != std::string::npos) penalty += W_PENALTY * std::abs(sp.x[i]);
            if (loc.find('R') != std::string::npos) penalty += W_PENALTY * std::abs((sp.x[i] + W[i]) - chip_w);
        }
        return penalty;
    }
};