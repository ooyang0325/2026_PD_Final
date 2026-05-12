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

    // Cache the indices of edge blocks to avoid repeated checks inside the SA loop
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

    // High-frequency call inside SA: only evaluate positions without writing back to d.blocks
    std::pair<double,double> eval() {
        return sp.evaluate(W, H);
    }

    // Used for final packing: evaluate, force edge blocks to snap to the boundary, and write coordinates back to d.blocks
    std::pair<double,double> pack() {
        auto [tw, th] = sp.evaluate(W, H);
        
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            d.blocks[i].lx = sp.x[i];
            d.blocks[i].ly = sp.y[i];
            d.blocks[i].width  = W[i];
            d.blocks[i].height = H[i];
        }

        // Force edge blocks to the boundary
        for (int i : edge_block_idx) {
            auto& b = d.blocks[i];
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];
            
            if (loc.find('T') != std::string::npos) b.ly = th - b.height;
            if (loc.find('B') != std::string::npos) b.ly = 0.0;
            if (loc.find('L') != std::string::npos) b.lx = 0.0;
            if (loc.find('R') != std::string::npos) b.lx = tw - b.width;
        }

        return {tw, th};
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

        // Add boundary and overlap penalties
        cost += edge_block_penalty(chip_w, chip_h);
        return cost;
    }

    // Virtual snapping and collision testing (Virtual Snapping & Overlap Penalty)
    double edge_block_penalty(double chip_w, double chip_h) const {
        double penalty = 0;
        const double OVERLAP_W = 1e6; // Very high penalty for overlap
        const double DIST_W = 1e1;    // Mild distance penalty to guide the SP toward the target

        int n = d.blocks.size();
        
        // Use static arrays to avoid memory allocation inside the SA hot loop (problem guarantees BLOCK < 50)
        double sx[256];
        double sy[256];
        for (int i = 0; i < n && i < 256; i++) {
            sx[i] = sp.x[i];
            sy[i] = sp.y[i];
        }

        // 1. Precompute the virtual snapped coordinates for each edge block
        for (int i : edge_block_idx) {
            int li = std::min(active_loc[i], (int)d.blocks[i].locations.size() - 1);
            const std::string& loc = d.blocks[i].locations[li];
            
            // Strictly test each character so TL satisfies both T and L
            if (loc.find('T') != std::string::npos) sy[i] = chip_h - H[i];
            if (loc.find('B') != std::string::npos) sy[i] = 0.0;
            if (loc.find('L') != std::string::npos) sx[i] = 0.0;
            if (loc.find('R') != std::string::npos) sx[i] = chip_w - W[i];

            // Add a small pulling penalty to encourage SA to place it near the target position and keep the topology reasonable
            penalty += DIST_W * (std::abs(sx[i] - sp.x[i]) + std::abs(sy[i] - sp.y[i]));
        }

        // 2. Perform O(N^2) collision detection on the virtual snapped coordinates
        for (int i : edge_block_idx) {
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                // Prevent double counting between two edge blocks
                if (d.blocks[j].type == BlockType::EDGE && i >= j) continue;

                double ox = std::min(sx[i] + W[i], sx[j] + W[j]) - std::max(sx[i], sx[j]);
                double oy = std::min(sy[i] + H[i], sy[j] + H[j]) - std::max(sy[i], sy[j]);
                
                // If overlap occurs after virtual snapping, give a destructive penalty
                if (ox > 1e-6 && oy > 1e-6) {
                    penalty += OVERLAP_W * (ox * oy);
                }
            }
        }
        return penalty;
    }
};