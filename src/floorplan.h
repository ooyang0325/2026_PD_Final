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

    Floorplan(Design& d_) : d(d_), sp((int)d_.blocks.size()),
        W(d_.blocks.size()), H(d_.blocks.size()),
        rotatable(d_.blocks.size(), false),
        ft_nets(d_.blocks.size(), 0),
        active_loc(d_.blocks.size(), 0) {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            W[i] = d.blocks[i].width;
            H[i] = d.blocks[i].height;
            rotatable[i] = (d.blocks[i].type == BlockType::SOFT);
        }
    }

    void flip_location(int i) {
        if (d.blocks[i].type != BlockType::EDGE) return;
        int nlocs = (int)d.blocks[i].locations.size();
        if (nlocs <= 1) return;
        active_loc[i] = (active_loc[i] + 1) % nlocs;
    }

    // Update soft block dimensions based on estimated feedthrough nets.
    // If reset_zero is true, restore base area when no FT nets are present.
    void apply_ft_areas(bool reset_zero = true) {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT) continue;
            double base = d.blocks[i].area;
            
            if (ft_nets[i] <= 0) {
                if (reset_zero) { 
                    // 保證基礎面積不會因進位而不足
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
            
            // 使用無條件進位至小數點後兩位，確保 W * H >= target_area
            W[i] = std::ceil(raw_w * 100.0) / 100.0;
            H[i] = std::ceil((target_area / W[i]) * 100.0) / 100.0;
        }
    }

    // Evaluate the sequence pair and write packed coordinates back to blocks.
    std::pair<double,double> pack() {
        auto [tw, th] = sp.evaluate(W, H);
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            d.blocks[i].lx = sp.x[i];
            d.blocks[i].ly = sp.y[i];
            d.blocks[i].width  = W[i];
            d.blocks[i].height = H[i];
        }
        // Leave boundary snapping to SA penalties to avoid biasing topology.
        return {tw, th};
    }

    double compute_hpwl() const {
        double total = 0;
        for (auto& conn : d.connections) {
            auto& ba = d.blocks[conn.from];
            auto& bb = d.blocks[conn.to];
            double cx_a = ba.lx + ba.width/2, cy_a = ba.ly + ba.height/2;
            double cx_b = bb.lx + bb.width/2, cy_b = bb.ly + bb.height/2;
            total += conn.nets * (std::abs(cx_b - cx_a) + std::abs(cy_b - cy_a));
        }
        return total;
    }

    // Area + alpha*HPWL with outline and edge-location penalties.
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

    // Enforce edge block location constraints via strong penalties.
    double edge_block_penalty(double chip_w, double chip_h) const {
        double penalty = 0;
        const double W_PENALTY = 1e6; // 極高懲罰值，強迫 SP 拓樸將 Edge 放邊界
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            auto& b = d.blocks[i];
            if (b.type != BlockType::EDGE || b.locations.empty()) continue;
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];
            char side = loc[0]; 
            if (side == 'T') {
                penalty += W_PENALTY * std::abs((b.ly + b.height) - chip_h);
            } else if (side == 'B') {
                penalty += W_PENALTY * std::abs(b.ly);
            } else if (side == 'L') {
                penalty += W_PENALTY * std::abs(b.lx);
            } else if (side == 'R') {
                penalty += W_PENALTY * std::abs((b.lx + b.width) - chip_w);
            }
        }
        return penalty;
    }
};