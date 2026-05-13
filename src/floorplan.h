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

    // 評估時強制加入 2.0um 的 Halo，保證 Block 之間絕對會產生可繞線的 Channel
    std::pair<double,double> eval() {
        const double HALO = 2.0;
        std::vector<double> pad_W = W, pad_H = H;
        for(size_t i=0; i<W.size(); i++) { pad_W[i] += HALO; pad_H[i] += HALO; }
        return sp.evaluate(pad_W, pad_H);
    }

    std::pair<double,double> pack() {
        auto [tw, th] = eval();
        commit();
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

    // 最後定案時，將 Edge Block 強制貼齊大會設定的 MAX OUTLINE
    void finalize_edge_blocks(double max_w, double max_h) {
        for (int i : edge_block_idx) {
            auto& b = d.blocks[i];
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];
            if (loc.find('T') != std::string::npos) b.ly = max_h - b.height;
            if (loc.find('B') != std::string::npos) b.ly = 0.0;
            if (loc.find('L') != std::string::npos) b.lx = 0.0;
            if (loc.find('R') != std::string::npos) b.lx = max_w - b.width;
        }
    }

    // Push edge blocks to extreme positions in the sequence pair so that
    // finalize_edge_blocks() never creates block overlaps.  The rules are:
    //   'B' (bottom snap) → first in gp, last in gm  → block is packed at y=0
    //   'T' (top snap)    → last  in gp, first in gm → block is packed at y=max
    //   'L' only          → first in both gp and gm  → block is packed at x=0
    //   'R' only          → last  in both gp and gm  → block is packed at x=max
    // After the snap, all non-edge blocks are strictly inside the snapped
    // edge block's zone (guaranteed by the HALO gap in the sequence-pair packing).
    void enforce_edge_block_extremes() {
        auto move_to_front = [](std::vector<int>& v, int val) {
            auto it = std::find(v.begin(), v.end(), val);
            if (it != v.begin()) std::rotate(v.begin(), it, it + 1);
        };
        auto move_to_back = [](std::vector<int>& v, int val) {
            auto it = std::find(v.begin(), v.end(), val);
            if (it != v.end() - 1) std::rotate(it, it + 1, v.end());
        };

        for (int i : edge_block_idx) {
            int li = std::min(active_loc[i], (int)d.blocks[i].locations.size() - 1);
            const std::string& loc = d.blocks[i].locations[li];
            bool has_B = (loc.find('B') != std::string::npos);
            bool has_T = (loc.find('T') != std::string::npos);
            bool has_L = (loc.find('L') != std::string::npos);
            bool has_R = (loc.find('R') != std::string::npos);

            if (has_B) {
                move_to_front(sp.gp, i);
                move_to_back (sp.gm, i);
            } else if (has_T) {
                move_to_back (sp.gp, i);
                move_to_front(sp.gm, i);
            } else if (has_L) {
                move_to_front(sp.gp, i);
                move_to_front(sp.gm, i);
            } else if (has_R) {
                move_to_back (sp.gp, i);
                move_to_back (sp.gm, i);
            }
        }
    }

    double compute_hpwl() const {
        double total = 0;
        for (auto& conn : d.connections) {
            int a = conn.from, b = conn.to;
            double cx_a = sp.x[a] + W[a] * 0.5, cy_a = sp.y[a] + H[a] * 0.5;
            double cx_b = sp.x[b] + W[b] * 0.5, cy_b = sp.y[b] + H[b] * 0.5;
            total += conn.nets * (std::abs(cx_b - cx_a) + std::abs(cy_b - cy_a));
        }
        return total;
    }

    double compute_cost(double chip_w, double chip_h, double alpha, double max_w, double max_h) const {
        double cost = (chip_w * chip_h) + alpha * compute_hpwl();
        if (chip_w > max_w) cost += 1e8 * (chip_w - max_w);
        if (chip_h > max_h) cost += 1e8 * (chip_h - max_h);
        cost += edge_block_penalty(chip_w, chip_h); // push edge blocks to actual chip boundary
        return cost;
    }

    double edge_block_penalty(double max_w, double max_h) const {
        double penalty = 0;
        const double DIST_W = 50.0;
        const double OVERLAP_W = 1e6;

        int n = d.blocks.size();
        std::vector<double> sx(n), sy(n);
        for (int i = 0; i < n; i++) { sx[i] = sp.x[i]; sy[i] = sp.y[i]; }

        // 預判將 Edge Block 推向 MAX 邊緣
        for (int i : edge_block_idx) {
            int li = std::min(active_loc[i], (int)d.blocks[i].locations.size() - 1);
            const std::string& loc = d.blocks[i].locations[li];
            
            if (loc.find('T') != std::string::npos) sy[i] = max_h - H[i];
            if (loc.find('B') != std::string::npos) sy[i] = 0.0;
            if (loc.find('L') != std::string::npos) sx[i] = 0.0;
            if (loc.find('R') != std::string::npos) sx[i] = max_w - W[i];

            penalty += DIST_W * (std::abs(sx[i] - sp.x[i]) + std::abs(sy[i] - sp.y[i]));
        }

        // 碰撞測試：確保推到邊緣後，不會跟其他 Block 撞在一起
        for (int i : edge_block_idx) {
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                if (d.blocks[j].type == BlockType::EDGE && i >= j) continue;
                double ox = std::min(sx[i] + W[i], sx[j] + W[j]) - std::max(sx[i], sx[j]);
                double oy = std::min(sy[i] + H[i], sy[j] + H[j]) - std::max(sy[i], sy[j]);
                if (ox > 1e-6 && oy > 1e-6) penalty += OVERLAP_W * (ox * oy);
            }
        }
        return penalty;
    }
};