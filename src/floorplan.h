#pragma once
#include "types.h"
#include "sequence_pair.h"
#include "channel.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <queue>

class Floorplan {
public:
    Design& d;
    SequencePair sp;

    std::vector<double> W, H;
    std::vector<bool> rotatable;
    std::vector<int> ft_nets;
    std::vector<int> active_loc;

    std::vector<int> edge_block_idx;
    std::vector<int> critical_conn_idx;
    double edge_penalty_coeff = 2e6;
    bool edge_use_best_location = true;

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
        critical_conn_idx.resize(d.connections.size());
        std::iota(critical_conn_idx.begin(), critical_conn_idx.end(), 0);
        std::sort(critical_conn_idx.begin(), critical_conn_idx.end(), [&](int a, int b) {
            return d.connections[a].nets > d.connections[b].nets;
        });
    }

    void flip_location(int i) {
        if (d.blocks[i].type != BlockType::EDGE) return;
        int nlocs = (int)d.blocks[i].locations.size();
        if (nlocs <= 1) return;
        active_loc[i] = (active_loc[i] + 1) % nlocs;
    }

    void choose_best_edge_locations(double chip_w, double chip_h) {
        for (int i : edge_block_idx) {
            const auto& b = d.blocks[i];
            if (b.locations.empty()) continue;
            int best_li = 0;
            double best_pen = 1e100;
            for (int li = 0; li < (int)b.locations.size(); li++) {
                double p = edge_option_penalty(i, b.locations[li], chip_w, chip_h);
                if (p < best_pen) {
                    best_pen = p;
                    best_li = li;
                }
            }
            active_loc[i] = best_li;
        }
    }

    // Snap EDGE blocks onto outline edges (evaluator requires exact alignment).
    void snap_edge_blocks(double chip_w, double chip_h) {
        if (edge_block_idx.empty() || chip_w < 1e-6 || chip_h < 1e-6) return;
        choose_best_edge_locations(chip_w, chip_h);
        for (int i : edge_block_idx) {
            const auto& b = d.blocks[i];
            if (b.locations.empty()) continue;
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];
            double x = sp.x[i];
            double y = sp.y[i];
            if (loc.find('L') != std::string::npos) x = 0.0;
            if (loc.find('R') != std::string::npos) x = chip_w - W[i];
            if (loc.find('B') != std::string::npos) y = 0.0;
            if (loc.find('T') != std::string::npos) y = chip_h - H[i];
            x = std::max(0.0, std::min(chip_w - W[i], x));
            y = std::max(0.0, std::min(chip_h - H[i], y));
            sp.x[i] = x;
            sp.y[i] = y;
        }
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
        double tw = 0.0, th = 0.0;
        // Outline can grow after snapping; re-snap to the latest extents.
        for (int iter = 0; iter < 3; iter++) {
            std::tie(tw, th) = sp.evaluate(W, H);
            snap_edge_blocks(tw, th);
        }
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
                        double max_w, double max_h,
                        double conn_weight = 0.0, int conn_top_k = 30,
                        double cong_weight = 0.0, int cong_bins = 10,
                        double hard_center_weight = 0.0) const {
        double area = chip_w * chip_h;
        double hpwl = compute_hpwl();
        double cost = area + alpha * hpwl;

        if (chip_w > max_w) cost += 1e8 * (chip_w - max_w);
        if (chip_h > max_h) cost += 1e8 * (chip_h - max_h);

        // Use a large penalty so SA naturally places edge blocks at the boundary
        cost += edge_block_penalty(chip_w, chip_h);
        if (conn_weight > 0.0) cost += conn_weight * connectivity_penalty(chip_w, chip_h, conn_top_k);
        if (cong_weight > 0.0) cost += cong_weight * congestion_penalty(chip_w, chip_h, cong_bins);
        if (hard_center_weight > 0.0) cost += hard_center_weight * hard_center_penalty(chip_w, chip_h);
        return cost;
    }

    // Penalize HARD_MACRO blocks sitting near the chip center (encourage periphery).
    double hard_center_penalty(double chip_w, double chip_h) const {
        if (chip_w < 1e-6 || chip_h < 1e-6) return 0.0;
        double cx = chip_w * 0.5, cy = chip_h * 0.5;
        double norm = chip_w * chip_w + chip_h * chip_h;
        double pen = 0.0;
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::HARD_MACRO) continue;
            double ccx = sp.x[i] + W[i] * 0.5;
            double ccy = sp.y[i] + H[i] * 0.5;
            double dx = (ccx - cx) / chip_w;
            double dy = (ccy - cy) / chip_h;
            pen += dx * dx + dy * dy;
        }
        return pen * norm;
    }

    double edge_block_penalty(double chip_w, double chip_h) const {
        double penalty = 0;
        for (int i : edge_block_idx) {
            const auto& b = d.blocks[i];
            if (b.locations.empty()) continue;
            double edge_dist = 0.0;
            if (edge_use_best_location) {
                edge_dist = 1e100;
                for (const auto& loc : b.locations) {
                    edge_dist = std::min(edge_dist, edge_option_penalty(i, loc, chip_w, chip_h));
                }
            } else {
                int li = std::min(active_loc[i], (int)b.locations.size() - 1);
                edge_dist = edge_option_penalty(i, b.locations[li], chip_w, chip_h);
            }
            penalty += edge_penalty_coeff * edge_dist;
        }
        return penalty;
    }

    double edge_option_penalty(int i, const std::string& loc, double chip_w, double chip_h) const {
        double p = 0.0;
        if (loc.find('T') != std::string::npos) p += std::abs((sp.y[i] + H[i]) - chip_h);
        if (loc.find('B') != std::string::npos) p += std::abs(sp.y[i]);
        if (loc.find('L') != std::string::npos) p += std::abs(sp.x[i]);
        if (loc.find('R') != std::string::npos) p += std::abs((sp.x[i] + W[i]) - chip_w);
        return p;
    }

    double connectivity_penalty(double chip_w, double chip_h, int top_k = 30) const {
        struct Rect {
            double lx, ly, w, h;
            bool traversable;
        };
        auto overlap = [](double a0, double a1, double b0, double b1) {
            return a0 < b1 - 1e-6 && b0 < a1 - 1e-6;
        };

        std::vector<Block> cur_blocks = d.blocks;
        for (int i = 0; i < (int)cur_blocks.size(); i++) {
            cur_blocks[i].lx = sp.x[i];
            cur_blocks[i].ly = sp.y[i];
            cur_blocks[i].width = W[i];
            cur_blocks[i].height = H[i];
        }
        auto channels = ChannelCalculator::compute(cur_blocks, chip_w, chip_h);

        std::vector<Rect> rects;
        rects.reserve(cur_blocks.size() + channels.size());
        for (int i = 0; i < (int)cur_blocks.size(); i++) {
            rects.push_back({cur_blocks[i].lx, cur_blocks[i].ly,
                             cur_blocks[i].width, cur_blocks[i].height,
                             cur_blocks[i].type == BlockType::SOFT});
        }
        for (auto& ch : channels) {
            rects.push_back({ch.lx, ch.ly, ch.width, ch.height, true});
        }

        int n = (int)rects.size();
        std::vector<std::vector<int>> adj(n);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                const auto& A = rects[i];
                const auto& B = rects[j];
                bool touch = false;
                if (std::abs((A.lx + A.w) - B.lx) < 1e-3 &&
                    overlap(A.ly, A.ly + A.h, B.ly, B.ly + B.h)) touch = true;
                if (std::abs((B.lx + B.w) - A.lx) < 1e-3 &&
                    overlap(A.ly, A.ly + A.h, B.ly, B.ly + B.h)) touch = true;
                if (std::abs((A.ly + A.h) - B.ly) < 1e-3 &&
                    overlap(A.lx, A.lx + A.w, B.lx, B.lx + B.w)) touch = true;
                if (std::abs((B.ly + B.h) - A.ly) < 1e-3 &&
                    overlap(A.lx, A.lx + A.w, B.lx, B.lx + B.w)) touch = true;
                if (!touch) continue;
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }

        auto is_reachable = [&](int src, int dst) {
            if (src == dst) return true;
            std::vector<char> vis(n, 0);
            std::queue<int> q;
            q.push(src);
            vis[src] = 1;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (int v : adj[u]) {
                    if (vis[v]) continue;
                    // Match router semantics: only source can leave non-traversable,
                    // only destination can be entered when non-traversable.
                    if (u != src && !rects[u].traversable) continue;
                    if (v != dst && !rects[v].traversable) continue;
                    if (v == dst) return true;
                    vis[v] = 1;
                    q.push(v);
                }
            }
            return false;
        };

        int use_k = top_k;
        if (use_k <= 0 || use_k > (int)critical_conn_idx.size())
            use_k = (int)critical_conn_idx.size();

        double penalty = 0.0;
        for (int rank = 0; rank < use_k; rank++) {
            const auto& conn = d.connections[critical_conn_idx[rank]];
            if (!is_reachable(conn.from, conn.to)) {
                penalty += std::pow((double)conn.nets, 1.2);
            }
        }
        return penalty;
    }

    double congestion_penalty(double chip_w, double chip_h, int bins = 10) const {
        if (bins < 2) bins = 2;
        double dx = chip_w / bins;
        double dy = chip_h / bins;
        if (dx < 1e-6 || dy < 1e-6) return 0.0;

        std::vector<double> demand_x((size_t)bins * bins, 0.0);
        std::vector<double> demand_y((size_t)bins * bins, 0.0);
        std::vector<double> blocked_area((size_t)bins * bins, 0.0);
        auto idx = [&](int ix, int iy) { return iy * bins + ix; };
        auto clampi = [&](int v) { return std::max(0, std::min(bins - 1, v)); };

        // Approximate bin blockage by overlapped block area.
        for (int bi = 0; bi < (int)d.blocks.size(); bi++) {
            double bx0 = sp.x[bi], bx1 = sp.x[bi] + W[bi];
            double by0 = sp.y[bi], by1 = sp.y[bi] + H[bi];
            int ix0 = clampi((int)std::floor(bx0 / dx));
            int ix1 = clampi((int)std::floor((bx1 - 1e-6) / dx));
            int iy0 = clampi((int)std::floor(by0 / dy));
            int iy1 = clampi((int)std::floor((by1 - 1e-6) / dy));
            for (int iy = iy0; iy <= iy1; iy++) {
                double cy0 = iy * dy, cy1 = (iy + 1) * dy;
                double oy = std::max(0.0, std::min(by1, cy1) - std::max(by0, cy0));
                if (oy <= 0.0) continue;
                for (int ix = ix0; ix <= ix1; ix++) {
                    double cx0 = ix * dx, cx1 = (ix + 1) * dx;
                    double ox = std::max(0.0, std::min(bx1, cx1) - std::max(bx0, cx0));
                    if (ox <= 0.0) continue;
                    blocked_area[idx(ix, iy)] += ox * oy;
                }
            }
        }

        const auto& spx = sp.x;
        const auto& spy = sp.y;
        for (auto& conn : d.connections) {
            int a = conn.from, b = conn.to;
            double ax = spx[a] + W[a] * 0.5, ay = spy[a] + H[a] * 0.5;
            double bx = spx[b] + W[b] * 0.5, by = spy[b] + H[b] * 0.5;
            double minx = std::min(ax, bx), maxx = std::max(ax, bx);
            double miny = std::min(ay, by), maxy = std::max(ay, by);
            double span_x = std::abs(bx - ax);
            double span_y = std::abs(by - ay);
            double span_sum = std::max(1e-6, span_x + span_y);
            // Direction-aware demand split (RUDY-like).
            double wx = 0.1 + 0.9 * (span_x / span_sum);
            double wy = 0.1 + 0.9 * (span_y / span_sum);
            int ix0 = clampi((int)std::floor(minx / dx));
            int ix1 = clampi((int)std::floor(maxx / dx));
            int iy0 = clampi((int)std::floor(miny / dy));
            int iy1 = clampi((int)std::floor(maxy / dy));
            int cells = std::max(1, (ix1 - ix0 + 1) * (iy1 - iy0 + 1));
            double unit = (double)conn.nets / cells;

            for (int iy = iy0; iy <= iy1; iy++) {
                for (int ix = ix0; ix <= ix1; ix++) {
                    demand_x[idx(ix, iy)] += unit * wx;
                    demand_y[idx(ix, iy)] += unit * wy;
                }
            }
        }

        double pen = 0.0;
        const double q = 1.8;
        const double bin_area = dx * dy;
        for (int iy = 0; iy < bins; iy++) {
            for (int ix = 0; ix < bins; ix++) {
                int id = idx(ix, iy);
                double free_ratio = 1.0 - blocked_area[id] / std::max(1e-6, bin_area);
                free_ratio = std::max(0.05, std::min(1.0, free_ratio));
                double cap_x = dy * 25.0 * free_ratio;
                double cap_y = dx * 25.0 * free_ratio;
                double ox = std::max(0.0, demand_x[id] / std::max(1e-6, cap_x) - 1.0);
                double oy = std::max(0.0, demand_y[id] / std::max(1e-6, cap_y) - 1.0);
                pen += std::pow(ox, q) + std::pow(oy, q);
            }
        }
        return pen;
    }
};