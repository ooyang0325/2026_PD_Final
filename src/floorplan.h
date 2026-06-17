#pragma once
#include "types.h"
#include "bstree.h"
#include "channel.h"
#include "ft_estimator.h"
#include "config.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// Floorplan built on a B*-Tree representation (replaces Sequence Pair).
// Public surface kept identical to the previous version so main.cpp and the
// output/channel pipeline are unaffected:
//   W, H, ft_nets, active_loc, edge_block_idx
//   pack() / eval() / commit()
//   apply_ft_areas(bool) / finalize_edge_blocks(max_w,max_h)
//   compute_cost(...) / compute_hpwl() / edge_block_penalty(...)
//   flip_location(i)

class Floorplan {
public:
    Design& d;
    BStarTree bst;

    std::vector<double> W, H;
    std::vector<bool> rotatable;
    std::vector<int> ft_nets;
    std::vector<int> active_loc;
    std::vector<int> edge_block_idx;

    // Cost normalization (set by the SA from sampling).  A well-conditioned
    // normalized cost is what makes the fixed-outline constraint actually
    // converge for large block counts (cf. the PA2 B*-tree floorplanner).
    double Anorm = 1.0;   // average sampled chip area
    double Wnorm = 1.0;   // average sampled HPWL
    double Fnorm = 1.0;   // average sampled feedthrough estimate
    double gamma = 100.0; // outline-violation weight (PA2 uses 100)
    double ftw   = 0.0;   // feedthrough penalty weight (0 until SA enables it)

    // padded scratch for HALO evaluation (reused, no per-eval alloc)
    std::vector<double> pad_W, pad_H;

    Floorplan(Design& d_) : d(d_), bst((int)d_.blocks.size()),
        W(d_.blocks.size()), H(d_.blocks.size()),
        rotatable(d_.blocks.size(), false),
        ft_nets(d_.blocks.size(), 0),
        active_loc(d_.blocks.size(), 0),
        pad_W(d_.blocks.size()), pad_H(d_.blocks.size()) {
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

    // Resize SOFT blocks to satisfy the (estimated) feedthrough area requirement.
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

    // Evaluate packing with a 2.0um HALO around every block so a routable
    // channel is always guaranteed between neighbors.  Does NOT commit.
    std::pair<double,double> eval() {
        const double HALO = cfg::HALO;
        for (size_t i = 0; i < W.size(); i++) { pad_W[i] = W[i] + HALO; pad_H[i] = H[i] + HALO; }
        return bst.pack(pad_W, pad_H);
    }

    std::pair<double,double> pack() {
        auto r = eval();
        commit();
        return r;
    }

    void commit() {
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            d.blocks[i].lx = bst.x[i];
            d.blocks[i].ly = bst.y[i];
            d.blocks[i].width  = W[i];
            d.blocks[i].height = H[i];
        }
    }

    // Snap edge blocks to the requested boundary (compact chip outline),
    // overlap-safely: a block is moved to its boundary only if the snapped
    // position does not collide with any other block.  Otherwise it keeps its
    // packed position.  This prevents snap-induced overlap FAILs when the edge
    // constraints are over-subscribed (e.g. two blocks both requiring the same
    // corner) — at worst one edge constraint is left unsatisfied rather than
    // producing a cascade of overlaps.
    void finalize_edge_blocks(double max_w, double max_h) {
        int nb = (int)d.blocks.size();
        for (int i : edge_block_idx) {
            auto& b = d.blocks[i];
            if (b.locations.empty()) continue;

            // LOCATION codes use AND semantics (QA A8): 1st letter = flush boundary
            // (T/B/L/R), 2nd letter = which equal-third along it (QA A10). Accumulate
            // the UNION of required thirds per axis so a block carrying several
            // same-edge codes (e.g. TL,TM) overlaps every one (QA A19/A20).
            bool fT=false, fB=false, fL=false, fR=false;
            double xlo=1e18, xhi=-1e18, ylo=1e18, yhi=-1e18;
            auto third_iv = [](char c, double L, bool xaxis, double& lo, double& hi) {
                double t3 = L / 3.0;
                if (xaxis) {
                    if      (c=='L') { lo=0;    hi=t3;   }
                    else if (c=='R') { lo=2*t3; hi=L;    }
                    else             { lo=t3;   hi=2*t3; }
                } else {
                    if      (c=='T') { lo=2*t3; hi=L;    }
                    else if (c=='B') { lo=0;    hi=t3;   }
                    else             { lo=t3;   hi=2*t3; }
                }
            };
            for (const auto& loc : b.locations) {
                if (loc.empty()) continue;
                char e = (char)std::toupper((unsigned char)loc[0]);
                char t = (loc.size() >= 2) ? (char)std::toupper((unsigned char)loc[1]) : 'M';
                double lo, hi;
                if (e=='T' || e=='B') { (e=='T'?fT:fB)=true; third_iv(t,max_w,true, lo,hi); xlo=std::min(xlo,lo); xhi=std::max(xhi,hi); }
                else if (e=='L' || e=='R') { (e=='L'?fL:fR)=true; third_iv(t,max_h,false, lo,hi); ylo=std::min(ylo,lo); yhi=std::max(yhi,hi); }
            }

            double ox=b.lx, oy=b.ly, nx=b.lx, ny=b.ly;
            if      (fT) ny = max_h - b.height;
            else if (fB) ny = 0.0;
            if      (fL) nx = 0.0;
            else if (fR) nx = max_w - b.width;
            if (!fL && !fR && xhi > xlo)
                nx = std::max(0.0, std::min(max_w - b.width,  (xlo+xhi)*0.5 - b.width *0.5));
            if (!fT && !fB && yhi > ylo)
                ny = std::max(0.0, std::min(max_h - b.height, (ylo+yhi)*0.5 - b.height*0.5));

            b.lx = nx; b.ly = ny;
            bool collide = false;
            for (int j = 0; j < nb && !collide; j++) {
                if (j == i) continue;
                const auto& o = d.blocks[j];
                double cx = std::min(b.lx + b.width,  o.lx + o.width)  - std::max(b.lx, o.lx);
                double cy = std::min(b.ly + b.height, o.ly + o.height) - std::max(b.ly, o.ly);
                if (cx > 0.01 && cy > 0.01) collide = true;
            }
            if (collide) { b.lx = ox; b.ly = oy; } // keep packed position
        }
    }

    // Pin = port-edge midpoint when the block declares a port edge, else center.
    std::pair<double,double> pin_xy(int i) const {
        double x = bst.x[i], y = bst.y[i], w = W[i], h = H[i];
        switch (d.blocks[i].port_edge) {
            case 1: return {x,           y + h * 0.5};
            case 2: return {x + w * 0.5, y + h};
            case 3: return {x + w,       y + h * 0.5};
            case 4: return {x + w * 0.5, y};
            default: return {x + w * 0.5, y + h * 0.5};
        }
    }

    double compute_hpwl() const {
        double total = 0;
        for (auto& conn : d.connections) {
            auto [ax, ay] = pin_xy(conn.from);
            auto [bx, by] = pin_xy(conn.to);
            total += conn.nets * (std::abs(bx - ax) + std::abs(by - ay));
        }
        return total;
    }

    // Routing-aware placement term: only HIGH-demand connections (those liable to
    // overflow a single channel) contribute, weighted by how much they exceed the
    // threshold times their separation.  Minimizing this pulls heavily-connected
    // blocks adjacent so their connecting channel is wide (high capacity) and the
    // net does not have to feed through / congest narrow channels — cutting both
    // channel and feedthrough overflow at the source.  O(#connections).
    double compute_congestion() const {
        double total = 0;
        for (auto& conn : d.connections) {
            double excess = conn.nets - cfg::PCONG_THRESH;
            if (excess <= 0) continue;
            int a = conn.from, b = conn.to;
            double cx_a = bst.x[a] + W[a] * 0.5, cy_a = bst.y[a] + H[a] * 0.5;
            double cx_b = bst.x[b] + W[b] * 0.5, cy_b = bst.y[b] + H[b] * 0.5;
            total += excess * (std::abs(cx_b - cx_a) + std::abs(cy_b - cy_a));
        }
        return total;
    }

    // Normalized SA cost (PA2-style): area and HPWL are scaled by their sampled
    // averages so all terms are O(1); the fixed outline is enforced by a strong
    // normalized violation penalty (gamma).  The true contest metric
    // (area + alpha*HPWL) is still used for final solution ranking in main.
    double compute_cost(double chip_w, double chip_h, double alpha, double max_w, double max_h) const {
        double area = chip_w * chip_h;
        double base = area / Anorm + alpha * (compute_hpwl() / Wnorm);

        double wv = std::max(0.0, chip_w - max_w) / max_w;
        double hv = std::max(0.0, chip_h - max_h) / max_h;
        double penalty = gamma * (wv + hv);

        if (ftw > 0.0) {
            double ft = ftest::cost(d.blocks, d.connections, bst.x, bst.y, W, H);
            penalty += ftw * (ft / Fnorm);
        }

        // Routing-aware: keep high-demand pairs close (wide connecting channel).
        if (cfg::PCONGW > 0.0)
            penalty += cfg::PCONGW * (compute_congestion() / Wnorm);

        penalty += edge_block_penalty(max_w, max_h);
        return base + penalty;
    }

    // Penalise edge blocks for not being at the chip boundary, plus a collision
    // term for the projected snap (so the SA avoids snaps that would overlap).
    // Normalized to the same O(1) scale as the base cost.
    double edge_block_penalty(double max_w, double max_h) const {
        double penalty = 0;
        const double W_DIST = 5.0;     // mild pull toward the boundary
        const double W_OVL  = 200.0;   // projected snap overlap is a hard FAIL -> heavy
                                       // (per-pair fixed term keeps even tiny
                                       //  overlaps in the hard-constraint tier)

        int n = (int)d.blocks.size();
        std::vector<double> sx(n), sy(n);
        for (int i = 0; i < n; i++) { sx[i] = bst.x[i]; sy[i] = bst.y[i]; }

        for (int i : edge_block_idx) {
            int li = std::min(active_loc[i], (int)d.blocks[i].locations.size() - 1);
            const std::string& loc = d.blocks[i].locations[li];

            if (loc.find('T') != std::string::npos) sy[i] = max_h - H[i];
            if (loc.find('B') != std::string::npos) sy[i] = 0.0;
            if (loc.find('L') != std::string::npos) sx[i] = 0.0;
            if (loc.find('R') != std::string::npos) sx[i] = max_w - W[i];

            penalty += W_DIST * (std::abs(sx[i] - bst.x[i]) / max_w +
                                 std::abs(sy[i] - bst.y[i]) / max_h);
        }

        for (int i : edge_block_idx) {
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                if (d.blocks[j].type == BlockType::EDGE && i >= j) continue;
                double ox = std::min(sx[i] + W[i], sx[j] + W[j]) - std::max(sx[i], sx[j]);
                double oy = std::min(sy[i] + H[i], sy[j] + H[j]) - std::max(sy[i], sy[j]);
                if (ox > 1e-6 && oy > 1e-6) penalty += W_OVL * (1.0 + (ox * oy) / Anorm);
            }
        }
        return penalty;
    }
};
