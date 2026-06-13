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

    // Shape-term weights, per-instance so the worker portfolio can run half
    // the searches with them and half without: they trade base packing
    // quality for artery/moat structure, which wins only on cases where that
    // structure is the binding constraint (e.g. b50u65), and the penalty-aware
    // final score picks the better strategy per case.
    double ftcw  = cfg::FTCW;
    double moatw = cfg::MOATW;

    // padded scratch for HALO evaluation (reused, no per-eval alloc)
    std::vector<double> pad_W, pad_H;

    // scratch for the FT-concentration and cut-overflow cost terms (reused;
    // mutable because compute_cost is const)
    mutable std::vector<double> ft_load_scratch;
    struct CutEv { double pos; double dcap; double ddem; };
    mutable std::vector<CutEv> cut_evs;

    // Routed-feedback artery marks: blocks the REAL router overloaded beyond
    // their absorbable capacity, recorded with the position where it happened.
    // The SA cost repels a marked block from its marked spot — an
    // estimate-free signal (set at checkpoints by the SA, see sa_optimizer.h).
    struct ArteryMark { int blk; double px, py; double sev; };
    std::vector<ArteryMark> artery_marks;

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
            int li = std::min(active_loc[i], (int)b.locations.size() - 1);
            const std::string& loc = b.locations[li];

            double ox = b.lx, oy = b.ly, nx = b.lx, ny = b.ly;
            if (loc.find('T') != std::string::npos) ny = max_h - b.height;
            if (loc.find('B') != std::string::npos) ny = 0.0;
            if (loc.find('L') != std::string::npos) nx = 0.0;
            if (loc.find('R') != std::string::npos) nx = max_w - b.width;

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

    double compute_hpwl() const {
        double total = 0;
        for (auto& conn : d.connections) {
            int a = conn.from, b = conn.to;
            double cx_a = bst.x[a] + W[a] * 0.5, cy_a = bst.y[a] + H[a] * 0.5;
            double cx_b = bst.x[b] + W[b] * 0.5, cy_b = bst.y[b] + H[b] * 0.5;
            total += conn.nets * (std::abs(cx_b - cx_a) + std::abs(cy_b - cy_a));
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
            auto fc = ftest::cost_conc(d.blocks, d.connections, bst.x, bst.y, W, H,
                                       cfg::FT_EST_CAP, ft_load_scratch);
            penalty += ftw * (fc.total / Fnorm);
            // Concentration: total FT being low is not enough — the same mass
            // piled onto one artery block is a guaranteed unfixable penalty,
            // spread across many blocks it is absorbed by in-place expansion.
            if (ftcw > 0.0)
                penalty += ftcw * (fc.excess / Fnorm);
        }

        // Routing-aware: keep high-demand pairs close (wide connecting channel).
        if (cfg::PCONGW > 0.0)
            penalty += cfg::PCONGW * (compute_congestion() / Wnorm);

        // Cut overflow: demand straddling a vertical/horizontal cut beyond the
        // cut's carrying capacity cannot be routed by ANY router — the blocks
        // must move to the same side.  Normalized per direction by the
        // full-die cut capacity so the term is O(1).
        if (moatw > 0.0) {
            double ox = cut_overflow_axis(true,  chip_w, chip_h);
            double oy = cut_overflow_axis(false, chip_w, chip_h);
            penalty += moatw * (ox / (25.0 * chip_h) + oy / (25.0 * chip_w));
        }

        // Routed-feedback repulsion: a marked block sitting near the spot
        // where the real router overloaded it keeps becoming an artery —
        // push it away; whoever replaces it has different size/capacity and
        // the next checkpoint re-judges the result.
        if (cfg::RFBW > 0.0 && !artery_marks.empty()) {
            double R = 0.25 * (max_w + max_h);
            for (const auto& m : artery_marks) {
                double cx = bst.x[m.blk] + W[m.blk] * 0.5;
                double cy = bst.y[m.blk] + H[m.blk] * 0.5;
                double dist = std::abs(cx - m.px) + std::abs(cy - m.py);
                if (dist < R) penalty += cfg::RFBW * m.sev * (1.0 - dist / R);
            }
        }

        penalty += edge_block_penalty(max_w, max_h);
        return base + penalty;
    }

    // Worst-cut overflow along one axis (vertical=true sweeps x-cuts that
    // left-right demand must cross).  Sweep events: a block crossing the cut
    // removes 25*span of channel capacity and (if SOFT) credits FT_EST_CAP of
    // absorbable feedthrough; a connection adds its nets between its two
    // block centers.  Returns max over cuts of (demand - capacity), >= 0.
    double cut_overflow_axis(bool vertical, double chip_w, double chip_h) const {
        int n = (int)d.blocks.size();
        double span_total = vertical ? chip_h : chip_w;
        cut_evs.clear();
        for (int i = 0; i < n; i++) {
            double lo   = vertical ? bst.x[i] : bst.y[i];
            double len  = vertical ? W[i] : H[i];
            double hgt  = vertical ? H[i] : W[i];
            double ftc  = (d.blocks[i].type == BlockType::SOFT) ? cfg::FT_EST_CAP : 0.0;
            double dcap = -25.0 * hgt + ftc; // crossing block: less channel, some FT
            cut_evs.push_back({lo,        dcap, 0.0});
            cut_evs.push_back({lo + len, -dcap, 0.0});
        }
        for (auto& c : d.connections) {
            int a = c.from, b = c.to;
            double ca = vertical ? bst.x[a] + W[a] * 0.5 : bst.y[a] + H[a] * 0.5;
            double cb = vertical ? bst.x[b] + W[b] * 0.5 : bst.y[b] + H[b] * 0.5;
            if (ca > cb) std::swap(ca, cb);
            cut_evs.push_back({ca,  0.0,  (double)c.nets});
            cut_evs.push_back({cb,  0.0, -(double)c.nets});
        }
        std::sort(cut_evs.begin(), cut_evs.end(),
                  [](const CutEv& a, const CutEv& b){ return a.pos < b.pos; });

        double cap_delta = 0, dem = 0, worst = 0;
        for (size_t k = 0; k < cut_evs.size(); k++) {
            cap_delta += cut_evs[k].dcap;
            dem       += cut_evs[k].ddem;
            // measure on the open segment after this event group
            if (k + 1 < cut_evs.size() && cut_evs[k+1].pos - cut_evs[k].pos < 1e-9)
                continue;
            double cap = 25.0 * span_total + cap_delta;
            worst = std::max(worst, dem - cap);
        }
        return worst;
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
