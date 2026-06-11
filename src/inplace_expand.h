#pragma once
#include "types.h"
#include "floorplan.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

// ─── In-place soft-block FT expansion ───────────────────────────────────────
//
// Grows undersized SOFT blocks into adjacent whitespace WITHOUT repacking the
// B*-tree.  The SA packs at base area and its best layouts already span the
// full max outline, so both repack-based expansion paths (legacy ft_iter loop,
// LegalizeLoop EXPAND) push the packing out of the outline and roll back —
// blocks ship at base area and every fed-through soft block costs a penalty.
// Growing in place consumes only the halo/whitespace beside the block, which
// is exactly the area the FT conversion is meant to claim.
//
// Coordinates become authoritative after this pass (same contract as the
// analytical legalizer): callers must not fp.pack() afterwards.
//
// All produced coords/dims sit on the 0.01 output grid so the 2-dp rounding
// in OutputWriter cannot create an overlap or shave the satisfied area.

namespace inplace {

inline double ceil2(double v)  { return std::ceil (v * 100.0 - 1e-7) / 100.0; }
inline double floor2(double v) { return std::floor(v * 100.0 + 1e-7) / 100.0; }

// Free room in each direction before hitting a neighbor or the outline, for
// the current rect of block i.  Guides candidate generation only; every
// candidate is overlap-validated exactly before being committed.
struct Slack { double L = 0, R = 0, D = 0, U = 0; };

inline Slack slack_of(const Design& d, int i) {
    const Block& b = d.blocks[i];
    Slack s;
    s.L = b.lx;
    s.R = d.outline.cur_width  - (b.lx + b.width);
    s.D = b.ly;
    s.U = d.outline.cur_height - (b.ly + b.height);
    for (int j = 0; j < (int)d.blocks.size(); j++) {
        if (j == i) continue;
        const Block& o = d.blocks[j];
        bool yo = o.ly < b.ly + b.height - 1e-9 && o.ly + o.height > b.ly + 1e-9;
        bool xo = o.lx < b.lx + b.width  - 1e-9 && o.lx + o.width  > b.lx + 1e-9;
        if (yo) {
            if (o.lx + o.width <= b.lx + 1e-9)      s.L = std::min(s.L, b.lx - (o.lx + o.width));
            else if (o.lx >= b.lx + b.width - 1e-9) s.R = std::min(s.R, o.lx - (b.lx + b.width));
        }
        if (xo) {
            if (o.ly + o.height <= b.ly + 1e-9)      s.D = std::min(s.D, b.ly - (o.ly + o.height));
            else if (o.ly >= b.ly + b.height - 1e-9) s.U = std::min(s.U, o.ly - (b.ly + b.height));
        }
    }
    s.L = std::max(0.0, s.L); s.R = std::max(0.0, s.R);
    s.D = std::max(0.0, s.D); s.U = std::max(0.0, s.U);
    return s;
}

// Exact validation: candidate rect stays inside the outline and overlaps no
// other block (diagonal neighbors included, which the per-axis slack misses).
inline bool rect_ok(const Design& d, int i, double lx, double ly, double w, double h) {
    if (lx < -1e-9 || ly < -1e-9) return false;
    if (lx + w > d.outline.cur_width  + 1e-9) return false;
    if (ly + h > d.outline.cur_height + 1e-9) return false;
    for (int j = 0; j < (int)d.blocks.size(); j++) {
        if (j == i) continue;
        const Block& o = d.blocks[j];
        double ox = std::min(lx + w, o.lx + o.width)  - std::max(lx, o.lx);
        double oy = std::min(ly + h, o.ly + o.height) - std::max(ly, o.ly);
        if (ox > 1e-6 && oy > 1e-6) return false;
    }
    return true;
}

// Grow block i to at least `tgt` area.  Tries a few height candidates
// (widen-only, keep-aspect, tallest) × two shift distributions and commits
// the first one that validates.  Never shrinks either dimension.
inline bool grow_one(Design& d, Floorplan& fp, int i, double tgt) {
    Block& b = d.blocks[i];
    Slack s = slack_of(d, i);
    double maxW = b.width  + s.L + s.R;
    double maxH = b.height + s.D + s.U;
    if (maxW * maxH < tgt - 1e-9) return false; // not reachable here

    double cur_ar = (b.height > 0) ? b.width / b.height : 1.0;
    double cands[3] = { b.height,                                 // widen only
                        std::sqrt(tgt / std::max(cur_ar, 1e-9)),  // keep aspect
                        maxH };                                   // tallest
    for (double hc : cands) {
        double nh = ceil2(std::min(std::max(hc, b.height), maxH));
        if (nh < b.height) nh = b.height;
        double nw = ceil2(tgt / nh);
        if (nw < b.width) nw = b.width;
        if (nw > maxW + 1e-9) continue;
        double ar = nw / nh;
        if (ar < b.min_ar - 1e-6 || ar > b.max_ar + 1e-6) continue;

        for (int mode = 0; mode < 2; mode++) {
            // mode 0: extend right/up first; mode 1: left/down first.
            double gw = nw - b.width, gh = nh - b.height;
            double dr = (mode == 0) ? std::min(s.R, gw) : gw - std::min(s.L, gw);
            double dl = gw - dr;
            double du = (mode == 0) ? std::min(s.U, gh) : gh - std::min(s.D, gh);
            double dd = gh - du;
            if (dl > s.L + 1e-9 || dr > s.R + 1e-9) continue;
            if (dd > s.D + 1e-9 || du > s.U + 1e-9) continue;
            double nlx = floor2(b.lx - dl);
            double nly = floor2(b.ly - dd);
            if (!rect_ok(d, i, nlx, nly, nw, nh)) continue;
            b.lx = nlx; b.ly = nly; b.width = nw; b.height = nh;
            fp.W[i] = nw; fp.H[i] = nh;
            return true;
        }
    }
    return false;
}

// One pass over all undersized SOFT blocks, smallest deficit first so the
// cheap guaranteed wins land before the big ones eat the whitespace.
// Returns the number of blocks grown.
inline int expand_pass(Design& d, Floorplan& fp) {
    std::vector<std::pair<double,int>> work; // (deficit, block idx)
    for (int i = 0; i < (int)d.blocks.size(); i++) {
        if (d.blocks[i].type != BlockType::SOFT || fp.ft_nets[i] <= 0) continue;
        double tgt = d.blocks[i].get_target_area(fp.ft_nets[i]);
        double act = d.blocks[i].width * d.blocks[i].height;
        if (act < tgt - 1.0) work.push_back({tgt - act, i});
    }
    std::sort(work.begin(), work.end());
    int grown = 0;
    for (auto& [deficit, i] : work) {
        double tgt = d.blocks[i].get_target_area(fp.ft_nets[i]);
        if (grow_one(d, fp, i, tgt + 1.0)) grown++; // +1 area eps for rounding
    }
    return grown;
}

} // namespace inplace
