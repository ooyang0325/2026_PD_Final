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

// When direct growth fails because block i is boxed in, try shoving the single
// blocking neighbor in one direction into its own free room, then regrow.
// EDGE blocks are never pushed (they must stay snapped to their boundary).
// Single-level only — no recursive push chains — and the shifted rect is
// overlap-validated exactly, so the caller's score gate is the only judge of
// whether the result is actually better.
inline bool push_and_grow(Design& d, Floorplan& fp, int i, double tgt) {
    const Block& b = d.blocks[i];
    for (int dir = 0; dir < 4; dir++) { // 0=R 1=L 2=U 3=D
        bool horiz = (dir < 2);
        Slack s = slack_of(d, i);
        double maxW = b.width  + s.L + s.R;
        double maxH = b.height + s.D + s.U;
        // Room still missing in this axis if the other axis is maxed out.
        double need = horiz ? (tgt / maxH - maxW) : (tgt / maxW - maxH);
        if (need <= 0) continue;       // axis not the bottleneck
        need = ceil2(need) + 0.02;     // grid + abutment margin

        // The unique neighbor walling off this direction (within i's span).
        int wall = -1;
        for (int j = 0; j < (int)d.blocks.size(); j++) {
            if (j == i) continue;
            const Block& o = d.blocks[j];
            bool span_ovl = horiz
                ? (o.ly < b.ly + b.height - 1e-9 && o.ly + o.height > b.ly + 1e-9)
                : (o.lx < b.lx + b.width  - 1e-9 && o.lx + o.width  > b.lx + 1e-9);
            if (!span_ovl) continue;
            bool on_side =
                (dir == 0) ? (o.lx >= b.lx + b.width - 1e-9 &&
                              o.lx < b.lx + b.width + s.R + need + 1e-9) :
                (dir == 1) ? (o.lx + o.width <= b.lx + 1e-9 &&
                              o.lx + o.width > b.lx - s.L - need - 1e-9) :
                (dir == 2) ? (o.ly >= b.ly + b.height - 1e-9 &&
                              o.ly < b.ly + b.height + s.U + need + 1e-9)
                           : (o.ly + o.height <= b.ly + 1e-9 &&
                              o.ly + o.height > b.ly - s.D - need - 1e-9);
            if (!on_side) continue;
            if (wall >= 0) { wall = -2; break; } // more than one: too complex
            wall = j;
        }
        if (wall < 0) continue;
        Block& w = d.blocks[wall];
        if (w.type == BlockType::EDGE) continue;

        double nlx = w.lx + ((dir == 0) ? need : (dir == 1) ? -need : 0.0);
        double nly = w.ly + ((dir == 2) ? need : (dir == 3) ? -need : 0.0);
        nlx = (dir == 1) ? floor2(nlx) : nlx;
        nly = (dir == 3) ? floor2(nly) : nly;
        if (!rect_ok(d, wall, nlx, nly, w.width, w.height)) continue;

        double olx = w.lx, oly = w.ly;
        w.lx = nlx; w.ly = nly;
        if (grow_one(d, fp, i, tgt)) return true;
        w.lx = olx; w.ly = oly; // growth still failed: undo the shove
    }
    return false;
}

// rect_ok, but ignoring two resident blocks (used when evaluating a swap).
inline bool rect_ok2(const Design& d, int skip_a, int skip_b,
                     double lx, double ly, double w, double h) {
    if (lx < -1e-9 || ly < -1e-9) return false;
    if (lx + w > d.outline.cur_width  + 1e-9) return false;
    if (ly + h > d.outline.cur_height + 1e-9) return false;
    for (int j = 0; j < (int)d.blocks.size(); j++) {
        if (j == skip_a || j == skip_b) continue;
        const Block& o = d.blocks[j];
        double ox = std::min(lx + w, o.lx + o.width)  - std::max(lx, o.lx);
        double oy = std::min(ly + h, o.ly + o.height) - std::max(ly, o.ly);
        if (ox > 1e-6 && oy > 1e-6) return false;
    }
    return true;
}

// Demand mass whose endpoint-to-endpoint line crosses the given rect —
// a cheap stand-in for "how much traffic would feed through a block here".
inline double traffic_at(const Design& d, int self,
                         double lx, double ly, double w, double h) {
    double t = 0;
    for (const auto& c : d.connections) {
        if (c.from == self || c.to == self) continue;
        const Block& a = d.blocks[c.from];
        const Block& b = d.blocks[c.to];
        if (ftest::seg_rect(a.lx + a.width * 0.5, a.ly + a.height * 0.5,
                            b.lx + b.width * 0.5, b.ly + b.height * 0.5,
                            lx, ly, w, h))
            t += c.nets;
    }
    return t;
}

// Last resort for an artery block whose FT deficit is far beyond any local
// growth: swap places with the block whose spot has the least crossing
// traffic (ties broken toward roomier spots).  Coordinate-space only; the
// caller's reroute + strict score gate decides whether the move survives.
inline bool relocate_artery(Design& d, Floorplan& fp, int i) {
    const Block& b = d.blocks[i];
    double cur_traffic = traffic_at(d, i, b.lx, b.ly, b.width, b.height);
    int best = -1;
    double best_key = cur_traffic; // must strictly beat staying put
    for (int j = 0; j < (int)d.blocks.size(); j++) {
        if (j == i || d.blocks[j].type == BlockType::EDGE) continue;
        const Block& v = d.blocks[j];
        if (!rect_ok2(d, i, j, v.lx, v.ly, b.width, b.height)) continue;
        if (!rect_ok2(d, i, j, b.lx, b.ly, v.width, v.height)) continue;
        double t = traffic_at(d, i, v.lx, v.ly, b.width, b.height);
        if (t < best_key - 1e-9) { best_key = t; best = j; }
    }
    if (best < 0) return false;
    std::swap(d.blocks[i].lx, d.blocks[best].lx);
    std::swap(d.blocks[i].ly, d.blocks[best].ly);
    return true;
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
        if (grow_one(d, fp, i, tgt + 1.0)) { grown++; continue; }
        if (push_and_grow(d, fp, i, tgt + 1.0)) { grown++; continue; }
        // Hopeless in place: move the artery off the traffic spine and let
        // the next grow→reroute round size it in its new spot.
        if (relocate_artery(d, fp, i)) grown++;
    }
    return grown;
}

} // namespace inplace
