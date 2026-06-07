#pragma once
#include "types.h"
#include <vector>
#include <cmath>
#include <algorithm>

// ─── FTAFP-style Feedthrough Estimation ─────────────────────────────────────
// Adapted from "FTAFP: A Feedthrough-Aware Floorplanner for Hierarchical Design
// of Large-Scale SoCs" (Li et al., ASPDAC '25), Section 3.1.
//
// The paper operates on multi-pin nets and clusters adjacent modules into
// sub-nets (Algorithm 1, union-find on the "common edge" adjacency relation),
// then estimates feedthrough along an MST over the sub-nets (Algorithm 2, a
// greedy directional traversal that counts the SOFT modules a route passes
// THROUGH).  This project's netlist is pairwise (2-pin nets), so:
//
//   * Net simplification collapses to a single rule: if the two endpoints share
//     a common edge (are directly adjacent), the connection needs no feedthrough
//     (paper rule 1).  Connections across whitespace also need no feedthrough
//     of a module that the route can go *around* (paper rule 2) — captured by
//     only charging SOFT blocks the direct route actually crosses.
//
//   * Feedthrough estimation reduces to: for a non-adjacent endpoint pair, the
//     route is the (Manhattan-ish) line between the two module centers; every
//     SOFT module that line passes through is "fed through" and accrues the
//     net's load.  This is the estimation model that guides the SA; the real
//     channel router later produces the authoritative feedthrough numbers.
//
// This is strictly tighter than a bounding-box overlap test (the previous
// heuristic), because a soft block sitting beside — but not on — the route is
// no longer charged.

namespace ftest {

// Two modules are "directly adjacent" (share a common edge) if they touch /
// nearly touch on one axis while overlapping on the other.  The packing leaves
// a ~HALO gap between neighbors, so a small tolerance is allowed.
inline bool adjacent(double ax, double ay, double aw, double ah,
                     double bx, double by, double bw, double bh,
                     double gap) {
    double ox = std::min(ax + aw, bx + bw) - std::max(ax, bx);
    double oy = std::min(ay + ah, by + bh) - std::max(ay, by);
    double xgap = std::max(ax, bx) - std::min(ax + aw, bx + bw);
    double ygap = std::max(ay, by) - std::min(ay + ah, by + bh);
    if (oy > 1e-6 && xgap <= gap + 1e-6 && xgap >= -1e-6) return true; // horizontal neighbors
    if (ox > 1e-6 && ygap <= gap + 1e-6 && ygap >= -1e-6) return true; // vertical neighbors
    return false;
}

// Segment (x0,y0)-(x1,y1) vs axis-aligned rect [rx,rx+rw] x [ry,ry+rh].
// Liang-Barsky clipping; true iff they intersect.
inline bool seg_rect(double x0, double y0, double x1, double y1,
                     double rx, double ry, double rw, double rh) {
    double dx = x1 - x0, dy = y1 - y0;
    double t0 = 0.0, t1 = 1.0;
    double p[4] = { -dx, dx, -dy, dy };
    double q[4] = { x0 - rx, (rx + rw) - x0, y0 - ry, (ry + rh) - y0 };
    for (int i = 0; i < 4; i++) {
        if (std::abs(p[i]) < 1e-12) {
            if (q[i] < 0) return false;          // parallel and outside slab
        } else {
            double r = q[i] / p[i];
            if (p[i] < 0) { if (r > t1) return false; if (r > t0) t0 = r; }
            else          { if (r < t0) return false; if (r < t1) t1 = r; }
        }
    }
    return t0 <= t1;
}

// Weighted feedthrough estimate for the SA cost: sum over non-adjacent
// connections of nets * (number of SOFT blocks the direct route crosses).
// Minimizing this drives heavily-connected blocks to sit adjacent (so their
// nets route through fat common-edge channels) instead of feeding through soft
// modules — which is what keeps real routing feedthrough (and the consequent
// soft-block area growth) small enough to stay inside the outline.
inline double cost(const std::vector<Block>& blocks,
                   const std::vector<Connection>& conns,
                   const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<double>& W, const std::vector<double>& H) {
    int n = (int)blocks.size();
    const double GAP = 3.0;
    double total = 0.0;
    for (const auto& c : conns) {
        int A = c.from, B = c.to;
        if (adjacent(x[A], y[A], W[A], H[A], x[B], y[B], W[B], H[B], GAP)) continue;
        double cax = x[A] + W[A] * 0.5, cay = y[A] + H[A] * 0.5;
        double cbx = x[B] + W[B] * 0.5, cby = y[B] + H[B] * 0.5;
        int crossed = 0;
        for (int i = 0; i < n; i++) {
            if (i == A || i == B) continue;
            if (blocks[i].type != BlockType::SOFT) continue;
            if (seg_rect(cax, cay, cbx, cby, x[i], y[i], W[i], H[i])) crossed++;
        }
        if (crossed > 0) total += (double)c.nets * crossed;
    }
    return total;
}

// Estimate per-block feedthrough net loads.  Writes ft_nets[i] for SOFT blocks.
inline void estimate(const std::vector<Block>& blocks,
                     const std::vector<Connection>& conns,
                     const std::vector<double>& x, const std::vector<double>& y,
                     const std::vector<double>& W, const std::vector<double>& H,
                     std::vector<int>& ft_nets) {
    int n = (int)blocks.size();
    std::fill(ft_nets.begin(), ft_nets.end(), 0);
    const double GAP = 3.0; // ~HALO tolerance for the adjacency (common-edge) test

    for (const auto& c : conns) {
        int A = c.from, B = c.to;
        if (adjacent(x[A], y[A], W[A], H[A], x[B], y[B], W[B], H[B], GAP))
            continue; // rule 1: directly adjacent -> no feedthrough

        double cax = x[A] + W[A] * 0.5, cay = y[A] + H[A] * 0.5;
        double cbx = x[B] + W[B] * 0.5, cby = y[B] + H[B] * 0.5;

        for (int i = 0; i < n; i++) {
            if (i == A || i == B) continue;
            if (blocks[i].type != BlockType::SOFT) continue;
            if (seg_rect(cax, cay, cbx, cby, x[i], y[i], W[i], H[i]))
                ft_nets[i] += c.nets;
        }
    }
}

} // namespace ftest
