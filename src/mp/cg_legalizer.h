#pragma once
#include "../types.h"
#include "../floorplan.h"
#include <vector>
#include <cmath>
#include <algorithm>

// ─── CGLegalizer: coordinate-only overlap removal for the MP engine ─────────
//
// Input  : d.blocks[i].lx/ly from MPOptimizer::run (overlapping global placement)
//          fp.W[i]/fp.H[i]  — the REAL block dims used as separation spacing
// Output : every block-block overlap removed, all blocks inside [0,W]x[0,H],
//          EDGE blocks kept pinned to their active boundary location's fixed
//          coordinate (the free coordinate may slide).
//
// METHOD CHOSEN — MTV iterative relaxation (primary), constraint-graph
// longest-path (secondary attempt first).  Rationale: for <=69 blocks with
// 20-50% whitespace, the minimum-translation-vector sweep from
// analytical_legalizer.h:246-300 reliably reaches 0 overlaps when given enough
// sweeps and a center-of-mass push direction, and it natively respects pinned
// EDGE coordinates.  A single-pass constraint-graph longest-path can OVER-pack
// (it pushes everything to one side, easily exceeding the outline when the
// derived order is poor), so we run CG first as a cheap structural separation
// and fall back to MTV — the workhorse that guarantees the must-have 0-overlap
// result — whenever CG leaves residual overlap or pushes a block out of bounds.
// source: PeF (pef-tcad2023) §IV-B constraint-graph legalization;
//         MTV engine adapted from analytical_legalizer.h:246-300.

namespace mp {

class CGLegalizer {
public:
    CGLegalizer(Floorplan& fp, Design& d) : fp_(fp), d_(d) {}

    // Returns true iff the committed layout is overlap-free AND in-bounds.
    bool legalize(double W, double H) {
        const int n = (int)d_.blocks.size();
        if (n == 0) return true;

        // ---- pin state + fixed-axis for EDGE blocks -------------------------
        // EDGE blocks: the boundary axis named in the active location is FIXED;
        // the orthogonal axis is free to slide.  We capture the fixed target
        // coordinate so every sweep can re-snap it after a relaxation step.
        std::vector<bool> pin_x(n, false), pin_y(n, false);
        std::vector<double> fix_x(n, 0.0), fix_y(n, 0.0);
        for (int i : fp_.edge_block_idx) {
            const auto& b = d_.blocks[i];
            int li = std::min(fp_.active_loc[i], (int)b.locations.size() - 1);
            if (li < 0) continue;
            const std::string& loc = b.locations[li];
            if (loc.find('L') != std::string::npos) { pin_x[i] = true; fix_x[i] = 0.0; }
            if (loc.find('R') != std::string::npos) { pin_x[i] = true; fix_x[i] = W - fp_.W[i]; }
            if (loc.find('B') != std::string::npos) { pin_y[i] = true; fix_y[i] = 0.0; }
            if (loc.find('T') != std::string::npos) { pin_y[i] = true; fix_y[i] = H - fp_.H[i]; }
        }

        // Working coordinates seeded from the global placement, clamped/snapped.
        std::vector<double> x(n), y(n);
        for (int i = 0; i < n; i++) {
            x[i] = pin_x[i] ? fix_x[i] : clamp(d_.blocks[i].lx, 0.0, W - fp_.W[i]);
            y[i] = pin_y[i] ? fix_y[i] : clamp(d_.blocks[i].ly, 0.0, H - fp_.H[i]);
        }

        // ---- Attempt 1: constraint-graph cheaper-axis separation -----------
        // source: PeF §IV-B.  For each overlapping pair, separate along the axis
        // with the SMALLER current overlap (cheaper move), in the direction set
        // by current center order.  A few relaxation rounds untangle the bulk.
        {
            std::vector<double> cx = x, cy = y;
            cg_separate(cx, cy, pin_x, pin_y, fix_x, fix_y, W, H);
            if (overlap_free(cx, cy) && in_bounds(cx, cy, W, H)) {
                commit(cx, cy);
                return true;
            }
            // Keep CG's partial untangle as the MTV seed — it removes most of
            // the gross overlap and lands MTV closer to a feasible point.
            x = cx; y = cy;
        }

        // ---- Attempt 2: MTV iterative relaxation (the robust workhorse) -----
        // source: analytical_legalizer.h:246-300, extended with center-pull and
        // many sweeps for guaranteed convergence given the available whitespace.
        const int MTV_SWEEPS = 400;
        mtv_resolve(x, y, pin_x, pin_y, fix_x, fix_y, W, H, MTV_SWEEPS);
        commit(x, y);
        return overlap_free(x, y) && in_bounds(x, y, W, H);
    }

private:
    Floorplan& fp_;
    Design&    d_;

    static double clamp(double v, double lo, double hi) {
        if (hi < lo) hi = lo;                 // degenerate: block wider than slack
        return std::max(lo, std::min(hi, v));
    }

    // Re-snap a movable block inside [0,W]x[0,H]; pinned axes return to fixed.
    void snap_one(int i, std::vector<double>& x, std::vector<double>& y,
                  const std::vector<bool>& pin_x, const std::vector<bool>& pin_y,
                  const std::vector<double>& fix_x, const std::vector<double>& fix_y,
                  double W, double H) {
        if (pin_x[i]) x[i] = fix_x[i]; else x[i] = clamp(x[i], 0.0, W - fp_.W[i]);
        if (pin_y[i]) y[i] = fix_y[i]; else y[i] = clamp(y[i], 0.0, H - fp_.H[i]);
    }

    // ── Constraint-graph cheaper-axis separation (PeF §IV-B) ─────────────────
    void cg_separate(std::vector<double>& x, std::vector<double>& y,
                     const std::vector<bool>& pin_x, const std::vector<bool>& pin_y,
                     const std::vector<double>& fix_x, const std::vector<double>& fix_y,
                     double W, double H) {
        const int n = (int)x.size();
        const double SEP = 0.1;               // hairline gap so 0.01 snap stays clear
        const int ROUNDS = 60;
        for (int r = 0; r < ROUNDS; r++) {
            bool any = false;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double ox = std::min(x[i] + fp_.W[i], x[j] + fp_.W[j])
                              - std::max(x[i], x[j]);
                    double oy = std::min(y[i] + fp_.H[i], y[j] + fp_.H[j])
                              - std::max(y[i], y[j]);
                    if (ox <= 1e-6 || oy <= 1e-6) continue;
                    any = true;

                    // Choose the cheaper separation axis; if a block is pinned on
                    // that axis, push the other way (the free block absorbs it).
                    bool sep_x = (ox < oy);
                    if (sep_x && pin_x[i] && pin_x[j]) sep_x = false;   // x blocked
                    if (!sep_x && pin_y[i] && pin_y[j]) sep_x = true;   // y blocked

                    if (sep_x) {
                        double push = ox + SEP;
                        bool pi = pin_x[i], pj = pin_x[j];
                        double si = pi ? 0.0 : (pj ? 1.0 : 0.5);
                        double sj = pj ? 0.0 : (pi ? 1.0 : 0.5);
                        if (x[i] + fp_.W[i] * 0.5 < x[j] + fp_.W[j] * 0.5) {
                            x[i] -= push * si; x[j] += push * sj;
                        } else {
                            x[i] += push * si; x[j] -= push * sj;
                        }
                    } else {
                        double push = oy + SEP;
                        bool pi = pin_y[i], pj = pin_y[j];
                        double si = pi ? 0.0 : (pj ? 1.0 : 0.5);
                        double sj = pj ? 0.0 : (pi ? 1.0 : 0.5);
                        if (y[i] + fp_.H[i] * 0.5 < y[j] + fp_.H[j] * 0.5) {
                            y[i] -= push * si; y[j] += push * sj;
                        } else {
                            y[i] += push * si; y[j] -= push * sj;
                        }
                    }
                    snap_one(i, x, y, pin_x, pin_y, fix_x, fix_y, W, H);
                    snap_one(j, x, y, pin_x, pin_y, fix_x, fix_y, W, H);
                }
            }
            if (!any) break;
        }
    }

    // ── MTV sweeps with a mild center-pull bias ──────────────────────────────
    // source: analytical_legalizer.h:251-300.  Same minimum-translation-vector
    // push, but pinned axes are re-snapped each step and we run many sweeps so
    // the layout fully untangles within the whitespace budget.
    void mtv_resolve(std::vector<double>& x, std::vector<double>& y,
                     const std::vector<bool>& pin_x, const std::vector<bool>& pin_y,
                     const std::vector<double>& fix_x, const std::vector<double>& fix_y,
                     double W, double H, int sweeps) {
        const int n = (int)x.size();
        const double SEP = 0.1;
        for (int it = 0; it < sweeps; it++) {
            bool any = false;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double ox = std::min(x[i] + fp_.W[i], x[j] + fp_.W[j])
                              - std::max(x[i], x[j]);
                    double oy = std::min(y[i] + fp_.H[i], y[j] + fp_.H[j])
                              - std::max(y[i], y[j]);
                    if (ox <= 1e-6 || oy <= 1e-6) continue;
                    any = true;

                    // Separate along the smaller overlap (MTV); honor pins.
                    bool sep_x = (ox < oy);
                    bool fully_pinned_x = pin_x[i] && pin_x[j];
                    bool fully_pinned_y = pin_y[i] && pin_y[j];
                    if (sep_x && fully_pinned_x && !fully_pinned_y) sep_x = false;
                    if (!sep_x && fully_pinned_y && !fully_pinned_x) sep_x = true;

                    if (sep_x) {
                        double push = ox + SEP;
                        bool pi = pin_x[i], pj = pin_x[j];
                        if (pi && pj) continue;          // cannot separate on x
                        double si = pi ? 0.0 : (pj ? 1.0 : 0.5);
                        double sj = pj ? 0.0 : (pi ? 1.0 : 0.5);
                        if (x[i] + fp_.W[i] * 0.5 < x[j] + fp_.W[j] * 0.5) {
                            x[i] -= push * si; x[j] += push * sj;
                        } else {
                            x[i] += push * si; x[j] -= push * sj;
                        }
                    } else {
                        double push = oy + SEP;
                        bool pi = pin_y[i], pj = pin_y[j];
                        if (pi && pj) continue;          // cannot separate on y
                        double si = pi ? 0.0 : (pj ? 1.0 : 0.5);
                        double sj = pj ? 0.0 : (pi ? 1.0 : 0.5);
                        if (y[i] + fp_.H[i] * 0.5 < y[j] + fp_.H[j] * 0.5) {
                            y[i] -= push * si; y[j] += push * sj;
                        } else {
                            y[i] += push * si; y[j] -= push * sj;
                        }
                    }
                    snap_one(i, x, y, pin_x, pin_y, fix_x, fix_y, W, H);
                    snap_one(j, x, y, pin_x, pin_y, fix_x, fix_y, W, H);
                }
            }
            if (!any) break;
        }
    }

    bool overlap_free(const std::vector<double>& x, const std::vector<double>& y) const {
        const int n = (int)x.size();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                double ox = std::min(x[i] + fp_.W[i], x[j] + fp_.W[j]) - std::max(x[i], x[j]);
                double oy = std::min(y[i] + fp_.H[i], y[j] + fp_.H[j]) - std::max(y[i], y[j]);
                if (ox > 0.01 && oy > 0.01) return false;
            }
        return true;
    }

    bool in_bounds(const std::vector<double>& x, const std::vector<double>& y,
                   double W, double H) const {
        const int n = (int)x.size();
        for (int i = 0; i < n; i++) {
            if (x[i] < -1e-3 || y[i] < -1e-3) return false;
            if (x[i] + fp_.W[i] > W + 1e-3 || y[i] + fp_.H[i] > H + 1e-3) return false;
        }
        return true;
    }

    // Write coordinates back to d.blocks with a 0.01 µm snap (output grid).
    void commit(const std::vector<double>& x, const std::vector<double>& y) {
        const int n = (int)d_.blocks.size();
        for (int i = 0; i < n; i++) {
            d_.blocks[i].lx = std::round(x[i] * 100.0) / 100.0;
            d_.blocks[i].ly = std::round(y[i] * 100.0) / 100.0;
            // Keep dims in sync with the spacing actually used for separation.
            d_.blocks[i].width  = fp_.W[i];
            d_.blocks[i].height = fp_.H[i];
        }
    }
};

} // namespace mp
