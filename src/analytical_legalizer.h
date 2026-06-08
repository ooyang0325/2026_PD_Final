#pragma once
#include "types.h"
#include "floorplan.h"
#include "channel.h"
#include "router.h"
#include "config.h"
#include <vector>
#include <cmath>
#include <algorithm>

// ─── Final-pass density-driven analytical legalizer ─────────────────────────
//
// Runs AFTER the legalize loop has plateaued.  Breaks the B*-tree topology by
// design (operates on (lx, ly) directly) — user-approved.
//
// Two phases:
//
//   PHASE 1 — Density-driven spreading.
//     A 2D grid is laid over the current compact outline.  Each cell's density
//     counts (a) block-area intersection in the cell, (b) overflow magnitude
//     of any channel passing through.  Each non-pinned block moves OPPOSITE
//     the local density gradient — i.e. toward less-crowded territory — by a
//     bounded step (≤ 20% of cell width per iter).  This OPENS UP free area
//     around currently-cramped soft blocks without trying to also minimize
//     HPWL (the SA + legalize loop already did the HPWL work).
//
//   PHASE 2 — Per-block expansion (legalisation of FT deficit).
//     With newly-freed space adjacent to undersized soft blocks, try resizing
//     each one to its FT-implied target.  Per-block strict rollback: if the
//     resize overlaps a neighbour or exits the outline, revert just that
//     block's W/H.  Catches the structural penalty mode this case suffers from
//     ("16 of 25 penalties are soft-block-undersize" — the legalize loop can't
//     expand because there's no immediately-adjacent free space; spreading
//     creates it).
//
// Pinning:
//   • EDGE blocks: pinned (user direction).
//   • HARD_MACRO: pinned (fixed dimensions — moving them often re-creates the
//                same hot spot somewhere else).
//   • SOFT: free to move and resize.
//
// Caller-side strict rollback covers the whole pass — if total penalty doesn't
// strictly improve, the snapshot is restored and ANA's work is discarded.

class AnalyticalLegalizer {
public:
    Floorplan& fp;
    Design& d;

    int    iterations = 25;     // spreading iterations
    double step_frac  = 0.2;    // fraction of cell size moved per iter
    int    ovl_iters  = 8;      // MTV overlap sweeps per spread iter
    int    grid_w     = 24;     // density-grid resolution
    int    grid_h     = 24;
    double cong_w     = 100.0;  // congestion-source weight in density grid
    double repel_w    = 3.0;    // unused, kept so existing knob wiring doesn't break

    AnalyticalLegalizer(Floorplan& fp_, Design& d_) : fp(fp_), d(d_) {}

    // Returns true if any block moved or resized.
    bool run() {
        const int n = (int)d.blocks.size();
        if (n == 0) return false;

        // Working coords.
        std::vector<double> x(n), y(n);
        std::vector<bool>   pinned(n, false);
        for (int i = 0; i < n; i++) { x[i] = d.blocks[i].lx; y[i] = d.blocks[i].ly; }
        for (int i : fp.edge_block_idx) pinned[i] = true;
        for (int i = 0; i < n; i++)
            if (d.blocks[i].type == BlockType::HARD_MACRO) pinned[i] = true;

        // Use the CURRENT compact outline as the world.  Edge blocks pinned at
        // their current boundary positions; spreading clipped within this box.
        const double mw = (d.outline.cur_width  > 0) ? d.outline.cur_width  : d.outline.max_width;
        const double mh = (d.outline.cur_height > 0) ? d.outline.cur_height : d.outline.max_height;
        const double cell_w = mw / grid_w;
        const double cell_h = mh / grid_h;
        if (cell_w <= 0 || cell_h <= 0) return false;
        const double max_move = step_frac * std::min(cell_w, cell_h);

        // Congestion sources for density-map weighting: each overflowed
        // channel contributes mag × cong_w to its centre cell.
        struct CongSrc { double cx, cy, mag; };
        std::vector<CongSrc> srcs;
        for (auto& ch : d.channels) {
            double ox = std::max(0.0, ch.nets_x - ch.cap_x());
            double oy = std::max(0.0, ch.nets_y - ch.cap_y());
            double mag = ox + oy;
            if (mag <= 1e-3) continue;
            srcs.push_back({ch.lx + ch.width  * 0.5,
                            ch.ly + ch.height * 0.5,
                            mag});
        }

        // Any soft-block deficit? — if not, nothing for phase 2 to do either.
        bool has_deficit = false;
        for (int i = 0; i < n; i++) {
            if (d.blocks[i].type != BlockType::SOFT) continue;
            if (fp.ft_nets[i] <= 0) continue;
            double req = d.blocks[i].get_target_area(fp.ft_nets[i]);
            double act = fp.W[i] * fp.H[i];
            if (act < req - 1.0) { has_deficit = true; break; }
        }
        if (srcs.empty() && !has_deficit) return false;

        bool any_moved = false;
        std::vector<double> density;

        // ── PHASE 1: density-driven spreading ───────────────────────────────
        for (int it = 0; it < iterations; it++) {
            // (1) Build density grid: sum of block-area intersections per cell.
            density.assign(grid_w * grid_h, 0.0);
            for (int i = 0; i < n; i++) {
                double bx0 = x[i],            bx1 = x[i] + fp.W[i];
                double by0 = y[i],            by1 = y[i] + fp.H[i];
                int    gx0 = std::max(0, (int)std::floor(bx0 / cell_w));
                int    gx1 = std::min(grid_w - 1, (int)std::floor((bx1 - 1e-9) / cell_w));
                int    gy0 = std::max(0, (int)std::floor(by0 / cell_h));
                int    gy1 = std::min(grid_h - 1, (int)std::floor((by1 - 1e-9) / cell_h));
                for (int gy = gy0; gy <= gy1; gy++) {
                    for (int gx = gx0; gx <= gx1; gx++) {
                        double cx0 = gx * cell_w, cx1 = (gx + 1) * cell_w;
                        double cy0 = gy * cell_h, cy1 = (gy + 1) * cell_h;
                        double iw = std::min(bx1, cx1) - std::max(bx0, cx0);
                        double ih = std::min(by1, cy1) - std::max(by0, cy0);
                        if (iw > 0 && ih > 0) density[gy * grid_w + gx] += iw * ih;
                    }
                }
            }
            // (2) Add congestion-source weight at overflowed channel centres.
            for (auto& s : srcs) {
                int gx = (int)std::max(0.0, std::min((double)(grid_w - 1), s.cx / cell_w));
                int gy = (int)std::max(0.0, std::min((double)(grid_h - 1), s.cy / cell_h));
                density[gy * grid_w + gx] += s.mag * cong_w;
            }

            // (3) For each non-pinned block, move opposite the local density
            // gradient.  Step capped to step_frac × cell_size.
            bool moved_this_iter = false;
            for (int i = 0; i < n; i++) {
                if (pinned[i]) continue;
                double cx = x[i] + fp.W[i] * 0.5;
                double cy = y[i] + fp.H[i] * 0.5;
                int gx = (int)std::max(0.0, std::min((double)(grid_w - 1), cx / cell_w));
                int gy = (int)std::max(0.0, std::min((double)(grid_h - 1), cy / cell_h));

                double gxl = (gx > 0)            ? density[gy * grid_w + (gx - 1)] : density[gy * grid_w + gx];
                double gxr = (gx < grid_w - 1)   ? density[gy * grid_w + (gx + 1)] : density[gy * grid_w + gx];
                double gyl = (gy > 0)            ? density[(gy - 1) * grid_w + gx] : density[gy * grid_w + gx];
                double gyr = (gy < grid_h - 1)   ? density[(gy + 1) * grid_w + gx] : density[gy * grid_w + gx];

                double gdx = gxr - gxl;
                double gdy = gyr - gyl;
                double gmag = std::sqrt(gdx * gdx + gdy * gdy);
                if (gmag < 1e-9) continue;
                double ux = -gdx / gmag, uy = -gdy / gmag; // unit vector down gradient
                double nx = std::max(0.0, std::min(mw - fp.W[i], x[i] + ux * max_move));
                double ny = std::max(0.0, std::min(mh - fp.H[i], y[i] + uy * max_move));
                if (std::abs(nx - x[i]) > 1e-6 || std::abs(ny - y[i]) > 1e-6) {
                    moved_this_iter = true;
                    any_moved = true;
                }
                x[i] = nx; y[i] = ny;
            }

            // (4) MTV-based overlap resolution (same as before).
            resolve_overlaps(x, y, pinned, mw, mh);

            if (!moved_this_iter) break; // converged
        }

        // Commit phase-1 positions to d.blocks (with 0.01 µm snap).
        for (int i = 0; i < n; i++) {
            if (pinned[i]) continue;
            d.blocks[i].lx = std::round(x[i] * 100.0) / 100.0;
            d.blocks[i].ly = std::round(y[i] * 100.0) / 100.0;
            if (d.blocks[i].lx < 0) d.blocks[i].lx = 0;
            if (d.blocks[i].ly < 0) d.blocks[i].ly = 0;
            if (d.blocks[i].lx + d.blocks[i].width  > mw)
                d.blocks[i].lx = mw - d.blocks[i].width;
            if (d.blocks[i].ly + d.blocks[i].height > mh)
                d.blocks[i].ly = mh - d.blocks[i].height;
        }

        // ── PHASE 2: per-block expansion of undersized soft blocks ──────────
        // For each soft block that's still undersized, compute its FT-implied
        // target W/H and TRY to set fp.W/H + d.blocks.width/height to it.  If
        // the new size overlaps a neighbour or exits the outline, revert just
        // this block.  Otherwise keep.  Pass over all blocks in deficit order
        // so the worst-undersized blocks get first dibs on the freed space.
        std::vector<std::pair<double, int>> by_deficit;
        for (int i = 0; i < n; i++) {
            if (d.blocks[i].type != BlockType::SOFT) continue;
            if (fp.ft_nets[i] <= 0) continue;
            double req = d.blocks[i].get_target_area(fp.ft_nets[i]);
            double act = fp.W[i] * fp.H[i];
            double def = req - act;
            if (def > 1.0) by_deficit.push_back({def, i});
        }
        std::sort(by_deficit.begin(), by_deficit.end(),
                  [](auto& a, auto& b){ return a.first > b.first; });

        for (auto& kv : by_deficit) {
            int i = kv.second;
            double req = d.blocks[i].get_target_area(fp.ft_nets[i]);
            double cur_ar = (fp.H[i] > 0) ? fp.W[i] / fp.H[i] : 1.0;
            cur_ar = std::max(d.blocks[i].min_ar,
                       std::min(d.blocks[i].max_ar, cur_ar));
            double new_w = std::ceil(std::sqrt(req * cur_ar) * 100.0) / 100.0;
            double new_h = std::ceil((req / new_w) * 100.0) / 100.0;

            double oldW = fp.W[i], oldH = fp.H[i];
            double oldBW = d.blocks[i].width, oldBH = d.blocks[i].height;
            fp.W[i] = new_w; fp.H[i] = new_h;
            d.blocks[i].width = new_w; d.blocks[i].height = new_h;

            // Check overlap with any other block and outline.
            bool ok = true;
            if (d.blocks[i].lx + new_w > mw + 0.01) ok = false;
            if (d.blocks[i].ly + new_h > mh + 0.01) ok = false;
            if (ok) {
                for (int j = 0; j < n && ok; j++) {
                    if (j == i) continue;
                    double ox = std::min(d.blocks[i].lx + d.blocks[i].width,
                                         d.blocks[j].lx + d.blocks[j].width)
                              - std::max(d.blocks[i].lx, d.blocks[j].lx);
                    double oy = std::min(d.blocks[i].ly + d.blocks[i].height,
                                         d.blocks[j].ly + d.blocks[j].height)
                              - std::max(d.blocks[i].ly, d.blocks[j].ly);
                    if (ox > 0.01 && oy > 0.01) ok = false;
                }
            }
            if (!ok) {
                fp.W[i] = oldW; fp.H[i] = oldH;
                d.blocks[i].width = oldBW; d.blocks[i].height = oldBH;
            } else {
                any_moved = true;
            }
        }

        return any_moved;
    }

private:
    void resolve_overlaps(std::vector<double>& x, std::vector<double>& y,
                          const std::vector<bool>& pinned,
                          double mw, double mh) {
        const int n = (int)x.size();
        const double SEP = 0.1;
        for (int it = 0; it < ovl_iters; it++) {
            bool any = false;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double ox = std::min(x[i] + fp.W[i], x[j] + fp.W[j])
                              - std::max(x[i],          x[j]);
                    double oy = std::min(y[i] + fp.H[i], y[j] + fp.H[j])
                              - std::max(y[i],          y[j]);
                    if (ox <= 1e-6 || oy <= 1e-6) continue;
                    any = true;

                    bool pi = pinned[i], pj = pinned[j];
                    if (pi && pj) continue;

                    if (ox < oy) {
                        double push = ox + SEP;
                        double si = pi ? 0.0 : (pj ? 1.0 : 0.5);
                        double sj = pj ? 0.0 : (pi ? 1.0 : 0.5);
                        if (x[i] + fp.W[i] * 0.5 < x[j] + fp.W[j] * 0.5) {
                            if (!pi) x[i] -= push * si;
                            if (!pj) x[j] += push * sj;
                        } else {
                            if (!pi) x[i] += push * si;
                            if (!pj) x[j] -= push * sj;
                        }
                    } else {
                        double push = oy + SEP;
                        double si = pi ? 0.0 : (pj ? 1.0 : 0.5);
                        double sj = pj ? 0.0 : (pi ? 1.0 : 0.5);
                        if (y[i] + fp.H[i] * 0.5 < y[j] + fp.H[j] * 0.5) {
                            if (!pi) y[i] -= push * si;
                            if (!pj) y[j] += push * sj;
                        } else {
                            if (!pi) y[i] += push * si;
                            if (!pj) y[j] -= push * sj;
                        }
                    }
                    if (!pi) {
                        x[i] = std::max(0.0, std::min(mw - fp.W[i], x[i]));
                        y[i] = std::max(0.0, std::min(mh - fp.H[i], y[i]));
                    }
                    if (!pj) {
                        x[j] = std::max(0.0, std::min(mw - fp.W[j], x[j]));
                        y[j] = std::max(0.0, std::min(mh - fp.H[j], y[j]));
                    }
                }
            }
            if (!any) break;
        }
    }
};
