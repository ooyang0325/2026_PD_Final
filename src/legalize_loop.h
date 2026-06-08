#pragma once
#include "types.h"
#include "floorplan.h"
#include "bstree.h"
#include "channel.h"
#include "router.h"
#include "config.h"
#include <vector>
#include <set>
#include <chrono>
#include <random>
#include <functional>
#include <algorithm>
#include <cmath>

// ─── Route-driven legalize loop ─────────────────────────────────────────────
//
// Post-SA refinement.  Given a packed layout that may carry channel-overflow
// or undersized-soft-block penalties, iteratively shrink the penalty using a
// cheap-first action cascade with STRICT rollback:
//
//   per iteration:
//     1. route + collect hotspots (overflowed channels + undersized soft blocks)
//     2. pick the worst hotspot by magnitude
//     3. try actions in order: ROTATE → DISPLACE → EXPAND
//        - apply, finalize, re-route, score
//        - accept only if score STRICTLY decreases AND layout stays valid
//        - else restore snapshot, try next action
//     4. if no action helps, blacklist this hotspot and try the next worst
//     5. if all remaining hotspots are blacklisted → plateau, terminate
//
// Output guarantees:
//   - Never returns a layout worse than the one passed in (snapshot+rollback).
//   - Never returns an invalid layout (validity check before accepting).
//
// Scoring contract: the caller provides `route_and_score()` returning the same
// scalar used to compare candidate layouts elsewhere in the pipeline
// (compute_final_cost + FAILW*edge_fails + PEN*overflow_count + magnitude).
// The loop never re-implements scoring — it just minimizes whatever is given.

class LegalizeLoop {
public:
    Floorplan& fp;
    Design& d;
    std::function<void()> finalize;       // pack + snap edge blocks + recompute channels
    std::function<double()> score;        // route + full penalty-aware score
    std::function<bool()> is_valid;       // in-bounds + non-overlapping
    std::mt19937 rng;

    // Tunables (env-overridable in config.h via cfg::LEG_*).
    int    max_iters       = 30;
    double max_seconds     = 30.0;
    int    displace_tries  = 8;     // K random perturbations per displace attempt
    double improve_eps     = 1.0;   // require this much absolute score drop to accept

    LegalizeLoop(Floorplan& fp_, Design& d_,
                 std::function<void()> finalize_,
                 std::function<double()> score_,
                 std::function<bool()> is_valid_,
                 unsigned seed = 7777)
        : fp(fp_), d(d_),
          finalize(std::move(finalize_)),
          score(std::move(score_)),
          is_valid(std::move(is_valid_)),
          rng(seed) {}

    // Runs to convergence or budget.  Returns the score of the best-and-current
    // layout (also left in fp/d).  Caller decides what to do if score signals
    // edge-fail (≥ FAILW): the loop itself never accepts such a state once a
    // valid one has been seen.
    double run() {
        auto t0 = std::chrono::steady_clock::now();
        auto elapsed = [&]{
            return std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
        };

        // Baseline: finalize + score the input layout.
        finalize();
        if (!is_valid()) return score();  // caller handles emergency fallback
        double current = score();

        // Blacklists prevent re-attacking a hotspot we've already exhausted.
        std::set<int> bl_ch_x, bl_ch_y, bl_soft;

        for (int it = 0; it < max_iters; it++) {
            if (elapsed() >= max_seconds) break;

            auto hotspots = collect_hotspots(bl_ch_x, bl_ch_y, bl_soft);
            if (hotspots.empty()) break;

            // Worst-first.
            std::sort(hotspots.begin(), hotspots.end(),
                      [](const Hotspot& a, const Hotspot& b){
                          return a.magnitude > b.magnitude;
                      });

            bool any_improvement_this_iter = false;
            for (const Hotspot& h : hotspots) {
                bool improved = false;
                double before = current;

                // Pick the block this hotspot points at (the soft block itself,
                // or a well-chosen adjacent block for a channel hotspot).
                int target = pick_target(h);
                if (target < 0) {
                    mark_blacklist(h, bl_ch_x, bl_ch_y, bl_soft);
                    continue;
                }

                // ── 1) ROTATE (cheapest) ────────────────────────────────────
                improved = try_with_rollback(current, [&]{
                    return action_rotate(target);
                });

                // ── 2) DISPLACE (medium) ────────────────────────────────────
                if (!improved) {
                    improved = try_displace(target, current);
                }

                // ── 3) EXPAND (most expensive) ──────────────────────────────
                if (!improved && h.kind == Hotspot::SOFT_UNDERSIZE) {
                    improved = try_with_rollback(current, [&]{
                        return action_expand(target);
                    });
                }

                if (improved && current < before - improve_eps) {
                    any_improvement_this_iter = true;
                    break; // re-pick worst hotspot on the new layout
                }

                // Nothing worked for this hotspot — blacklist it and move on.
                mark_blacklist(h, bl_ch_x, bl_ch_y, bl_soft);
            }

            if (!any_improvement_this_iter) break; // plateau
        }

        return current;
    }

private:
    // ── Hotspot identification ─────────────────────────────────────────────
    struct Hotspot {
        enum Kind { CHANNEL_X, CHANNEL_Y, SOFT_UNDERSIZE };
        Kind kind;
        int idx;
        double magnitude;
    };

    std::vector<Hotspot> collect_hotspots(const std::set<int>& bl_ch_x,
                                          const std::set<int>& bl_ch_y,
                                          const std::set<int>& bl_soft) const {
        std::vector<Hotspot> out;
        for (int i = 0; i < (int)d.channels.size(); i++) {
            const auto& ch = d.channels[i];
            double ox = ch.nets_x - ch.cap_x();
            double oy = ch.nets_y - ch.cap_y();
            if (ox > 1e-3 && !bl_ch_x.count(i))
                out.push_back({Hotspot::CHANNEL_X, i, ox});
            if (oy > 1e-3 && !bl_ch_y.count(i))
                out.push_back({Hotspot::CHANNEL_Y, i, oy});
        }
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (d.blocks[i].type != BlockType::SOFT) continue;
            if (fp.ft_nets[i] <= 0) continue;
            if (bl_soft.count(i)) continue;
            double req = d.blocks[i].get_target_area(fp.ft_nets[i]);
            double act = fp.W[i] * fp.H[i];
            if (act < req - 1.0)
                out.push_back({Hotspot::SOFT_UNDERSIZE, i, req - act});
        }
        return out;
    }

    static void mark_blacklist(const Hotspot& h,
                                std::set<int>& bl_ch_x,
                                std::set<int>& bl_ch_y,
                                std::set<int>& bl_soft) {
        switch (h.kind) {
            case Hotspot::CHANNEL_X:      bl_ch_x.insert(h.idx); break;
            case Hotspot::CHANNEL_Y:      bl_ch_y.insert(h.idx); break;
            case Hotspot::SOFT_UNDERSIZE: bl_soft.insert(h.idx); break;
        }
    }

    // For a SOFT_UNDERSIZE hotspot the target is the block itself.  For a
    // channel hotspot pick the adjacent SOFT block with the most through-FT
    // (most "responsible" for the overflow); fall back to any adjacent
    // non-EDGE block; finally, anything adjacent.
    int pick_target(const Hotspot& h) const {
        if (h.kind == Hotspot::SOFT_UNDERSIZE) return h.idx;

        const auto& ch = d.channels[h.idx];
        int best = -1, best_score = -1;
        for (int i = 0; i < (int)d.blocks.size(); i++) {
            if (!touches(d.blocks[i], ch)) continue;
            int s;
            if (d.blocks[i].type == BlockType::SOFT)
                s = 10000 + fp.ft_nets[i];
            else if (d.blocks[i].type == BlockType::HARD_MACRO)
                s = 500;
            else // EDGE — least preferred (snap constraint)
                s = 1;
            if (s > best_score) { best_score = s; best = i; }
        }
        return best;
    }

    static bool touches(const Block& b, const Channel& ch) {
        // Shares an edge with the channel (within tolerance).
        bool tx = std::abs(b.lx + b.width - ch.lx) < 1e-3 ||
                  std::abs(ch.lx + ch.width - b.lx) < 1e-3;
        bool ty = std::abs(b.ly + b.height - ch.ly) < 1e-3 ||
                  std::abs(ch.ly + ch.height - b.ly) < 1e-3;
        bool yovl = b.ly < ch.ly + ch.height - 1e-6 && b.ly + b.height > ch.ly + 1e-6;
        bool xovl = b.lx < ch.lx + ch.width  - 1e-6 && b.lx + b.width  > ch.lx + 1e-6;
        return (tx && yovl) || (ty && xovl);
    }

    // ── Snapshot / restore ─────────────────────────────────────────────────
    struct Snapshot {
        BStarTree::State bst_state;
        std::vector<double> W, H;
        std::vector<int> active_loc;
        std::vector<int> ft_nets;
        std::vector<Block> blocks;
        std::vector<Channel> channels;
        std::vector<RoutePath> paths;
        Outline outline;
    };

    Snapshot take_snapshot() const {
        Snapshot s;
        s.bst_state = fp.bst.save();
        s.W = fp.W;
        s.H = fp.H;
        s.active_loc = fp.active_loc;
        s.ft_nets = fp.ft_nets;
        s.blocks = d.blocks;
        s.channels = d.channels;
        s.paths = d.paths;
        s.outline = d.outline;
        return s;
    }

    void restore_snapshot(const Snapshot& s) {
        fp.bst.restore(s.bst_state);
        fp.W = s.W;
        fp.H = s.H;
        fp.active_loc = s.active_loc;
        fp.ft_nets = s.ft_nets;
        d.blocks = s.blocks;
        d.channels = s.channels;
        d.paths = s.paths;
        d.outline = s.outline;
    }

    // ── Generic try-with-rollback ───────────────────────────────────────────
    // Apply f() (which mutates fp/d "in the small"), then finalize + re-route +
    // score.  Accept iff the layout is valid AND the score strictly improves.
    // Otherwise restore.
    template<class F>
    bool try_with_rollback(double& current, F&& f) {
        Snapshot snap = take_snapshot();
        if (!f()) { restore_snapshot(snap); return false; }
        finalize();
        if (!is_valid()) { restore_snapshot(snap); return false; }
        double s = score();
        if (s < current - improve_eps) {
            current = s;
            return true;
        }
        restore_snapshot(snap);
        return false;
    }

    // ── Actions ─────────────────────────────────────────────────────────────

    // Swap W ↔ H for a SOFT block, respecting aspect-ratio limits.
    // Returns false if the block isn't SOFT or rotation is out of range.
    bool action_rotate(int blk) {
        if (d.blocks[blk].type != BlockType::SOFT) return false;
        double nw = fp.H[blk], nh = fp.W[blk];
        double ar = (nh > 0) ? nw / nh : 1.0;
        if (ar < d.blocks[blk].min_ar - 1e-6) return false;
        if (ar > d.blocks[blk].max_ar + 1e-6) return false;
        std::swap(fp.W[blk], fp.H[blk]);
        return true;
    }

    // Grow a SOFT block to exactly satisfy its observed FT demand (or 10% more
    // if it already meets the demand — used as a "make some room" knob).
    bool action_expand(int blk) {
        if (d.blocks[blk].type != BlockType::SOFT) return false;
        double tgt = d.blocks[blk].get_target_area(fp.ft_nets[blk]);
        double cur_area = fp.W[blk] * fp.H[blk];
        if (cur_area >= tgt - 1.0) tgt = cur_area * 1.10; // bump

        double cur_ar = (fp.H[blk] > 0) ? fp.W[blk] / fp.H[blk] : 1.0;
        cur_ar = std::max(d.blocks[blk].min_ar,
                  std::min(d.blocks[blk].max_ar, cur_ar));
        double raw_w = std::sqrt(tgt * cur_ar);
        fp.W[blk] = std::ceil(raw_w * 100.0) / 100.0;
        fp.H[blk] = std::ceil((tgt / fp.W[blk]) * 100.0) / 100.0;
        return true;
    }

    // Mini-SA at T=0 around block `blk`.  Tries up to K random perturbations
    // (swap or move) involving the node currently holding `blk` and keeps the
    // best valid candidate.  Accept iff that best STRICTLY improves over
    // `current`; otherwise restore.
    bool try_displace(int blk, double& current) {
        if (fp.bst.n < 2) return false;

        Snapshot start = take_snapshot();
        Snapshot best  = start;
        double best_s  = current; // must beat to win
        bool found     = false;

        std::uniform_int_distribution<int> nd(0, fp.bst.n - 1);
        for (int t = 0; t < displace_tries; t++) {
            restore_snapshot(start);

            // Re-locate the node currently holding `blk` (delete-to-leaf
            // shuffles blk[]; after each restore the topology is the same so
            // the node id is stable, but we look it up cleanly anyway).
            int node = -1;
            for (int n = 0; n < fp.bst.n; n++)
                if (fp.bst.blk[n] == blk) { node = n; break; }
            if (node < 0) continue;

            if (rng() & 1u) {
                // Swap with a random other node.
                int other = nd(rng);
                while (other == node) other = nd(rng);
                fp.bst.swap_nodes(node, other);
            } else {
                // Delete the node holding `blk` to a leaf, then re-insert at a
                // random parent.  Migrating `blk` down the tree relocates it
                // into a different neighborhood in the packing.
                int freed = fp.bst.delete_to_leaf(node, rng);
                int p;
                int guard = 0;
                do { p = nd(rng); ++guard; } while (p == freed && guard < 32);
                if (p == freed) continue;
                fp.bst.insert_at(freed, p, (rng() & 1u), rng);
            }

            finalize();
            if (!is_valid()) continue;
            double s = score();
            if (s < best_s - improve_eps) {
                best_s = s;
                best   = take_snapshot();
                found  = true;
            }
        }

        if (found) {
            restore_snapshot(best);
            current = best_s;
            return true;
        }
        restore_snapshot(start);
        return false;
    }
};
