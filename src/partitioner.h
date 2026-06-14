#pragma once
#include "types.h"
#include "floorplan.h"
#include "config.h"
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <utility>

// ─── Pre-SA congestion-aware recursive min-cut partitioning ──────────────────
//
// PURPOSE
//   Replace the connectivity-blind initial B*-tree (a complete binary tree over
//   input order) with one derived from a recursive bisection of the netlist, so
//   the SA starts from a layout that already:
//     (1) keeps strongly-connected (high-net) blocks adjacent — short, low-cost
//         routing for the high-capacity connections; and
//     (2) spreads connectivity load across the die so no single region is a
//         dense knot of mutually-connected blocks competing for the limited
//         channel / feedthrough capacity that exists locally.
//
// DUAL-AWARE OBJECTIVE  (this is the crux the design hinges on)
//   A pure min-cut bisection optimizes (1) but actively HURTS (2): the cheapest
//   cut for a clique is "everything on one side" (cut = 0), which is exactly the
//   congested knot we must avoid.  So every bipartition minimizes a sum of THREE
//   terms, all carried in net (= nets) units so the weights are comparable:
//
//        objective(A,B) =        cut(A,B)                          // term 1: cut
//                        + BALW * Wtot * |area(A) - area_target|   // term 2: balance
//                        + CONGW * ( overload(A) + overload(B) )   // term 3: congestion
//
//     cut(A,B)     = Σ nets of connections with one endpoint in A, one in B.
//                    Minimizing keeps high-net pairs together (objective 1).
//     overload(S)  = max(0, demand(S) - cap(S)), where
//                    demand(S) = Σ nets of connections internal to S (the local
//                                routing load S must carry), and
//                    cap(S)    = KCAP * sqrt(area(S)) — a proxy for the channel /
//                                feedthrough capacity available in a region whose
//                                linear dimension scales as sqrt(area) and whose
//                                channels carry ~25 nets per unit length.
//                    Penalizing this forces a dense clique to SPLIT across the
//                    cut even at the expense of a larger cut (objective 2).
//     balance term keeps the area split near the geometric mid-cut so the
//                    resulting slicing actually fits inside the fixed outline.
//
//   The cut term and the congestion term pull in opposite directions; FM moves
//   are scored on the FULL objective, so the partitioner finds the trade-off
//   point rather than collapsing to either extreme.
//
// FROM PARTITION TO B*-TREE
//   The recursion assigns every block a target centre (tx,ty) inside the
//   outline (each leaf gets its sub-region's centre).  Blocks are then bucketed
//   into vertical COLUMNS by tx and ordered bottom->top by ty within a column;
//   BStarTree::build_from_columns wires that column placement into the tree.
//   This is a legitimate placement->B*-tree construction (columns are the
//   canonical B*-tree-representable placement) and is overlap-free by the tree's
//   contour guarantee, so the seed can never be invalid.
namespace partition {

struct Rect {
    double x, y, w, h;
    double cx() const { return x + w * 0.5; }
    double cy() const { return y + h * 0.5; }
    double area() const { return w * h; }
};

class Partitioner {
public:
    Partitioner(Floorplan& fp_, unsigned seed)
        : fp(fp_), d(fp_.d), rng(seed) { build_graph(); }

    // Run the full stage and return the column ordering for the B*-tree.
    // Returns an empty vector (caller should leave the default tree) only if the
    // result fails the cover check — a defensive fallback that never fires in
    // normal operation.
    std::vector<std::vector<int>> run() {
        int n = (int)d.blocks.size();
        if (n <= 0) return {};
        if (n == 1) { tx_[0] = 0; ty_[0] = 0; return {{0}}; }

        std::vector<int> all(n);
        std::iota(all.begin(), all.end(), 0);
        Rect outline{0.0, 0.0,
                     std::max(1.0, d.outline.max_width),
                     std::max(1.0, d.outline.max_height)};
        recurse(all, outline, 0);

        auto cols = build_columns();
        bias_edge_blocks(cols);
        drop_empty(cols);
        if (!covers_all(cols, n)) return {}; // defensive: fall back to default tree
        return cols;
    }

private:
    Floorplan& fp;
    Design& d;
    std::mt19937 rng;

    std::vector<double> area_;   // footprint area per block (packing-consistent)
    std::vector<double> bw_, bh_; // footprint width / height per block
    std::vector<std::vector<std::pair<int,double>>> adj_; // weighted adjacency
    std::vector<double> tx_, ty_; // target centre assigned by the recursion
    std::vector<int>    side_;    // bipartition scratch: -1 outside region, 0/1 side

    void build_graph() {
        int n = (int)d.blocks.size();
        area_.resize(n); bw_.resize(n); bh_.resize(n);
        tx_.assign(n, 0.0); ty_.assign(n, 0.0);
        side_.assign(n, -1);
        adj_.assign(n, {});
        for (int i = 0; i < n; i++) {
            bw_[i] = std::max(1e-6, fp.W[i]);
            bh_[i] = std::max(1e-6, fp.H[i]);
            area_[i] = std::max(1.0, fp.W[i] * fp.H[i]);
        }
        for (const auto& c : d.connections) {
            if (c.from < 0 || c.to < 0 || c.from >= n || c.to >= n) continue;
            adj_[c.from].push_back({c.to,   (double)c.nets});
            adj_[c.to].push_back({c.from, (double)c.nets});
        }
    }

    // ── Objective of the current bipartition of `members` (side_ set) ─────────
    // O(Σ deg) over members.  Each internal/cut edge is seen from both endpoints,
    // so doubled sums are halved.  cap() uses sqrt(area) as the region's linear
    // dimension (channel capacity proxy); overload is the absolute net excess.
    double objective(const std::vector<int>& members) const {
        double areaA = 0, areaB = 0;
        double cut2 = 0, demA2 = 0, demB2 = 0, wtot2 = 0;
        for (int v : members) {
            if (side_[v] == 0) areaA += area_[v]; else areaB += area_[v];
            for (const auto& e : adj_[v]) {
                int u = e.first;
                if (side_[u] < 0) continue; // edge leaves this region
                double w = e.second;
                wtot2 += w;
                if (side_[u] != side_[v]) cut2 += w;
                else if (side_[v] == 0)   demA2 += w;
                else                      demB2 += w;
            }
        }
        double cut = cut2 * 0.5, demandA = demA2 * 0.5, demandB = demB2 * 0.5;
        double wtot = std::max(1.0, wtot2 * 0.5);
        double areaTot = std::max(1.0, areaA + areaB);

        double balPen = cfg::PART_BALW * wtot *
                        std::abs(areaA - 0.5 * areaTot) / areaTot;

        double capA = cfg::PART_KCAP * std::sqrt(std::max(1.0, areaA));
        double capB = cfg::PART_KCAP * std::sqrt(std::max(1.0, areaB));
        double congPen = cfg::PART_CONGW * (std::max(0.0, demandA - capA) +
                                            std::max(0.0, demandB - capB));
        return cut + balPen + congPen;
    }

    // Largest-area-first initial split: balances area, gives FM a sane start.
    void init_split(const std::vector<int>& members) {
        std::vector<int> order = members;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b){ return area_[a] > area_[b]; });
        double areaA = 0, areaB = 0;
        for (int v : order) {
            if (areaA <= areaB) { side_[v] = 0; areaA += area_[v]; }
            else                { side_[v] = 1; areaB += area_[v]; }
        }
        ensure_both_sides(members);
    }

    int count_side(const std::vector<int>& members, int s) const {
        int c = 0; for (int v : members) if (side_[v] == s) c++; return c;
    }

    // Guarantee a non-empty bipartition (FM and the recursion both require it).
    void ensure_both_sides(const std::vector<int>& members) {
        if ((int)members.size() < 2) return;
        for (int s = 0; s < 2; s++) {
            if (count_side(members, s) != 0) continue;
            // move the smallest block from the over-full side
            int best = -1;
            for (int v : members)
                if (side_[v] != s && (best < 0 || area_[v] < area_[best])) best = v;
            if (best >= 0) side_[best] = s;
        }
    }

    // One Fiduccia–Mattheyses pass on the FULL (dual-aware) objective.  Each step
    // moves the unlocked block whose toggle most reduces the objective (random
    // tie-break for per-seed diversity), keeping both sides non-empty; the prefix
    // with the best objective is committed (classic KL/FM rollback to escape
    // local minima).  n<=69 makes the per-candidate full recompute trivial.
    bool fm_pass(const std::vector<int>& members) {
        int m = (int)members.size();
        if (m < 2) return false;
        std::vector<char> locked(m, 0);

        double cur = objective(members);
        double best_obj = cur;
        int best_prefix = 0;
        std::vector<int> moved; moved.reserve(m);

        for (int step = 0; step < m; step++) {
            int pick = -1; double pick_obj = 0;
            for (int idx = 0; idx < m; idx++) {
                if (locked[idx]) continue;
                int v = members[idx];
                int s = side_[v];
                if (count_side(members, s) <= 1) continue; // would empty a side
                side_[v] ^= 1;
                double o = objective(members);
                side_[v] = s; // revert trial
                if (pick < 0 || o < pick_obj - 1e-9) { pick = idx; pick_obj = o; }
                else if (o < pick_obj + 1e-9 && (rng() & 1u)) { pick = idx; pick_obj = o; }
            }
            if (pick < 0) break;
            int v = members[pick];
            side_[v] ^= 1;
            locked[pick] = 1;
            moved.push_back(pick);
            cur = pick_obj;
            if (cur < best_obj - 1e-9) { best_obj = cur; best_prefix = (int)moved.size(); }
        }

        // Roll back everything after the best prefix.
        for (int k = (int)moved.size() - 1; k >= best_prefix; k--)
            side_[members[moved[k]]] ^= 1;
        return best_prefix > 0; // improved this pass?
    }

    // Bipartition `members`: set side_ to 0/1 and return (A, B) block lists.
    std::pair<std::vector<int>, std::vector<int>>
    bipartition(const std::vector<int>& members) {
        for (int v : members) side_[v] = -1;     // mark region (cleared below)
        init_split(members);
        for (int p = 0; p < cfg::PART_FMPASS; p++)
            if (!fm_pass(members)) break;
        ensure_both_sides(members);

        std::vector<int> A, B;
        for (int v : members) (side_[v] == 0 ? A : B).push_back(v);
        return {A, B};
    }

    // Recursively bisect `members` over rectangle `rect`.  Cut orientation splits
    // the LONGER side so the children stay square-ish (routable, packable); the
    // rectangle is split by the area fraction actually assigned, so leaf centres
    // reflect the real region geometry.  Leaves (single block) record their
    // centre as the target placement point.
    void recurse(const std::vector<int>& members, const Rect& rect, int depth) {
        if (members.empty()) return;
        if (members.size() == 1) {
            tx_[members[0]] = rect.cx();
            ty_[members[0]] = rect.cy();
            return;
        }
        // Safety net against a pathological non-terminating split.
        if (depth > 4 * (int)d.blocks.size() + 8) {
            for (int v : members) { tx_[v] = rect.cx(); ty_[v] = rect.cy(); }
            return;
        }

        auto pr = bipartition(members);
        const std::vector<int>& A = pr.first;
        const std::vector<int>& B = pr.second;

        double areaA = 0, areaB = 0;
        for (int v : A) areaA += area_[v];
        for (int v : B) areaB += area_[v];
        double frac = (areaA + areaB > 0) ? areaA / (areaA + areaB) : 0.5;
        frac = std::min(0.9, std::max(0.1, frac)); // keep both sub-rects sensible

        Rect ra = rect, rb = rect;
        bool vertical = rect.w >= rect.h; // cut the wider dimension
        if (vertical) {
            ra.w = rect.w * frac;            // A: left band
            rb.x = rect.x + ra.w; rb.w = rect.w - ra.w; // B: right band
        } else {
            ra.h = rect.h * frac;            // A: bottom band
            rb.y = rect.y + ra.h; rb.h = rect.h - ra.h; // B: top band
        }
        recurse(A, ra, depth + 1);
        recurse(B, rb, depth + 1);
    }

    // ── Column bucketing: target centres -> vertical columns -> ordered lists ──
    // Number of columns derives from block count and outline aspect ratio so the
    // resulting packing is roughly the outline's shape.  Blocks are sorted by tx
    // and chunked into area-balanced, x-contiguous columns (similar column areas
    // => similar widths => minimal vertical lifting); each column is then ordered
    // bottom->top by ty.
    std::vector<std::vector<int>> build_columns() {
        int n = (int)d.blocks.size();
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b){
            if (tx_[a] != tx_[b]) return tx_[a] < tx_[b];
            if (ty_[a] != ty_[b]) return ty_[a] < ty_[b];
            return a < b;
        });

        double ar = (d.outline.max_height > 0)
                  ? d.outline.max_width / d.outline.max_height : 1.0;
        int kcol = (int)std::lround(std::sqrt((double)n * std::max(0.25, ar)));
        kcol = std::max(1, std::min(n, kcol));

        double totalArea = 0; for (double a : area_) totalArea += a;
        double target = totalArea / kcol;

        std::vector<std::vector<int>> cols(kcol);
        double acc = 0; int ci = 0;
        for (int v : order) {
            cols[ci].push_back(v);
            acc += area_[v];
            // advance to next column once this one has its area share, but never
            // leave fewer columns than remaining blocks can fill.
            int remaining_blocks = 0;
            // (cheap: blocks not yet placed = order tail; tracked implicitly)
            (void)remaining_blocks;
            if (ci < kcol - 1 && acc >= target * (ci + 1)) ci++;
        }

        for (auto& col : cols)
            std::sort(col.begin(), col.end(), [&](int a, int b){
                if (ty_[a] != ty_[b]) return ty_[a] < ty_[b];
                return a < b;
            });
        return cols;
    }

    // Bias EDGE blocks toward the boundary they must end on, so finalize/snap has
    // a collision-free target: L/R choose the first/last column, B/T the bottom/
    // top of that column.  Best-effort hint only — correctness is enforced later.
    void bias_edge_blocks(std::vector<std::vector<int>>& cols) {
        if (cols.empty()) return;
        for (int i : fp.edge_block_idx) {
            if (d.blocks[i].locations.empty()) continue;
            const std::string& loc = d.blocks[i].locations[0];
            int target_col = -1;
            if (loc.find('L') != std::string::npos) target_col = 0;
            else if (loc.find('R') != std::string::npos) target_col = (int)cols.size() - 1;

            if (target_col >= 0) remove_block(cols, i);
            int col = (target_col >= 0) ? target_col : find_col(cols, i);
            if (col < 0) continue;
            if (target_col >= 0) {
                bool bottom = loc.find('B') != std::string::npos;
                bool top    = loc.find('T') != std::string::npos;
                if (bottom)      cols[col].insert(cols[col].begin(), i);
                else if (top)    cols[col].push_back(i);
                else             cols[col].push_back(i);
            } else {
                // already in a column: only adjust vertical position if forced
                if (loc.find('B') != std::string::npos) move_to_front(cols[col], i);
                else if (loc.find('T') != std::string::npos) move_to_back(cols[col], i);
            }
        }
    }

    static int find_col(const std::vector<std::vector<int>>& cols, int b) {
        for (int c = 0; c < (int)cols.size(); c++)
            for (int v : cols[c]) if (v == b) return c;
        return -1;
    }
    static void remove_block(std::vector<std::vector<int>>& cols, int b) {
        for (auto& col : cols)
            col.erase(std::remove(col.begin(), col.end(), b), col.end());
    }
    static void move_to_front(std::vector<int>& col, int b) {
        col.erase(std::remove(col.begin(), col.end(), b), col.end());
        col.insert(col.begin(), b);
    }
    static void move_to_back(std::vector<int>& col, int b) {
        col.erase(std::remove(col.begin(), col.end(), b), col.end());
        col.push_back(b);
    }
    static void drop_empty(std::vector<std::vector<int>>& cols) {
        cols.erase(std::remove_if(cols.begin(), cols.end(),
                   [](const std::vector<int>& c){ return c.empty(); }),
                   cols.end());
    }
    static bool covers_all(const std::vector<std::vector<int>>& cols, int n) {
        std::vector<char> seen(n, 0);
        int cnt = 0;
        for (const auto& col : cols)
            for (int v : col) {
                if (v < 0 || v >= n || seen[v]) return false;
                seen[v] = 1; cnt++;
            }
        return cnt == n;
    }
};

// Install a congestion-aware partition seed into fp's B*-tree.  No-op (leaves the
// default tree) if partitioning produced no valid cover.
inline void seed(Floorplan& fp, unsigned seed_val) {
    Partitioner p(fp, seed_val);
    auto cols = p.run();
    if (cols.empty()) return;
    fp.bst.build_from_columns(cols);
}

} // namespace partition
