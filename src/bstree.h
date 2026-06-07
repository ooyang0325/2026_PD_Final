#pragma once
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
#include <algorithm>

// ─── B*-Tree floorplan representation ───────────────────────────────────────
// Replaces Sequence Pair.  Each tree NODE holds a BLOCK id (blk[node]).
//   - left child  : block placed immediately to the RIGHT of parent (x = px + pw)
//   - right child : block placed at the SAME x as parent, stacked ABOVE (contour)
// Packing is an O(n) DFS with a skyline contour (here an O(n) scan per node,
// which is trivially fast for n <= 50).  Coordinates x[],y[] are indexed by
// BLOCK id (not node id) so the rest of the floorplanner can read them directly.
//
// Perturbations (classic Chang et al. B*-tree moves):
//   swap_nodes(a,b)  : exchange the block ids held by two nodes        (O(1))
//   move_random()    : delete a node to a leaf and reinsert elsewhere  (O(h))
// Rotation / aspect-ratio changes are handled by the Floorplan (W/H arrays).

class BStarTree {
public:
    int n;
    int root;
    std::vector<int> par, lc, rc;   // topology over node ids 0..n-1
    std::vector<int> blk;           // blk[node] = block id occupying this node
    std::vector<double> x, y;       // packed lower-left, indexed by BLOCK id

    // reusable scratch (avoids per-pack allocation in the SA hot loop)
    std::vector<double> nodex;
    std::vector<int> stk, placed;

    BStarTree() : n(0), root(-1) {}

    explicit BStarTree(int n_)
        : n(n_), root(n_ > 0 ? 0 : -1),
          par(n_, -1), lc(n_, -1), rc(n_, -1), blk(n_),
          x(n_, 0.0), y(n_, 0.0), nodex(n_, 0.0) {
        // Complete binary tree over node ids for a square-ish initial packing.
        std::iota(blk.begin(), blk.end(), 0);
        for (int i = 0; i < n; i++) {
            int l = 2 * i + 1, r = 2 * i + 2;
            lc[i] = (l < n) ? l : -1;
            rc[i] = (r < n) ? r : -1;
            if (l < n) par[l] = i;
            if (r < n) par[r] = i;
        }
        stk.reserve(n);
        placed.reserve(n);
    }

    // Pack with the given block widths/heights (indexed by BLOCK id).
    // Returns (total_width, total_height).  Fills x[],y[].
    std::pair<double,double> pack(const std::vector<double>& W,
                                  const std::vector<double>& H) {
        stk.clear();
        placed.clear();
        if (root < 0) return {0.0, 0.0};

        stk.push_back(root);
        nodex[root] = 0.0;

        while (!stk.empty()) {
            int nd = stk.back(); stk.pop_back();
            int b = blk[nd];
            double bx = nodex[nd];
            double w = W[b];

            // y = skyline height over [bx, bx+w] across already-placed blocks
            double by = 0.0;
            for (int pb : placed) {
                if (x[pb] < bx + w - 1e-9 && x[pb] + W[pb] > bx + 1e-9) {
                    double top = y[pb] + H[pb];
                    if (top > by) by = top;
                }
            }
            x[b] = bx;
            y[b] = by;
            placed.push_back(b);

            // Push right first so left is processed first (node, left, right DFS).
            if (rc[nd] != -1) { nodex[rc[nd]] = bx;      stk.push_back(rc[nd]); }
            if (lc[nd] != -1) { nodex[lc[nd]] = bx + w;  stk.push_back(lc[nd]); }
        }

        double tw = 0, th = 0;
        for (int i = 0; i < n; i++) {
            if (x[i] + W[i] > tw) tw = x[i] + W[i];
            if (y[i] + H[i] > th) th = y[i] + H[i];
        }
        return {tw, th};
    }

    // ── Moves ───────────────────────────────────────────────────────────────
    void swap_nodes(int a, int b) { std::swap(blk[a], blk[b]); }

    // Remove the content originally at node 'a' by migrating it down to a leaf,
    // then detaching that leaf.  Returns the freed (now-detached) node id.
    int delete_to_leaf(int a, std::mt19937& rng) {
        int cur = a;
        while (lc[cur] != -1 || rc[cur] != -1) {
            int c;
            if (rc[cur] == -1)      c = lc[cur];
            else if (lc[cur] == -1) c = rc[cur];
            else                    c = (rng() & 1u) ? lc[cur] : rc[cur];
            std::swap(blk[cur], blk[c]);
            cur = c;
        }
        int p = par[cur];
        if (p != -1) { if (lc[p] == cur) lc[p] = -1; else rc[p] = -1; }
        par[cur] = -1; lc[cur] = -1; rc[cur] = -1;
        return cur;
    }

    // Attach freed node as a child of p on the given side, pushing any existing
    // child down to become a child of the freed node.
    void insert_at(int freed, int p, bool left, std::mt19937& rng) {
        int old = left ? lc[p] : rc[p];
        if (left) lc[p] = freed; else rc[p] = freed;
        par[freed] = p; lc[freed] = -1; rc[freed] = -1;
        if (old != -1) {
            if (rng() & 1u) lc[freed] = old; else rc[freed] = old;
            par[old] = freed;
        }
    }

    void move_random(std::mt19937& rng) {
        if (n < 2) return;
        std::uniform_int_distribution<int> nd(0, n - 1);
        int a = nd(rng);
        int freed = delete_to_leaf(a, rng);   // freed != root for n >= 2
        int p;
        do { p = nd(rng); } while (p == freed);
        insert_at(freed, p, (rng() & 1u), rng);
    }

    // ── Save / restore (selective use in SA) ─────────────────────────────────
    struct State { int root; std::vector<int> par, lc, rc, blk; };
    void save(State& s) const { s.root = root; s.par = par; s.lc = lc; s.rc = rc; s.blk = blk; }
    State save() const { State s; save(s); return s; }
    void restore(const State& s) { root = s.root; par = s.par; lc = s.lc; rc = s.rc; blk = s.blk; }
};
