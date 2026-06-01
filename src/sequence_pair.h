#pragma once
#include "types.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>

// Sequence Pair floorplan representation.
// (Gamma+, Gamma-): block i is left of j iff i appears before j in both.
// Block i is below j iff i appears before j in Gamma+ but after in Gamma-.
//
// Evaluation uses a Fenwick (BIT) prefix-max over inv_gm position, giving
// O(n log n) per pack instead of O(n^2). All scratch buffers are pre-allocated
// to avoid heap churn inside the SA hot loop.

class SequencePair {
public:
    int n;
    std::vector<int> gp, gm;   // Gamma+, Gamma- (permutations of 0..n-1)
    std::vector<double> x, y;  // block lower-left coords after eval

    // Pre-allocated scratch (avoid alloc per evaluate())
    std::vector<int> inv_gm;
    std::vector<double> ft_x, ft_y; // Fenwick trees (1-indexed, size n+1)

    SequencePair() : n(0) {}
    explicit SequencePair(int n_) : n(n_), gp(n_), gm(n_),
        x(n_, 0), y(n_, 0), inv_gm(n_),
        ft_x(n_ + 1, 0.0), ft_y(n_ + 1, 0.0) {
        std::iota(gp.begin(), gp.end(), 0);
        std::iota(gm.begin(), gm.end(), 0);
    }

    // Evaluate packing given block widths/heights.
    // Returns (total_width, total_height).
    std::pair<double,double> evaluate(const std::vector<double>& W,
                                      const std::vector<double>& H) {
        // inv_gm[a] = position of block a in Gamma-
        for (int i = 0; i < n; i++) inv_gm[gm[i]] = i;

        // ----- X sweep -----
        // x[b] = max over a left-of b of (x[a] + W[a])
        // a is left of b iff inv_gm[a] < inv_gm[b] AND a precedes b in gp.
        // Fenwick indexed by 1..n with pos = inv_gm[b]+1.
        std::fill(ft_x.begin(), ft_x.end(), 0.0);
        for (int k = 0; k < n; k++) {
            int b = gp[k];
            int pos = inv_gm[b] + 1;
            // prefix max [1..pos-1]
            double best_x = 0;
            for (int i = pos - 1; i > 0; i -= i & -i)
                if (ft_x[i] > best_x) best_x = ft_x[i];
            x[b] = best_x;
            // update Fenwick at pos with x[b] + W[b]
            double v = best_x + W[b];
            for (int i = pos; i <= n; i += i & -i)
                if (v > ft_x[i]) ft_x[i] = v;
        }

        // ----- Y sweep -----
        // y[b] = max over a below b of (y[a] + H[a])
        // a is below b iff inv_gm[a] > inv_gm[b] AND a precedes b in gp.
        // Reverse-map: r = n - inv_gm[b], then condition becomes r(a) < r(b),
        // and we use the same prefix-max Fenwick.
        std::fill(ft_y.begin(), ft_y.end(), 0.0);
        for (int k = 0; k < n; k++) {
            int b = gp[k];
            int r = n - inv_gm[b];
            double best_y = 0;
            for (int i = r - 1; i > 0; i -= i & -i)
                if (ft_y[i] > best_y) best_y = ft_y[i];
            y[b] = best_y;
            double v = best_y + H[b];
            for (int i = r; i <= n; i += i & -i)
                if (v > ft_y[i]) ft_y[i] = v;
        }

        double tw = 0, th = 0;
        for (int i = 0; i < n; i++) {
            double xr = x[i] + W[i];
            double yt = y[i] + H[i];
            if (xr > tw) tw = xr;
            if (yt > th) th = yt;
        }
        return {tw, th};
    }

    // --- Perturbation moves (O(1)) ---
    void swap_gp(int i, int j) { std::swap(gp[i], gp[j]); }
    void swap_gm(int i, int j) { std::swap(gm[i], gm[j]); }
    void swap_both(int i, int j) { swap_gp(i,j); swap_gm(i,j); }

    int random_move(std::mt19937& rng) {
        std::uniform_int_distribution<int> move_type(1, 3);
        std::uniform_int_distribution<int> idx(0, n-1);
        int mt = move_type(rng);
        int i = idx(rng), j = idx(rng);
        while (i == j) j = idx(rng);
        if (mt == 1) swap_gp(i, j);
        else if (mt == 2) swap_gm(i, j);
        else swap_both(i, j);
        return mt;
    }

    // Save/restore (used only for best-state checkpoints, not per-iteration)
    struct State {
        std::vector<int> gp, gm;
    };
    State save() const { return {gp, gm}; }
    void restore(const State& s) { gp = s.gp; gm = s.gm; }
};
