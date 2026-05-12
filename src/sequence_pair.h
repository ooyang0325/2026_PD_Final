#pragma once
#include "types.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>

// Sequence Pair floorplan representation with O(n^2) evaluation.
// Sequence pair (Gamma+, Gamma-): block i is left of j iff i appears before j in both Gamma+, Gamma-.
// Block i is below j iff i appears before j in Gamma+ but after in Gamma-.
// Evaluation: longest path in horizontal/vertical constraint graph = block positions.

class SequencePair {
public:
    int n;
    std::vector<int> gp, gm;   // Gamma+, Gamma- (permutations of 0..n-1)
    std::vector<double> x, y;  // block lower-left coords after eval

    SequencePair() : n(0) {}
    explicit SequencePair(int n_) : n(n_), gp(n_), gm(n_), x(n_, 0), y(n_, 0) {
        std::iota(gp.begin(), gp.end(), 0);
        std::iota(gm.begin(), gm.end(), 0);
    }

    // Evaluate packing given block widths/heights.
    // Returns (total_width, total_height).
    std::pair<double,double> evaluate(const std::vector<double>& W,
                                       const std::vector<double>& H) {
        // Build inverse of gm
        std::vector<int> inv_gm(n);
        for (int i = 0; i < n; i++) inv_gm[gm[i]] = i;

        // Longest path in H-graph: x[j] = max over all i left-of j of (x[i]+W[i])
        // i is left of j iff pos_gp[i] < pos_gp[j] AND pos_gm[i] < pos_gm[j]
        // Process in Gamma+ order; for each block in gp, check already-placed blocks in gm order.
        // O(n^2) sweep:
        x.assign(n, 0);
        y.assign(n, 0);

        // Position in Gamma+
        std::vector<int> pos_gp(n);
        for (int i = 0; i < n; i++) pos_gp[gp[i]] = i;

        // Horizontal: process blocks in Gamma+ order
        // For block gp[k], all blocks gp[0..k-1] that appear before gp[k] in Gamma- are "left-of"
        for (int k = 0; k < n; k++) {
            int b = gp[k];
            double best_x = 0;
            for (int p = 0; p < k; p++) {
                int a = gp[p];
                // a is before b in Gamma+ (p < k). Check if a before b in Gamma-
                if (inv_gm[a] < inv_gm[b]) {
                    // a is left-of b
                    best_x = std::max(best_x, x[a] + W[a]);
                }
            }
            x[b] = best_x;
        }

        // Vertical: i is below j iff pos_gp[i] < pos_gp[j] AND pos_gm[i] > pos_gm[j]
        // Process in Gamma+ order
        for (int k = 0; k < n; k++) {
            int b = gp[k];
            double best_y = 0;
            for (int p = 0; p < k; p++) {
                int a = gp[p];
                if (inv_gm[a] > inv_gm[b]) {
                    // a is below b
                    best_y = std::max(best_y, y[a] + H[a]);
                }
            }
            y[b] = best_y;
        }

        double tw = 0, th = 0;
        for (int i = 0; i < n; i++) {
            tw = std::max(tw, x[i] + W[i]);
            th = std::max(th, y[i] + H[i]);
        }
        return {tw, th};
    }

    // --- Perturbation moves ---

    // M1: swap two blocks in Gamma+
    void swap_gp(int i, int j) { std::swap(gp[i], gp[j]); }
    // M2: swap two blocks in Gamma-
    void swap_gm(int i, int j) { std::swap(gm[i], gm[j]); }
    // M3: swap same element in both (changes topology without rotation)
    void swap_both(int i, int j) { swap_gp(i,j); swap_gm(i,j); }

    // Random move: returns move type (1,2,3) and indices
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

    // Undo move (swap is its own inverse)
    void undo_move(int mt, int i, int j) {
        if (mt == 1) swap_gp(i, j);
        else if (mt == 2) swap_gm(i, j);
        else swap_both(i, j);
    }

    // Save/restore
    struct State {
        std::vector<int> gp, gm;
    };
    State save() const { return {gp, gm}; }
    void restore(const State& s) { gp = s.gp; gm = s.gm; }
};
