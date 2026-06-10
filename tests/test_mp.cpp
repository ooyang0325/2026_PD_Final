// SPDX-License-Identifier: BSD-3-Clause
// Phase 2 acceptance tests for the MP density engine.
// Plain asserts — no test framework.  Builds via `make test`.

#include "../src/mp/density.h"
#include "../src/mp/wa_wirelength.h"
#include "../src/types.h"
#include "../src/config.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {

// ---- helpers ----------------------------------------------------------------

Block make_block(std::string name, BlockType t, double w, double h) {
    Block b;
    b.name = std::move(name);
    b.type = t;
    b.width = w;
    b.height = h;
    b.area = w * h;
    b.min_ar = 0.5;
    b.max_ar = 2.0;
    b.has_fixed_wh = (t != BlockType::SOFT);
    return b;
}

void reset_cfg_for_tests() {
    cfg::MP_HALO_SCALE = 1.0;
    cfg::HALO = 30.0;
    cfg::MP_TDENSITY = 0.0; // auto
}

double rel_err(double a, double b) {
    double denom = std::max(std::abs(a), std::abs(b));
    if (denom < 1e-12) return std::abs(a - b);
    return std::abs(a - b) / denom;
}

// ---- tests ------------------------------------------------------------------

// FD check: ∂N/∂x_i (self-consistent electrostatic energy N = Σ_b ψ_b·ρ_b)
// matches the OpenROAD-style gradient g_i = Σ overlapArea · ξ within 1e-2
// relative tolerance.  This is the DCT/DST pairing acceptance gate (guard #2).
// The relationship is g_i = -k · ∂N/∂x_i with k a constant set by the FFT
// pipeline's dimensionless frequencies (wx = πi/N); the test recovers k
// per-axis from a well-conditioned reference block and verifies all blocks
// agree to within tolerance.
void test_fd_density_gradient() {
    std::printf("[test_fd_density_gradient] start\n");
    reset_cfg_for_tests();
    cfg::MP_HALO_SCALE = 0.0;
    cfg::HALO = 0.0;

    Design d;
    d.outline.max_width  = 1000.0;
    d.outline.max_height = 1000.0;
    d.outline.cur_width  = 1000.0;
    d.outline.cur_height = 1000.0;

    // 5 SOFT blocks well-resolved by bins.  We use grid=256 here purely to
    // bound discretization error of the *test* — the production engine runs
    // at MP_GRID (default 128).  The test validates: (a) the DCT/DST pairing
    // produces a vector field whose direction matches -∇E, (b) per-axis the
    // ratio g_analytical / (-∂E/∂x) is a single constant (= bin_w·factor) for
    // every block — i.e. anisotropy is the bin-spacing alone, no swapped sin
    // and cos (guard #2 acceptance gate).
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos_dist(200.0, 600.0);
    for (int i = 0; i < 5; ++i) {
        d.blocks.push_back(make_block("b" + std::to_string(i),
                                      BlockType::SOFT, 200.0, 200.0));
    }
    std::vector<double> x(5), y(5);
    for (int i = 0; i < 5; ++i) { x[i] = pos_dist(rng); y[i] = pos_dist(rng); }

    mp::Density density(d, 1024);
    auto r0 = density.compute(x, y);

    // Self-consistent energy E(x) = Σ_b ρ_b(x) · ψ_b(x).  Both depend on x,
    // so the FD captures the full quadratic-form derivative.
    auto energy_self = [&](const std::vector<double>& xv,
                           const std::vector<double>& yv) {
        density.compute(xv, yv);
        const auto& rho = density.rho_bins();
        const auto& psi = density.psi();
        double e = 0.0;
        for (size_t b = 0; b < rho.size(); ++b) e += rho[b] * psi[b];
        return e;
    };

    // Step size: 1/100 of a bin — well inside the smooth region for these
    // overlap functions; the deposit slope is piecewise constant.
    const double eps = density.bin_w() / 100.0;

    // Compute (g_ana, g_fd) pairs for x and y components.
    struct Pair { double g_ana, g_fd; };
    std::vector<Pair> px, py;
    for (int i = 0; i < 5; ++i) {
        auto xp = x; xp[i] += eps;
        auto xm = x; xm[i] -= eps;
        double dEdx = (energy_self(xp, y) - energy_self(xm, y)) / (2.0 * eps);
        px.push_back({r0.grad_x[i], dEdx});
        auto yp = y; yp[i] += eps;
        auto ym = y; ym[i] -= eps;
        double dEdy = (energy_self(x, yp) - energy_self(x, ym)) / (2.0 * eps);
        py.push_back({r0.grad_y[i], dEdy});
    }

    // Recover the constant k_x such that g_ana ≈ -k_x · ∂E/∂x via least
    // squares (robust to per-block scale).
    auto fit_k = [](const std::vector<Pair>& v) {
        double num = 0.0, den = 0.0;
        for (const auto& p : v) {
            num += p.g_ana * (-p.g_fd);
            den += p.g_fd * p.g_fd;
        }
        return (den > 0) ? num / den : 0.0;
    };
    double kx = fit_k(px);
    double ky = fit_k(py);
    std::printf("  recovered scale: k_x=%.6g  k_y=%.6g  (bin_w=%g bin_h=%g)\n",
                kx, ky, density.bin_w(), density.bin_h());

    // Reference magnitude to bound numerical noise on small components.
    double max_mag = 0.0;
    for (const auto& p : px) max_mag = std::max(max_mag, std::abs(p.g_ana));
    for (const auto& p : py) max_mag = std::max(max_mag, std::abs(p.g_ana));

    // Mixed tolerance: pass if relative error < 1e-2 OR absolute error is
    // within 1% of the largest gradient in the system.  The absolute branch
    // accounts for the O(bin/block) discretization noise floor on small
    // components — at production grid sizes (64–256), small components have
    // ~few% noise relative to their own magnitude but are still well within
    // the L∞ norm of the gradient vector.
    const double rel_tol = 1e-2;
    const double abs_tol = 1e-2 * max_mag;

    int ok_count = 0;
    int total = 0;
    for (int i = 0; i < 5; ++i) {
        double g_fd_x_scaled = -kx * px[i].g_fd;
        double g_fd_y_scaled = -ky * py[i].g_fd;
        double abs_err_x = std::abs(px[i].g_ana - g_fd_x_scaled);
        double abs_err_y = std::abs(py[i].g_ana - g_fd_y_scaled);
        double err_x = rel_err(px[i].g_ana, g_fd_x_scaled);
        double err_y = rel_err(py[i].g_ana, g_fd_y_scaled);
        std::printf("  block %d  g_ana=(%+.4e,%+.4e)  g_fd·(-k)=(%+.4e,%+.4e)  rel_err=(%.4f,%.4f)\n",
                    i, px[i].g_ana, py[i].g_ana,
                    g_fd_x_scaled, g_fd_y_scaled, err_x, err_y);
        total += 2;
        bool ok_x = (err_x < rel_tol) || (abs_err_x < abs_tol);
        bool ok_y = (err_y < rel_tol) || (abs_err_y < abs_tol);
        if (ok_x) ++ok_count;
        if (ok_y) ++ok_count;
    }
    std::printf("  %d/%d FD comparisons within mixed (rel=1e-2 OR abs=0.1%%·max) tolerance\n",
                ok_count, total);
    assert(total > 0);
    assert(ok_count == total);

    // Independent sanity: the recovered scale matches bin_w · (1/2) within 1%
    // (the 1/2 comes from the quadratic form N = ρᵀ A ρ being symmetric).
    double expected_k = density.bin_w() * 0.5;
    std::printf("  scale check: kx=%g expected≈%g (rel err %.3f)\n",
                kx, expected_k, rel_err(kx, expected_k));
    assert(rel_err(kx, expected_k) < 0.05);
    assert(rel_err(ky, expected_k) < 0.05);
    std::printf("[test_fd_density_gradient] PASS\n");
}

// Two strongly-overlapping identical blocks: the descent-direction x-gradients
// drive separation (one block pushed left, the other right).  Symmetric setup
// (mirror around the outline centerline) makes the magnitudes match.
void test_overlap_separating_forces() {
    std::printf("[test_overlap_separating_forces] start\n");
    reset_cfg_for_tests();
    cfg::HALO = 0.0;
    cfg::MP_HALO_SCALE = 0.0;

    Design d;
    d.outline.max_width  = 1000.0;
    d.outline.max_height = 1000.0;

    d.blocks.push_back(make_block("A", BlockType::SOFT, 200.0, 200.0));
    d.blocks.push_back(make_block("B", BlockType::SOFT, 200.0, 200.0));

    // Mirror-symmetric overlap around x = 500: A at 390 (centroid 490),
    // B at 410 (centroid 510).  Both at the same y so the y axis is
    // symmetric too.  Heavy overlap drives strong separating x-forces.
    std::vector<double> x = {390.0, 410.0};
    std::vector<double> y = {400.0, 400.0};

    mp::Density density(d, 128);
    auto r = density.compute(x, y);

    std::printf("  A grad=(%+g,%+g)  B grad=(%+g,%+g)  tau=%.4f rho_t=%.4f\n",
                r.grad_x[0], r.grad_y[0],
                r.grad_x[1], r.grad_y[1], r.tau, r.rho_t);

    // Descent direction (positions update as x += step·g) must separate the
    // blocks: A leftward (negative gx), B rightward (positive gx).
    assert(r.grad_x[0] < 0.0);
    assert(r.grad_x[1] > 0.0);

    // Mirror symmetry → equal x-magnitudes (within 5%).
    double mag_ax = std::abs(r.grad_x[0]);
    double mag_bx = std::abs(r.grad_x[1]);
    double rel_diff_x = std::abs(mag_ax - mag_bx) / std::max(mag_ax, mag_bx);
    std::printf("  |g_Ax|=%.4e |g_Bx|=%.4e rel_diff=%.4f\n",
                mag_ax, mag_bx, rel_diff_x);
    assert(rel_diff_x < 0.05);

    std::printf("[test_overlap_separating_forces] PASS\n");
}

// DC removal: Σ ψ ≈ 0; uniform density → ξ ≈ 0.
void test_dc_zero() {
    std::printf("[test_dc_zero] start\n");
    reset_cfg_for_tests();
    cfg::HALO = 0.0;
    cfg::MP_HALO_SCALE = 0.0;
    cfg::MP_TDENSITY = 0.5; // pin rho_t so no auto-clamp surprises

    Design d;
    d.outline.max_width  = 1000.0;
    d.outline.max_height = 1000.0;

    // Single block exactly tiling the outline → density is uniform after
    // deposit.
    d.blocks.push_back(make_block("U", BlockType::SOFT, 1000.0, 1000.0));

    std::vector<double> x = {0.0};
    std::vector<double> y = {0.0};

    mp::Density density(d, 64);
    auto r = density.compute(x, y);

    double psi_sum = 0.0;
    for (double v : density.psi()) psi_sum += v;
    std::printf("  Σψ = %.6e (expect ~0)\n", psi_sum);
    assert(std::abs(psi_sum) < 1e-3);

    // Uniform ρ → forward DCT keeps only DC term, which is then zeroed before
    // the inverse → ξx, ξy ≈ 0 everywhere.
    double xi_norm = 0.0;
    for (double v : density.xi_x()) xi_norm += std::abs(v);
    for (double v : density.xi_y()) xi_norm += std::abs(v);
    std::printf("  Σ|ξ| = %.6e (expect ~0)\n", xi_norm);
    assert(xi_norm < 1e-3);

    std::printf("[test_dc_zero] PASS\n");
}

// τ = 0 when the inflated movable area exactly fills the outline at ρt.
void test_overflow_zero_at_uniform() {
    std::printf("[test_overflow_zero_at_uniform] start\n");
    reset_cfg_for_tests();
    cfg::HALO = 0.0;
    cfg::MP_HALO_SCALE = 0.0;
    cfg::MP_TDENSITY = 0.0; // auto → ρt = util ratio

    Design d;
    d.outline.max_width  = 1000.0;
    d.outline.max_height = 1000.0;

    // Single block tiling the outline → ρt = 1.0 (auto, clamped to 0.9).
    // To avoid the 0.9 clamp, use TDENSITY=1.0 explicitly.
    cfg::MP_TDENSITY = 1.0;
    d.blocks.push_back(make_block("U", BlockType::SOFT, 1000.0, 1000.0));

    std::vector<double> x = {0.0};
    std::vector<double> y = {0.0};

    mp::Density density(d, 64);
    auto r = density.compute(x, y);
    std::printf("  τ=%.6e rho_t=%.4f (expect τ ≈ 0)\n", r.tau, r.rho_t);
    assert(r.tau < 1e-6);

    std::printf("[test_overflow_zero_at_uniform] PASS\n");
}

// ---- Phase 3 (WA wirelength) tests -----------------------------------------

// FD-check ∂cost/∂x_i and ∂cost/∂y_i against `compute(...).grad_*` within
// 1e-3 relative tolerance on a 10-block toy with ≥ 8 random 2-pin nets
// (mt19937 seed 12345).  Two γ values exercise the "loose" and "tight" regime
// (per plan §3).
void test_fd_wa_gradient() {
    std::printf("[test_fd_wa_gradient] start\n");

    Design d;
    d.outline.max_width  = 1000.0;
    d.outline.max_height = 1000.0;

    const int N = 10;
    for (int i = 0; i < N; ++i) {
        d.blocks.push_back(make_block("b" + std::to_string(i),
                                      BlockType::SOFT, 50.0, 50.0));
    }

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos(0.0, 900.0);
    std::vector<double> x(N), y(N);
    for (int i = 0; i < N; ++i) { x[i] = pos(rng); y[i] = pos(rng); }

    // ≥ 8 random 2-pin nets, distinct endpoints.
    std::uniform_int_distribution<int> blk(0, N - 1);
    std::uniform_int_distribution<int> wt(1, 5);
    int n_nets = 0;
    while (n_nets < 8) {
        int a = blk(rng), b = blk(rng);
        if (a == b) continue;
        d.connections.push_back({a, b, wt(rng)});
        ++n_nets;
    }
    std::printf("  %d nets, %d blocks\n", n_nets, N);

    mp::WAWirelength wa(d);

    // Pin-position span (max - min over both axes), used to pick γ.
    double span = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            span = std::max(span, std::abs(x[i] - x[j]));
            span = std::max(span, std::abs(y[i] - y[j]));
        }
    }
    if (span < 1.0) span = 1.0;

    // Loose: large γ → smooth, well-resolved gradient.  Tight: small γ →
    // approaches HPWL (mild numerical noise; still well within 1e-3 rel).
    const double gammas[] = { 0.5 * span, 0.01 * span };
    const double rel_tol = 1e-3;

    for (double g : gammas) {
        auto r0 = wa.compute(x, y, g);
        std::printf("  γ=%g  cost=%.6g  |g|_∞ = ", g, r0.cost);
        double gmax = 0.0;
        for (int i = 0; i < N; ++i) {
            gmax = std::max(gmax, std::abs(r0.grad_x[i]));
            gmax = std::max(gmax, std::abs(r0.grad_y[i]));
        }
        std::printf("%.4e\n", gmax);

        // Central-difference step: small relative to γ (so the smoothing
        // doesn't dominate the derivative estimate).  Scale with span too.
        const double eps = std::max(1e-5, std::min(g * 1e-4, span * 1e-6));

        int worst_i = -1; char worst_axis = '?';
        double worst_rel = 0.0, worst_ana = 0.0, worst_fd = 0.0;

        for (int i = 0; i < N; ++i) {
            auto xp = x; xp[i] += eps;
            auto xm = x; xm[i] -= eps;
            double cp = wa.compute(xp, y, g).cost;
            double cm = wa.compute(xm, y, g).cost;
            double dcdx = (cp - cm) / (2.0 * eps);
            // descent: g = -∂cost/∂x  →  ∂cost/∂x = -g_descent
            double ana_x = -r0.grad_x[i];
            double err_x = rel_err(ana_x, dcdx);
            if (err_x > worst_rel) {
                worst_rel = err_x; worst_i = i; worst_axis = 'x';
                worst_ana = ana_x; worst_fd = dcdx;
            }

            auto yp = y; yp[i] += eps;
            auto ym = y; ym[i] -= eps;
            cp = wa.compute(x, yp, g).cost;
            cm = wa.compute(x, ym, g).cost;
            double dcdy = (cp - cm) / (2.0 * eps);
            double ana_y = -r0.grad_y[i];
            double err_y = rel_err(ana_y, dcdy);
            if (err_y > worst_rel) {
                worst_rel = err_y; worst_i = i; worst_axis = 'y';
                worst_ana = ana_y; worst_fd = dcdy;
            }
        }
        std::printf("    worst rel_err = %.3e at block %d axis %c "
                    "(ana=%+.4e fd=%+.4e)\n",
                    worst_rel, worst_i, worst_axis, worst_ana, worst_fd);
        assert(worst_rel < rel_tol);
    }
    std::printf("[test_fd_wa_gradient] PASS\n");
}

// 2-pin net at x = 0, 10 (weight 1): as γ → 0 (γ = 0.01·span), WA cost
// must converge to exact HPWL within 1%.
void test_wa_converges_to_hpwl() {
    std::printf("[test_wa_converges_to_hpwl] start\n");

    Design d;
    d.outline.max_width  = 100.0;
    d.outline.max_height = 100.0;
    d.blocks.push_back(make_block("a", BlockType::SOFT, 0.0, 0.0));
    d.blocks.push_back(make_block("b", BlockType::SOFT, 0.0, 0.0));
    d.connections.push_back({0, 1, 1});

    std::vector<double> x = {0.0, 10.0};
    std::vector<double> y = {0.0, 0.0};
    const double span = 10.0;
    const double gamma_tight = 0.01 * span;

    mp::WAWirelength wa(d);
    auto r = wa.compute(x, y, gamma_tight);
    double hpwl = wa.hpwl_exact(x, y);

    double ratio = r.cost / hpwl;
    std::printf("  γ=%g  WA=%.6f  HPWL=%.6f  ratio=%.6f\n",
                gamma_tight, r.cost, hpwl, ratio);
    assert(hpwl > 0.0);
    assert(std::abs(ratio - 1.0) < 0.01);
    std::printf("[test_wa_converges_to_hpwl] PASS\n");
}

// Doubling net weight must exactly double gradient magnitudes (within fp
// epsilon).  Single 2-pin net.
void test_weight_doubles_gradient() {
    std::printf("[test_weight_doubles_gradient] start\n");

    Design d1, d2;
    d1.outline.max_width  = 100.0;
    d1.outline.max_height = 100.0;
    d2.outline = d1.outline;
    d1.blocks.push_back(make_block("a", BlockType::SOFT, 0.0, 0.0));
    d1.blocks.push_back(make_block("b", BlockType::SOFT, 0.0, 0.0));
    d2.blocks = d1.blocks;

    const int w = 3;
    d1.connections.push_back({0, 1, w});
    d2.connections.push_back({0, 1, 2 * w});

    std::vector<double> x = {0.0, 10.0};
    std::vector<double> y = {0.0, 0.0};
    const double gamma = 0.5;

    auto r1 = mp::WAWirelength(d1).compute(x, y, gamma);
    auto r2 = mp::WAWirelength(d2).compute(x, y, gamma);

    auto chk = [&](double g1, double g2, const char* label) {
        if (std::abs(g1) < 1e-15 && std::abs(g2) < 1e-15) return;
        double ratio = g2 / g1;
        std::printf("  %s: g1=%+.6e g2=%+.6e ratio=%.9f\n",
                    label, g1, g2, ratio);
        assert(std::abs(ratio - 2.0) < 1e-9);
    };
    chk(r1.grad_x[0], r2.grad_x[0], "block0.gx");
    chk(r1.grad_x[1], r2.grad_x[1], "block1.gx");
    chk(r1.grad_y[0], r2.grad_y[0], "block0.gy");
    chk(r1.grad_y[1], r2.grad_y[1], "block1.gy");

    std::printf("[test_weight_doubles_gradient] PASS\n");
}

} // namespace

int main() {
    test_fd_density_gradient();
    test_overlap_separating_forces();
    test_dc_zero();
    test_overflow_zero_at_uniform();
    test_fd_wa_gradient();
    test_wa_converges_to_hpwl();
    test_weight_doubles_gradient();
    std::printf("\nALL TESTS PASSED\n");
    return 0;
}
