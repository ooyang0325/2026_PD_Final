// SPDX-License-Identifier: BSD-3-Clause
// Phase 3: weighted-average (WA) wirelength + gradient for 2-pin nets.
//
// SIGN CONVENTION (single source of truth — guard #1, matches density.h):
//   grad_x[i], grad_y[i] are *descent* directions.  Callers update positions
//   as `x[i] += α·g[i]` (no sign flip).  OpenROAD's
//   `getWireLengthGradientPinWA` returns `(gradientMinX - gradientMaxX)` —
//   i.e. -∂WA/∂x_pin — for exactly this reason
//   (nesterovBase.cpp:1523).
//
// Formulas:
//   - WA cost per net (ePlace Eq. 6, exponent-shifted, dreamplace/wa_functional.h:38-78):
//       WA_x = Σ x_p·exp((x_p-x_max)/γ) / Σ exp((x_p-x_max)/γ)
//            - Σ x_p·exp(-(x_p-x_min)/γ) / Σ exp(-(x_p-x_min)/γ)
//     Exponent floor: arg ≥ −300 (dreamplace/wa_functional.h:102-103 cutoff
//     pattern via minWireLengthForceBar in nesterovBase.cpp:1326).
//   - Per-pin gradient (closed form, dreamplace/wa_functional.h:177-185 and
//     nesterovBase.cpp:1460-1524):
//       d/dx_p (max_part) =
//         exp_x[p]·((1 + x_p/γ)·Σexp_x − (1/γ)·Σx·exp_x) / (Σexp_x)²
//       d/dx_p (min_part) =
//         exp_nx[p]·((1 − x_p/γ)·Σexp_nx + (1/γ)·Σx·exp_nx) / (Σexp_nx)²
//     ∂WA/∂x_p = d(max_part)/dx_p − d(min_part)/dx_p
//     descent g_p = −∂WA/∂x_p = d(min_part)/dx_p − d(max_part)/dx_p
//   - γ from overflow (ePlace Eq. 38 / nesterovPlace.cpp:1172-1186):
//     wireLengthCoef = 1/γ is scheduled; this function returns γ.
//
// Pin model: block center (x[i]+w/2, y[i]+h/2).  Since ∂pin/∂block = 1
// (constant offset), ∂WA/∂block = ∂WA/∂pin.
//
// Phase 3 supports 2-pin nets only (per plan §3) — multi-pin nets
// in Design::connections trigger an assert/abort.

#pragma once

#include "../types.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace mp {

struct WAResult {
    std::vector<double> grad_x;  // descent direction (see header note)
    std::vector<double> grad_y;
    double cost = 0.0;           // Σ w·(WA_x + WA_y) — HPWL surrogate
};

class WAWirelength {
public:
    explicit WAWirelength(const Design& d)
        : d_(d), nb_((int)d.blocks.size())
    {
        // Phase 3 is 2-pin only (per plan §3).  Connections are already 2-pin
        // by Design::Connection's {from,to} schema — but assert defensively in
        // case the data model grows.  Out-of-range or self-loop edges are
        // skipped during compute().
        for (const auto& c : d_.connections) {
            if (c.from == c.to) {
                std::fprintf(stderr,
                    "WAWirelength: self-loop net (from==to==%d) is undefined "
                    "for 2-pin WA; aborting (Phase 3 is 2-pin only).\n",
                    c.from);
                std::abort();
            }
        }
    }

    // Q1(b) plan-C accessors — route-feedback connection-weight boost.
    // set_boost(a, b, v) accumulates v into the symmetric pair (a, b).
    // After the route-feedback step, each subsequent compute() uses
    // w' = w · (1 + wl_boost_[a][b]).  v ≥ 0; the loop caps growth.
    void set_boost(int a, int b, double v) {
        ensure_boost();
        if (a < 0 || b < 0 || a >= nb_ || b >= nb_) return;
        wl_boost_[a][b] += v;
        wl_boost_[b][a] += v;
    }
    void clear_boost() {
        if ((int)wl_boost_.size() == nb_) {
            for (auto& row : wl_boost_) std::fill(row.begin(), row.end(), 0.0);
        }
    }

    WAResult compute(const std::vector<double>& x,
                     const std::vector<double>& y,
                     double gamma)
    {
        assert((int)x.size() == nb_ && (int)y.size() == nb_);
        assert(gamma > 0.0);
        const double inv_g = 1.0 / gamma;

        WAResult res;
        res.grad_x.assign(nb_, 0.0);
        res.grad_y.assign(nb_, 0.0);
        res.cost = 0.0;

        // Q1(b) plan-C: route-feedback per-connection weight boost set by the
        // outer RB loop.  Default boost is 0; the per-net weight becomes
        // w * (1 + wl_boost_[a][b]) when the matrix is allocated.  Plan-A's
        // global super-linear amplification was tested and regressed case4
        // (it over-concentrated multi-hub structures); the boost is now
        // targeted at the specific pairs whose routes overflow channels.
        for (const auto& c : d_.connections) {
            if (c.from < 0 || c.to < 0) continue;
            if (c.from >= nb_ || c.to >= nb_) continue;
            double w = (double)c.nets;
            if (w == 0.0) continue;
            if ((int)wl_boost_.size() == nb_) {
                w *= (1.0 + wl_boost_[c.from][c.to]);
            }

            const int a = c.from;
            const int b = c.to;

            // Pin centers (block center model).
            const double pax = x[a] + 0.5 * d_.blocks[a].width;
            const double pay = y[a] + 0.5 * d_.blocks[a].height;
            const double pbx = x[b] + 0.5 * d_.blocks[b].width;
            const double pby = y[b] + 0.5 * d_.blocks[b].height;

            // X axis.
            double gax_x = 0.0, gbx_x = 0.0; double cost_x = 0.0;
            wa_axis_2pin(pax, pbx, inv_g, cost_x, gax_x, gbx_x);

            // Y axis.
            double gay_y = 0.0, gby_y = 0.0; double cost_y = 0.0;
            wa_axis_2pin(pay, pby, inv_g, cost_y, gay_y, gby_y);

            res.cost += w * (cost_x + cost_y);
            // ∂pin/∂block = 1 (block center model), so propagate directly.
            res.grad_x[a] += w * gax_x;
            res.grad_x[b] += w * gbx_x;
            res.grad_y[a] += w * gay_y;
            res.grad_y[b] += w * gby_y;
        }
        return res;
    }

    // True bounding-box HPWL (weighted) — exposed for the γ→0 acceptance test.
    double hpwl_exact(const std::vector<double>& x,
                      const std::vector<double>& y) const
    {
        assert((int)x.size() == nb_ && (int)y.size() == nb_);
        double total = 0.0;
        for (const auto& c : d_.connections) {
            if (c.from < 0 || c.to < 0) continue;
            if (c.from >= nb_ || c.to >= nb_) continue;
            const int a = c.from;
            const int b = c.to;
            const double pax = x[a] + 0.5 * d_.blocks[a].width;
            const double pay = y[a] + 0.5 * d_.blocks[a].height;
            const double pbx = x[b] + 0.5 * d_.blocks[b].width;
            const double pby = y[b] + 0.5 * d_.blocks[b].height;
            total += (double)c.nets *
                     (std::abs(pax - pbx) + std::abs(pay - pby));
        }
        return total;
    }

private:
    // Exponent floor — clamp argument to ≥ −300 before std::exp (guard #2 /
    // dreamplace cutoff pattern).
    static double safe_exp(double arg) {
        if (arg < -300.0) arg = -300.0;
        return std::exp(arg);
    }

    // 2-pin WA along one axis.  Returns cost and the *descent* gradient
    // for each pin (g = -∂WA/∂pin); ∂pin/∂block = 1 so caller adds these
    // directly to per-block gradients.
    //
    // Closed form derived from dreamplace/wa_functional.h:177-185 and
    // nesterovBase.cpp:1460-1524 (the same expression, returned in the
    // "minX - maxX" descent form).
    static void wa_axis_2pin(double p0, double p1, double inv_g,
                             double& cost, double& g0, double& g1)
    {
        // Bounding box: pmax = max(p0,p1), pmin = min(p0,p1).
        const double pmax = std::max(p0, p1);
        const double pmin = std::min(p0, p1);

        // exp((p - pmax)/γ) for the max side: pin at pmax contributes 1.
        const double e0max = safe_exp((p0 - pmax) * inv_g);
        const double e1max = safe_exp((p1 - pmax) * inv_g);
        // exp(-(p - pmin)/γ) for the min side: pin at pmin contributes 1.
        const double e0min = safe_exp(-(p0 - pmin) * inv_g);
        const double e1min = safe_exp(-(p1 - pmin) * inv_g);

        const double S_max  = e0max + e1max;            // Σ exp_x
        const double S_min  = e0min + e1min;            // Σ exp_nx
        const double Sx_max = p0 * e0max + p1 * e1max;  // Σ x·exp_x
        const double Sx_min = p0 * e0min + p1 * e1min;  // Σ x·exp_nx

        const double max_part = Sx_max / S_max;
        const double min_part = Sx_min / S_min;
        cost = max_part - min_part;

        // ∂(max_part)/∂p_k = exp_x[k]·((1 + p_k/γ)·S_max - (1/γ)·Sx_max)
        //                    / S_max² .
        // ∂(min_part)/∂p_k = exp_nx[k]·((1 - p_k/γ)·S_min + (1/γ)·Sx_min)
        //                    / S_min² .
        // descent g_k = ∂(min_part) - ∂(max_part).
        auto dmax_dp = [&](double pk, double ek) {
            return ek * ((1.0 + inv_g * pk) * S_max - inv_g * Sx_max)
                   / (S_max * S_max);
        };
        auto dmin_dp = [&](double pk, double ek) {
            return ek * ((1.0 - inv_g * pk) * S_min + inv_g * Sx_min)
                   / (S_min * S_min);
        };
        const double d0_max = dmax_dp(p0, e0max);
        const double d1_max = dmax_dp(p1, e1max);
        const double d0_min = dmin_dp(p0, e0min);
        const double d1_min = dmin_dp(p1, e1min);
        g0 = d0_min - d0_max;
        g1 = d1_min - d1_max;
    }

    void ensure_boost() {
        if ((int)wl_boost_.size() != nb_) {
            wl_boost_.assign(nb_, std::vector<double>(nb_, 0.0));
        }
    }

    const Design& d_;
    int nb_;
    std::vector<std::vector<double>> wl_boost_;  // n×n, default empty
};

// γ schedule per ePlace Eq. 38 / nesterovPlace.cpp:1172-1186.  `base` is
// 0.25 / avg_bin_size (Replace.h default).  OpenROAD schedules
// wireLengthCoef = 1/γ; this helper returns γ itself.
//
//   τ > 1.0  →  1/γ = 0.1·base                          → γ = 10/base
//   τ < 0.1  →  1/γ = 10·base                           → γ = 0.1/base
//   else     →  1/γ = base / pow(10, (τ-0.1)·20/9 − 1)
//                                                       → γ = pow(10, (τ-0.1)·20/9 − 1) / base
//
inline double gamma_from_overflow(double tau, double avg_bin_size) {
    assert(avg_bin_size > 0.0);
    const double base = 0.25 / avg_bin_size;
    double inv_gamma;
    if (tau > 1.0) {
        inv_gamma = 0.1 * base;
    } else if (tau < 0.1) {
        inv_gamma = 10.0 * base;
    } else {
        inv_gamma = base / std::pow(10.0, (tau - 0.1) * 20.0 / 9.0 - 1.0);
    }
    return 1.0 / inv_gamma;
}

} // namespace mp
