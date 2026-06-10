// SPDX-License-Identifier: BSD-3-Clause
// Phase 4: ePlace-style Nesterov global placement loop.
//
// SIGN CONVENTION (single source of truth — guard #1, matches density.h and
// wa_wirelength.h):
//   Density::grad_x/grad_y and WAWirelength::grad_x/grad_y are *descent*
//   directions.  This loop also stores descent directions in `cur_grad_x/y`
//   and updates positions as `x_new = x + step * g`.  Total gradient per axis
//       g_i = g_w_i + lambda * q_i * g_d_i
//   where q_i = inflated block area (charge), matching ePlace Alg. 2 /
//   nesterovBase.cpp:2785-2786 (descent form).  No sign flip anywhere.
//
// Anti-pattern guards followed:
//   #1 sign convention — descent, documented here.
//   #4 do not under-set gamma — call gamma_from_overflow() every iter.
//   #7 divergence guards — best-tau snapshot, NaN/Inf -> shrink lambda0
//      and retry from x0 with x10 perturbation, max 1 retry.
//
// Citations: ePlace Alg. 2 + Eq. 29-36; OpenROAD
// nesterovPlace.cpp:871-928 (backtracking), 1064-1138 (main loop),
// 1140-1170 (initial lambda search / wirelength bootstrap),
// 560-640 (divergence guards); nesterovBase.cpp:1526-1530, 2572-2578,
// 2788-2801 (Jacobi preconditioner), 2701-2723 (lambda0 init Eq. 35),
// 2725-2743 (steplength), 2975-2986, 3098 (lambda * phiCoef).

#pragma once

#include "density.h"
#include "wa_wirelength.h"
#include "../types.h"
#include "../config.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

namespace mp {

struct NesterovParams {
    int    max_iter        = 1500;    // cfg::MP_MAX_ITER
    double target_overflow = 0.08;    // cfg::MP_TARGET_OVF
    double init_lambda     = 8e-5;    // cfg::MP_INIT_LAMBDA  (ePlace Eq. 35)
    double phi_min         = 0.95;
    double phi_max         = 1.05;
    int    backtrack_max   = 10;      // nesterovPlace.cpp:875 maxBackTrack
    double backtrack_floor = 0.01;    // floor on accepted step length
    double backtrack_accept = 0.95;   // newStep > 0.95 * step -> accept
};

struct NesterovResult {
    std::vector<double> x;   // final lower-left x per block
    std::vector<double> y;
    double tau_final  = 1.0;
    double hpwl_final = 0.0;
    int    iters      = 0;
    bool   converged  = false;   // tau <= target_overflow
    bool   diverged   = false;   // hit NaN/Inf path and could not recover
};

class Nesterov {
public:
    Nesterov(Density& dens, WAWirelength& wl, const Design& d, NesterovParams p)
        : dens_(dens), wl_(wl), d_(d), p_(p),
          nb_(static_cast<int>(d.blocks.size())),
          outline_w_(d.outline.max_width),
          outline_h_(d.outline.max_height)
    {
        assert(nb_ > 0);
        assert(outline_w_ > 0.0 && outline_h_ > 0.0);

        // Cache inflated block "charge" q_i = inflated block area.
        // ePlace Eq. 35 / nesterovBase.cpp:2701-2723 — q_i denotes the
        // movable charge used to scale the density gradient when forming
        // lambda0 = sum|wlGrad| / sum(q_i |denGrad_i|) * init_lambda.
        const auto& halo = dens_.halo();
        q_.assign(nb_, 0.0);
        is_fixed_.assign(nb_, false);
        for (int i = 0; i < nb_; ++i) {
            const Block& b = d_.blocks[i];
            double iw = b.width  + 2.0 * halo[i];
            double ih = b.height + 2.0 * halo[i];
            q_[i] = iw * ih;
            if (b.type == BlockType::EDGE) is_fixed_[i] = true;
        }

        // Pin counts per block (for Jacobi preconditioner — ePlace Eq. 30-33).
        pin_count_.assign(nb_, 0.0);
        for (const auto& c : d_.connections) {
            if (c.from >= 0 && c.from < nb_) pin_count_[c.from] += 1.0;
            if (c.to   >= 0 && c.to   < nb_) pin_count_[c.to]   += 1.0;
        }

        avg_bin_size_ = 0.5 * (dens_.bin_w() + dens_.bin_h());
        assert(avg_bin_size_ > 0.0);
    }

    NesterovResult run(const std::vector<double>& x0,
                       const std::vector<double>& y0)
    {
        assert(static_cast<int>(x0.size()) == nb_);
        assert(static_cast<int>(y0.size()) == nb_);

        const bool debug = (std::getenv("MP_DEBUG") != nullptr &&
                            std::getenv("MP_DEBUG")[0] != '\0' &&
                            std::getenv("MP_DEBUG")[0] != '0');

        // Up to 2 attempts: original lambda0 scale, then 0.1x with x10
        // perturbation if first attempt diverged.  Mirrors the
        // "revert + retry" pattern in nesterovPlace.cpp:560-640 collapsed
        // for the no-routability path in Phase 4.
        NesterovResult res;
        double lambda_scale  = p_.init_lambda;
        double perturb_scale = 1.0;
        for (int attempt = 0; attempt < 2; ++attempt) {
            res = run_once(x0, y0, lambda_scale, perturb_scale, debug);
            if (!res.diverged) {
                return res;
            }
            if (debug) {
                std::fprintf(stderr,
                    "[Nesterov] attempt %d diverged; retrying with "
                    "lambda*=0.1, perturb*=10\n", attempt);
            }
            lambda_scale  *= 0.1;
            perturb_scale *= 10.0;
        }
        return res;
    }

private:
    // ---- One full Nesterov attempt --------------------------------------------
    // OpenROAD-style bookkeeping (nesterovBase.cpp:3137-3171):
    //   cur, next     — bare descent positions
    //   cur_slp, next_slp — lookahead (SLP) positions where gradient is
    //                       evaluated; gradient at SLP drives the step.
    //   prev_slp      — previous lookahead (for steplength estimate).
    //   Trial:  next     = cur_slp + step * grad(cur_slp)
    //           next_slp = next + coeff * (next - cur)
    //   Then:   cur <- next, cur_slp <- next_slp, prev_slp <- cur_slp.
    NesterovResult run_once(const std::vector<double>& x0,
                            const std::vector<double>& y0,
                            double lambda_scale,
                            double perturb_scale,
                            bool debug)
    {
        // ==== Init bare and SLP positions (nesterovBase.cpp:2669) ====
        //   curSLPCoordi = prevSLPCoordi = curCoordi = initial pos.
        std::vector<double> x_cur     = x0, y_cur     = y0;
        std::vector<double> x_cur_slp = x0, y_cur_slp = y0;
        std::vector<double> x_prev_slp(nb_), y_prev_slp(nb_);

        // Synthetic prev-SLP: small deterministic perturbation so the initial
        // ||dx|| in compute_step is nonzero.  source: nesterovPlace.cpp:1140-1170
        // (prevSLP bootstrap via initialPrevCoordiUpdateCoef path).
        std::mt19937 rng(static_cast<unsigned>(nb_) * 7919u + 1u);
        std::uniform_real_distribution<double> pert(-1.0, 1.0);
        const double base_pert = 0.001 * avg_bin_size_ * perturb_scale;
        for (int i = 0; i < nb_; ++i) {
            x_prev_slp[i] = x_cur_slp[i] - base_pert * pert(rng);
            y_prev_slp[i] = y_cur_slp[i] - base_pert * pert(rng);
            if (is_fixed_[i]) {
                x_prev_slp[i] = x_cur_slp[i];
                y_prev_slp[i] = y_cur_slp[i];
            }
        }

        // ==== Gradient at cur_slp; lambda0 bootstrap (Eq. 35) ====
        std::vector<double> gx_cur_slp(nb_), gy_cur_slp(nb_);
        std::vector<double> gx_prev_slp(nb_), gy_prev_slp(nb_);

        auto den_cs = dens_.compute(x_cur_slp, y_cur_slp);
        double tau  = den_cs.tau;
        double gamma = gamma_from_overflow(tau, avg_bin_size_);
        auto wl_cs  = wl_.compute(x_cur_slp, y_cur_slp, gamma);

        // ePlace Eq. 35 / nesterovBase.cpp:2701-2723:
        //   lambda0 = wireLengthGradSum / densityGradSum * init_lambda.
        // The density gradient already integrates `overlapArea * field` (q_i
        // is baked in), matching OpenROAD's densityGradSum_ form at
        // nesterovBase.cpp:2810-2814.
        double sum_wl  = 0.0, sum_den = 0.0;
        for (int i = 0; i < nb_; ++i) {
            if (is_fixed_[i]) continue;
            sum_wl  += std::fabs(wl_cs.grad_x[i]) + std::fabs(wl_cs.grad_y[i]);
            sum_den += std::fabs(den_cs.grad_x[i]) + std::fabs(den_cs.grad_y[i]);
        }
        double lambda = (sum_den > 0.0)
                        ? (sum_wl / sum_den) * lambda_scale
                        : lambda_scale;
        if (!std::isfinite(lambda) || lambda <= 0.0) lambda = lambda_scale;

        combine_gradient(wl_cs.grad_x, wl_cs.grad_y,
                         den_cs.grad_x, den_cs.grad_y,
                         lambda, gx_cur_slp, gy_cur_slp);
        precondition(lambda, gx_cur_slp, gy_cur_slp);
        zero_fixed_grad(gx_cur_slp, gy_cur_slp);

        auto den_ps = dens_.compute(x_prev_slp, y_prev_slp);
        auto wl_ps  = wl_.compute(x_prev_slp, y_prev_slp, gamma);
        combine_gradient(wl_ps.grad_x, wl_ps.grad_y,
                         den_ps.grad_x, den_ps.grad_y,
                         lambda, gx_prev_slp, gy_prev_slp);
        precondition(lambda, gx_prev_slp, gy_prev_slp);
        zero_fixed_grad(gx_prev_slp, gy_prev_slp);

        // ==== Initial steplength = ||delta_coord|| / ||delta_grad|| ====
        // source: nesterovBase.cpp:2725-2743.
        double step = compute_step(x_prev_slp, y_prev_slp, x_cur_slp, y_cur_slp,
                                   gx_prev_slp, gy_prev_slp,
                                   gx_cur_slp,  gy_cur_slp);
        if (!std::isfinite(step) || step <= 0.0) step = 0.01 * avg_bin_size_;

        // ==== Best snapshot (guard #7) ====
        NesterovResult best;
        best.x = x_cur;
        best.y = y_cur;
        best.tau_final  = tau;
        best.hpwl_final = wl_.hpwl_exact(x_cur, y_cur);
        best.iters      = 0;
        best.converged  = false;
        best.diverged   = false;
        double best_tau = tau;

        // Nesterov sequence init.  ePlace Alg. 2: a_0 = 1.
        double a_k = 1.0;

        double prev_hpwl = best.hpwl_final;
        int    stall_cnt = 0;
        // OpenROAD uses a constant referenceHpwl (≈ design HPWL @ initial
        // placement).  Setting ref_hpwl = HPWL at iter 0 keeps the scaled
        // diff comparable across iters; OpenROAD's default 446e6 is design-
        // scale.  source: Replace.h:77, nesterovBase.cpp:3082-3083.
        double ref_hpwl = std::max(1.0, best.hpwl_final);

        if (debug) {
            std::fprintf(stderr,
                "[Nesterov] start  nb=%d lambda0=%.3e step0=%.3e gamma0=%.3e "
                "tau0=%.3f HPWL0=%.3e\n",
                nb_, lambda, step, gamma, tau, best.hpwl_final);
        }

        // Buffers for the inner loop.
        std::vector<double> x_next(nb_), y_next(nb_);
        std::vector<double> x_next_slp(nb_), y_next_slp(nb_);
        std::vector<double> gx_next_slp(nb_), gy_next_slp(nb_);

        int iter = 0;
        for (; iter < p_.max_iter; ++iter) {
            // ePlace Alg. 2: a_{k+1} = (1 + sqrt(4*a_k^2 + 1))/2.
            // source: nesterovPlace.cpp:1071.
            double a_kp1 = 0.5 * (1.0 + std::sqrt(4.0 * a_k * a_k + 1.0));
            double coeff = (a_k - 1.0) / a_kp1;   // nesterovPlace.cpp:1074

            // ==== Backtracking (nesterovPlace.cpp:871-928) ====
            bool   diverged_inner = false;
            double try_step  = step;
            double new_step  = step;
            double tau_new   = tau;
            double gamma_new = gamma;
            int    bt = 0;
            for (bt = 0; bt < p_.backtrack_max; ++bt) {
                // Trial move (nesterovBase.cpp:3156-3170):
                //   next     = cur_slp + step * grad(cur_slp)
                //   next_slp = next + coeff * (next - cur)
                for (int i = 0; i < nb_; ++i) {
                    if (is_fixed_[i]) {
                        x_next[i]     = x_cur[i];
                        y_next[i]     = y_cur[i];
                        x_next_slp[i] = x_cur_slp[i];
                        y_next_slp[i] = y_cur_slp[i];
                        continue;
                    }
                    double nx = x_cur_slp[i] + try_step * gx_cur_slp[i];
                    double ny = y_cur_slp[i] + try_step * gy_cur_slp[i];
                    nx = clamp_x(nx, d_.blocks[i].width);
                    ny = clamp_y(ny, d_.blocks[i].height);
                    x_next[i] = nx;
                    y_next[i] = ny;
                    double sx = nx + coeff * (nx - x_cur[i]);
                    double sy = ny + coeff * (ny - y_cur[i]);
                    x_next_slp[i] = clamp_x(sx, d_.blocks[i].width);
                    y_next_slp[i] = clamp_y(sy, d_.blocks[i].height);
                }

                // Gradient at next_slp (this is where the next-iter step is
                // measured).  source: nesterovPlace.cpp:879-882.
                auto den_ns = dens_.compute(x_next_slp, y_next_slp);
                if (!std::isfinite(den_ns.tau)) { diverged_inner = true; break; }
                gamma_new = gamma_from_overflow(den_ns.tau, avg_bin_size_);
                auto wl_ns = wl_.compute(x_next_slp, y_next_slp, gamma_new);
                combine_gradient(wl_ns.grad_x, wl_ns.grad_y,
                                 den_ns.grad_x, den_ns.grad_y,
                                 lambda, gx_next_slp, gy_next_slp);
                precondition(lambda, gx_next_slp, gy_next_slp);
                zero_fixed_grad(gx_next_slp, gy_next_slp);
                if (has_non_finite(gx_next_slp) || has_non_finite(gy_next_slp)) {
                    diverged_inner = true;
                    break;
                }

                // Measured steplength between cur_slp and next_slp.
                new_step = compute_step(x_cur_slp, y_cur_slp,
                                        x_next_slp, y_next_slp,
                                        gx_cur_slp, gy_cur_slp,
                                        gx_next_slp, gy_next_slp);
                if (!std::isfinite(new_step) || new_step <= 0.0) {
                    diverged_inner = true;
                    break;
                }
                tau_new = den_ns.tau;

                if (new_step > p_.backtrack_accept * try_step) {
                    break;   // accept
                }
                double next_try = 0.5 * try_step;
                if (next_try < p_.backtrack_floor * avg_bin_size_) {
                    break;   // floored — accept current measurement
                }
                try_step = next_try;
            }

            if (diverged_inner) {
                if (debug) {
                    std::fprintf(stderr,
                        "[Nesterov] iter %d: NaN/Inf in backtracking, "
                        "reverting to best snapshot (tau=%.3f).\n",
                        iter, best_tau);
                }
                best.diverged = true;
                best.iters    = iter;
                return best;
            }

            tau   = tau_new;
            gamma = gamma_new;
            step  = new_step;

            // ==== HPWL on the *bare* next position (evaluation point) ====
            double hpwl_new = wl_.hpwl_exact(x_next, y_next);

            // ==== lambda update (phiCoef) ====
            // source: nesterovBase.cpp:2975-2986, 3082-3083.
            double scaled_diff_hpwl = (hpwl_new - prev_hpwl) / ref_hpwl;
            double phi_coef = (scaled_diff_hpwl < 0.0)
                ? p_.phi_max
                : p_.phi_max * std::pow(p_.phi_max, -scaled_diff_hpwl);
            phi_coef = std::max(p_.phi_min, std::min(p_.phi_max, phi_coef));
            lambda *= phi_coef;
            if (!std::isfinite(lambda) || lambda <= 0.0) {
                if (debug) {
                    std::fprintf(stderr,
                        "[Nesterov] iter %d: lambda became invalid (%.3e), "
                        "reverting.\n", iter, lambda);
                }
                best.diverged = true;
                best.iters    = iter;
                return best;
            }

            // ==== Best snapshot ====
            if (tau < best_tau) {
                best_tau        = tau;
                best.x          = x_next;
                best.y          = y_next;
                best.tau_final  = tau;
                best.hpwl_final = hpwl_new;
                best.iters      = iter + 1;
            }

            // ==== Roll forward (swap pointers) ====
            // source: nesterovBase.cpp:2995-2999 (SLP swap).
            std::swap(x_prev_slp, x_cur_slp);
            std::swap(y_prev_slp, y_cur_slp);
            std::swap(gx_prev_slp, gx_cur_slp);
            std::swap(gy_prev_slp, gy_cur_slp);
            std::swap(x_cur,     x_next);
            std::swap(y_cur,     y_next);
            std::swap(x_cur_slp, x_next_slp);
            std::swap(y_cur_slp, y_next_slp);
            std::swap(gx_cur_slp, gx_next_slp);
            std::swap(gy_cur_slp, gy_next_slp);
            a_k = a_kp1;

            // ==== Stall detection ====
            double rel_dhpwl = (prev_hpwl > 1e-12)
                ? std::fabs(hpwl_new - prev_hpwl) / prev_hpwl
                : 0.0;
            if (rel_dhpwl < 1e-6) ++stall_cnt; else stall_cnt = 0;
            prev_hpwl = hpwl_new;

            // ==== Debug log (sparse) ====
            if (debug && (iter % 50 == 0 || iter < 5)) {
                std::fprintf(stderr,
                    "[Nesterov] iter %4d  tau=%.4f  HPWL=%.4e  "
                    "lambda=%.3e  step=%.3e  gamma=%.3e  phi=%.4f  bt=%d\n",
                    iter, tau, hpwl_new, lambda, step, gamma, phi_coef, bt);
            }

            // ==== Termination ====
            if (tau <= p_.target_overflow) {
                best.x          = x_cur;
                best.y          = y_cur;
                best.tau_final  = tau;
                best.hpwl_final = hpwl_new;
                best.iters      = iter + 1;
                best.converged  = true;
                if (debug) {
                    std::fprintf(stderr,
                        "[Nesterov] converged at iter %d (tau=%.4f <= %.4f)\n",
                        iter, tau, p_.target_overflow);
                }
                return best;
            }
            if (stall_cnt >= 50) {
                if (debug) {
                    std::fprintf(stderr,
                        "[Nesterov] stall at iter %d (tau=%.4f)\n", iter, tau);
                }
                best.iters = iter + 1;
                return best;
            }
        }

        best.iters = iter;
        return best;
    }

    // ---- helpers --------------------------------------------------------------

    void combine_gradient(const std::vector<double>& wlx,
                          const std::vector<double>& wly,
                          const std::vector<double>& dx,
                          const std::vector<double>& dy,
                          double lambda,
                          std::vector<double>& gx,
                          std::vector<double>& gy) const
    {
        // source: nesterovBase.cpp:2785-2786 (descent form):
        //   sumGrads = wireLengthGrads + densityPenalty * densityGrads
        // The q_i (block area) factor lives in the *preconditioner*
        // (nesterovBase.cpp:2572-2578 getDensityPreconditioner), NOT here.
        for (int i = 0; i < nb_; ++i) {
            gx[i] = wlx[i] + lambda * dx[i];
            gy[i] = wly[i] + lambda * dy[i];
        }
    }

    void precondition(double lambda,
                      std::vector<double>& gx,
                      std::vector<double>& gy) const
    {
        // Jacobi preconditioner: h_i = max(1, #pins_i + lambda * q_i).
        // source: ePlace Eq. 30-33; nesterovBase.cpp:1526-1530, 2572-2578,
        // 2788-2801 (sumPrecondi = wlPre + lam*denPre; floor at minPrecondi).
        for (int i = 0; i < nb_; ++i) {
            double h = pin_count_[i] + lambda * q_[i];
            if (h < 1.0) h = 1.0;
            gx[i] /= h;
            gy[i] /= h;
        }
    }

    void zero_fixed_grad(std::vector<double>& gx,
                         std::vector<double>& gy) const
    {
        // EDGE blocks are pinned at their initial location for Phase 4
        // (full snap-to-nearest comes in Phase 6 per plan §6).
        for (int i = 0; i < nb_; ++i) {
            if (is_fixed_[i]) {
                gx[i] = 0.0;
                gy[i] = 0.0;
            }
        }
    }

    static bool has_non_finite(const std::vector<double>& v) {
        for (double x : v) if (!std::isfinite(x)) return true;
        return false;
    }

    static double compute_step(const std::vector<double>& xp,
                               const std::vector<double>& yp,
                               const std::vector<double>& xc,
                               const std::vector<double>& yc,
                               const std::vector<double>& gp_x,
                               const std::vector<double>& gp_y,
                               const std::vector<double>& gc_x,
                               const std::vector<double>& gc_y)
    {
        // source: nesterovBase.cpp:2725-2743 — ||dx|| / ||dg||.
        double dx2 = 0.0, dg2 = 0.0;
        const size_t n = xp.size();
        for (size_t i = 0; i < n; ++i) {
            double a = xc[i] - xp[i];
            double b = yc[i] - yp[i];
            dx2 += a * a + b * b;
            double ga = gc_x[i] - gp_x[i];
            double gb = gc_y[i] - gp_y[i];
            dg2 += ga * ga + gb * gb;
        }
        if (dg2 <= 0.0) return 0.0;
        return std::sqrt(dx2) / std::sqrt(dg2);
    }

    double clamp_x(double v, double w) const {
        double hi = outline_w_ - w;
        if (hi < 0.0) hi = 0.0;
        if (v < 0.0) return 0.0;
        if (v > hi)  return hi;
        return v;
    }
    double clamp_y(double v, double h) const {
        double hi = outline_h_ - h;
        if (hi < 0.0) hi = 0.0;
        if (v < 0.0) return 0.0;
        if (v > hi)  return hi;
        return v;
    }

    Density&        dens_;
    WAWirelength&   wl_;
    const Design&   d_;
    NesterovParams  p_;

    int    nb_;
    double outline_w_, outline_h_;
    double avg_bin_size_;

    std::vector<double> q_;            // inflated block area (charge)
    std::vector<double> pin_count_;    // per-block pin count
    std::vector<bool>   is_fixed_;     // EDGE blocks: gradient zeroed
};

} // namespace mp
