// SPDX-License-Identifier: BSD-3-Clause
// Phase 2: electrostatic density layer (ePlace) for the MP engine.
//
// SIGN CONVENTION (single source of truth — guard #1):
//   We follow OpenROAD / ePlace: the spectral solve yields a potential ψ and a
//   field (ξx, ξy) = (ψ·wx, ψ·wy) per bin via DST/DCT pairings (see
//   docs/refs/openroad-gpl/fft.cpp:169–182).  Per-block we accumulate
//       grad_x[i] = Σ_b overlapArea(b,i) · ξx[b]
//       grad_y[i] = Σ_b overlapArea(b,i) · ξy[b]
//   (loop pattern from nesterovBase.cpp:2582–2601).  This is the
//   *descent direction*: callers update positions as `x[i] += α·g[i]` (no
//   sign flip).  Bin charge ρ' counts movable area only; fixed objects'
//   charge is scaled by ρt (ePlace §3.2 / guard #6).
//
// Formulas cited inline.  Vendored Poisson/FFT pipeline is in src/vendor/fft.* —
// we never re-derive its normalization (guard #3).

#pragma once

#include "../types.h"
#include "../config.h"
#include "../vendor/fft.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <vector>

namespace mp {

struct DensityResult {
    std::vector<double> grad_x; // descent direction per block (see header note)
    std::vector<double> grad_y;
    double tau = 0.0;           // ePlace Eq. 37: movable overflow ratio
    double rho_t = 0.0;         // target density used this call
};

class Density {
public:
    // grid: square count of bins per axis (MP_GRID); bin sizes follow the
    // outline aspect ratio (frequency rescaling handled inside fft.cpp:54–57).
    Density(const Design& d, int grid)
        : d_(d), nb_((int)d.blocks.size()), grid_(grid)
    {
        assert(grid_ > 0);
        if (cfg::MP_HALO_SCALE < 0.0) {
            std::fprintf(stderr,
                "Density: MP_HALO_SCALE must be >= 0 (got %g)\n",
                cfg::MP_HALO_SCALE);
            std::abort();
        }
        outline_w_ = d.outline.max_width;
        outline_h_ = d.outline.max_height;
        if (!(outline_w_ > 0.0) || !(outline_h_ > 0.0)) {
            std::fprintf(stderr,
                "Density: outline must have positive width/height (%g x %g)\n",
                outline_w_, outline_h_);
            std::abort();
        }
        bin_w_ = outline_w_ / grid_;
        bin_h_ = outline_h_ / grid_;
        bin_area_ = bin_w_ * bin_h_;

        // Per-block total incident nets (used by halo computation).
        nets_total_.assign(nb_, 0.0);
        for (const auto& c : d.connections) {
            if (c.from >= 0 && c.from < nb_) nets_total_[c.from] += c.nets;
            if (c.to   >= 0 && c.to   < nb_) nets_total_[c.to]   += c.nets;
        }

        // Halo[i] = per-side gap, ePlace-like inflation.  Demand bound from
        // architecture §: h_i = max(HALO/2, MP_HALO_SCALE · T_i / (2·k_i·25))
        halo_.assign(nb_, 0.0);
        for (int i = 0; i < nb_; ++i) {
            const Block& b = d.blocks[i];
            double k = (b.type == BlockType::EDGE) ? 3.0 : 4.0;
            double demand = cfg::MP_HALO_SCALE * nets_total_[i]
                            / (2.0 * k * 25.0);
            halo_[i] = std::max(cfg::HALO * 0.5, demand);
        }

        // ePlace §IV / inflated utilization: ρt default = Σ inflated movable
        // area / outline area, clamped to ≤ 0.9.
        movable_area_inflated_ = 0.0;
        for (int i = 0; i < nb_; ++i) {
            if (d.blocks[i].type == BlockType::EDGE) continue;
            double iw = d.blocks[i].width  + 2.0 * halo_[i];
            double ih = d.blocks[i].height + 2.0 * halo_[i];
            movable_area_inflated_ += iw * ih;
        }
        rho_t_ = cfg::MP_TDENSITY > 0.0
            ? cfg::MP_TDENSITY
            : std::min(0.9,
                       movable_area_inflated_ / (outline_w_ * outline_h_));

        // Vendored DCT/Poisson pipeline.
        fft_.reset(new gpl::FFT(grid_, grid_,
                                static_cast<float>(bin_w_),
                                static_cast<float>(bin_h_)));

        rho_bins_.assign(grid_ * grid_, 0.0);
        xi_x_.assign(grid_ * grid_, 0.0);
        xi_y_.assign(grid_ * grid_, 0.0);
        psi_.assign(grid_ * grid_, 0.0);
    }

    // Recompute ρt (e.g. after halo bumps).  rho_t == 0 → auto path.
    void update_target_density(double rho_t = 0.0) {
        movable_area_inflated_ = 0.0;
        for (int i = 0; i < nb_; ++i) {
            if (d_.blocks[i].type == BlockType::EDGE) continue;
            double iw = d_.blocks[i].width  + 2.0 * halo_[i];
            double ih = d_.blocks[i].height + 2.0 * halo_[i];
            movable_area_inflated_ += iw * ih;
        }
        rho_t_ = (rho_t > 0.0)
            ? rho_t
            : std::min(0.9,
                       movable_area_inflated_ / (outline_w_ * outline_h_));
    }

    DensityResult compute(const std::vector<double>& x,
                          const std::vector<double>& y)
    {
        assert((int)x.size() == nb_ && (int)y.size() == nb_);

        std::fill(rho_bins_.begin(), rho_bins_.end(), 0.0);

        double movable_area_sum = 0.0;

        // Deposit each block's charge into bins it overlaps.
        for (int i = 0; i < nb_; ++i) {
            const Block& b = d_.blocks[i];
            double iw = b.width  + 2.0 * halo_[i];
            double ih = b.height + 2.0 * halo_[i];
            double lx = x[i] - halo_[i];
            double ly = y[i] - halo_[i];

            // ePlace Eq. 25 (local smoothing): if the cell is sub-bin, inflate
            // its density footprint to √2·bin and scale the per-bin charge by
            // (cell_dim / (√2·bin_dim)) so total deposited area stays = cell.
            // (nesterovBase.cpp:2471–2501)
            const double kSqrt2 = 1.41421356237309504880;
            double dens_w = iw, dens_h = ih;
            double scale_x = 1.0, scale_y = 1.0;
            if (iw < kSqrt2 * bin_w_) {
                scale_x = iw / (kSqrt2 * bin_w_);
                dens_w = kSqrt2 * bin_w_;
            }
            if (ih < kSqrt2 * bin_h_) {
                scale_y = ih / (kSqrt2 * bin_h_);
                dens_h = kSqrt2 * bin_h_;
            }
            // Re-center the density box so the deposited centroid stays put.
            double dlx = lx + 0.5 * iw - 0.5 * dens_w;
            double dly = ly + 0.5 * ih - 0.5 * dens_h;
            double dux = dlx + dens_w;
            double duy = dly + dens_h;

            // Per-block density scale.  Fixed objects (EDGE blocks) are
            // weighted by ρt (ePlace §3.2 / guard #6) so they do not
            // over-repel and create dead halos.
            double charge_scale = scale_x * scale_y;
            if (b.type == BlockType::EDGE) {
                charge_scale *= rho_t_;
            } else {
                movable_area_sum += b.width * b.height; // raw area, not inflated
            }

            deposit(dlx, dly, dux, duy, charge_scale);
        }

        // Movable-only ρ used by both the overflow τ and the Poisson source.
        // ePlace Eq. 37 — compare per-bin movable density against ρt.
        double overflow_num = 0.0;
        for (int j = 0; j < grid_; ++j) {
            for (int i = 0; i < grid_; ++i) {
                double rho_mov = rho_bins_[bin_idx(i, j)] / bin_area_;
                if (rho_mov > rho_t_) {
                    overflow_num += (rho_mov - rho_t_) * bin_area_;
                }
            }
        }
        double tau = (movable_area_sum > 0.0)
                     ? overflow_num / movable_area_sum
                     : 0.0;

        // Hand bin density to FFT (already per-bin sum of area contributions —
        // doFFT() expects raw density values, and its normalization at lines
        // 109–120 handles the rest).  We pass area-weighted density (charge
        // per bin = deposited area) so that the resulting field ξ has units
        // consistent with ∂Φ/∂x in the bin-area integrand.
        for (int j = 0; j < grid_; ++j) {
            for (int i = 0; i < grid_; ++i) {
                fft_->updateDensity(i, j,
                    static_cast<float>(rho_bins_[bin_idx(i, j)]));
            }
        }
        fft_->doFFT();

        // Cache ψ and ξ for tests + per-block accumulation.
        for (int j = 0; j < grid_; ++j) {
            for (int i = 0; i < grid_; ++i) {
                psi_[bin_idx(i, j)] = fft_->getElectroPhi(i, j);
                auto e = fft_->getElectroField(i, j);
                xi_x_[bin_idx(i, j)] = e.first;
                xi_y_[bin_idx(i, j)] = e.second;
            }
        }

        // Per-block density gradient — nesterovBase.cpp:2582–2601 pattern.
        DensityResult res;
        res.grad_x.assign(nb_, 0.0);
        res.grad_y.assign(nb_, 0.0);
        res.tau = tau;
        res.rho_t = rho_t_;

        for (int idx = 0; idx < nb_; ++idx) {
            const Block& b = d_.blocks[idx];
            double iw = b.width  + 2.0 * halo_[idx];
            double ih = b.height + 2.0 * halo_[idx];
            double lx = x[idx] - halo_[idx];
            double ly = y[idx] - halo_[idx];

            const double kSqrt2 = 1.41421356237309504880;
            double dens_w = iw, dens_h = ih;
            double scale_x = 1.0, scale_y = 1.0;
            if (iw < kSqrt2 * bin_w_) {
                scale_x = iw / (kSqrt2 * bin_w_);
                dens_w = kSqrt2 * bin_w_;
            }
            if (ih < kSqrt2 * bin_h_) {
                scale_y = ih / (kSqrt2 * bin_h_);
                dens_h = kSqrt2 * bin_h_;
            }
            double dlx = lx + 0.5 * iw - 0.5 * dens_w;
            double dly = ly + 0.5 * ih - 0.5 * dens_h;
            double dux = dlx + dens_w;
            double duy = dly + dens_h;

            double charge_scale = scale_x * scale_y;
            if (b.type == BlockType::EDGE) {
                charge_scale *= rho_t_;
            }

            int i0 = clamp_idx((int)std::floor(dlx / bin_w_));
            int i1 = clamp_idx((int)std::ceil(dux / bin_w_));
            int j0 = clamp_idx((int)std::floor(dly / bin_h_));
            int j1 = clamp_idx((int)std::ceil(duy / bin_h_));

            double gx = 0.0, gy = 0.0;
            for (int j = j0; j < j1; ++j) {
                double by_lo = j * bin_h_;
                double by_hi = by_lo + bin_h_;
                double oy = std::max(0.0, std::min(duy, by_hi) - std::max(dly, by_lo));
                if (oy <= 0.0) continue;
                for (int i = i0; i < i1; ++i) {
                    double bx_lo = i * bin_w_;
                    double bx_hi = bx_lo + bin_w_;
                    double ox = std::max(0.0, std::min(dux, bx_hi) - std::max(dlx, bx_lo));
                    if (ox <= 0.0) continue;
                    double overlap = ox * oy * charge_scale;
                    gx += overlap * xi_x_[bin_idx(i, j)];
                    gy += overlap * xi_y_[bin_idx(i, j)];
                }
            }
            res.grad_x[idx] = gx;
            res.grad_y[idx] = gy;
        }
        return res;
    }

    // Test hooks.
    const std::vector<double>& rho_bins() const { return rho_bins_; }
    const std::vector<double>& xi_x() const { return xi_x_; }
    const std::vector<double>& xi_y() const { return xi_y_; }
    const std::vector<double>& psi() const { return psi_; }
    const std::vector<double>& halo() const { return halo_; }
    double rho_t() const { return rho_t_; }
    double bin_w() const { return bin_w_; }
    double bin_h() const { return bin_h_; }
    int grid() const { return grid_; }

private:
    int bin_idx(int i, int j) const { return j * grid_ + i; }

    int clamp_idx(int k) const {
        if (k < 0) return 0;
        if (k > grid_) return grid_;
        return k;
    }

    // Distribute charge_scale × area into all bins overlapped by [lx,ux]×[ly,uy].
    void deposit(double lx, double ly, double ux, double uy, double charge_scale) {
        if (charge_scale <= 0.0) return;
        int i0 = clamp_idx((int)std::floor(lx / bin_w_));
        int i1 = clamp_idx((int)std::ceil(ux / bin_w_));
        int j0 = clamp_idx((int)std::floor(ly / bin_h_));
        int j1 = clamp_idx((int)std::ceil(uy / bin_h_));
        for (int j = j0; j < j1; ++j) {
            double by_lo = j * bin_h_;
            double by_hi = by_lo + bin_h_;
            double oy = std::max(0.0, std::min(uy, by_hi) - std::max(ly, by_lo));
            if (oy <= 0.0) continue;
            for (int i = i0; i < i1; ++i) {
                double bx_lo = i * bin_w_;
                double bx_hi = bx_lo + bin_w_;
                double ox = std::max(0.0, std::min(ux, bx_hi) - std::max(lx, bx_lo));
                if (ox <= 0.0) continue;
                rho_bins_[bin_idx(i, j)] += ox * oy * charge_scale;
            }
        }
    }

    const Design& d_;
    int nb_;
    int grid_;
    double outline_w_ = 0.0, outline_h_ = 0.0;
    double bin_w_ = 0.0, bin_h_ = 0.0, bin_area_ = 0.0;
    double rho_t_ = 0.0;
    double movable_area_inflated_ = 0.0;

    std::vector<double> nets_total_;
    std::vector<double> halo_;

    std::vector<double> rho_bins_; // area-weighted charge per bin
    std::vector<double> psi_;
    std::vector<double> xi_x_;
    std::vector<double> xi_y_;

    std::unique_ptr<gpl::FFT> fft_;
};

} // namespace mp
