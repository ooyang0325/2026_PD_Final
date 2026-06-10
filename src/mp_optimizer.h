#pragma once
#include "types.h"
#include "config.h"
#include "floorplan.h"
#include "mp/density.h"
#include "mp/wa_wirelength.h"
#include "mp/nesterov.h"
#include <random>
#include <algorithm>
#include <string>
#include <vector>

// Phase 4: replace the random-scatter stub with the ePlace-style Nesterov
// solver.  Seeded random init still produces x0/y0; Nesterov::run then
// descends to the configured overflow target.  CG-legalize / routability
// loop arrive in Phases 5 and 6 (see docs/plans/2026-06-10-mp-floorplanner.md).
class MPOptimizer {
public:
    MPOptimizer(Floorplan& fp, Design& d, unsigned seed)
        : fp_(fp), d_(d), rng_(seed) {}

    void run() {
        const double W = d_.outline.max_width;
        const double H = d_.outline.max_height;
        d_.outline.cur_width  = W;
        d_.outline.cur_height = H;

        std::uniform_real_distribution<double> u01(0.0, 1.0);
        const int nb = static_cast<int>(d_.blocks.size());

        // ---- Seeded random scatter init (this becomes x0/y0) ----
        std::vector<double> x0(nb), y0(nb);
        for (int i = 0; i < nb; ++i) {
            auto& b = d_.blocks[i];
            double w = b.width, h = b.height;
            double maxx = std::max(0.0, W - w);
            double maxy = std::max(0.0, H - h);

            if (b.type == BlockType::EDGE && !b.locations.empty()) {
                std::uniform_int_distribution<int> pick(0, (int)b.locations.size() - 1);
                int li = pick(rng_);
                fp_.active_loc[i] = li;
                const std::string& loc = b.locations[li];

                double x = u01(rng_) * maxx;
                double y = u01(rng_) * maxy;
                if (loc.find('L') != std::string::npos) x = 0.0;
                if (loc.find('R') != std::string::npos) x = maxx;
                if (loc.find('B') != std::string::npos) y = 0.0;
                if (loc.find('T') != std::string::npos) y = maxy;
                x0[i] = x;
                y0[i] = y;
            } else {
                x0[i] = u01(rng_) * maxx;
                y0[i] = u01(rng_) * maxy;
            }
            fp_.W[i] = w;
            fp_.H[i] = h;
        }

        // ---- Nesterov global place ----
        mp::Density       dens(d_, cfg::MP_GRID);
        mp::WAWirelength  wl(d_);
        mp::NesterovParams params;
        params.max_iter        = cfg::MP_MAX_ITER;
        params.target_overflow = cfg::MP_TARGET_OVF;
        params.init_lambda     = cfg::MP_INIT_LAMBDA;
        params.phi_min         = cfg::MP_PHI_MIN;
        params.phi_max         = cfg::MP_PHI_MAX;
        mp::Nesterov nes(dens, wl, d_, params);
        auto result = nes.run(x0, y0);

        // Copy result back into the Design.  Even on divergence the returned
        // snapshot is the best-tau point encountered (guard #7).
        for (int i = 0; i < nb; ++i) {
            d_.blocks[i].lx = result.x[i];
            d_.blocks[i].ly = result.y[i];
        }
    }

private:
    Floorplan& fp_;
    Design&    d_;
    std::mt19937 rng_;
};
