#pragma once
#include "types.h"
#include "floorplan.h"
#include <random>
#include <algorithm>
#include <string>

// Phase 1 stub for the mathematical-programming (ePlace-style) engine.
// Just scatters blocks at seeded-random positions inside the max outline.
// Subsequent phases will replace run() with WA + density + Nesterov.
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

        int nb = (int)d_.blocks.size();
        for (int i = 0; i < nb; i++) {
            auto& b = d_.blocks[i];

            // Phase 1: leave width/height alone.  HARD already has fixed dims,
            // SOFT keeps current dims (sized later by the routability loop),
            // EDGE has been parsed with valid dims.
            double w = b.width, h = b.height;

            double maxx = std::max(0.0, W - w);
            double maxy = std::max(0.0, H - h);

            if (b.type == BlockType::EDGE && !b.locations.empty()) {
                // Pick a random allowed location string (T/B/L/R combination).
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
                b.lx = x;
                b.ly = y;
            } else {
                b.lx = u01(rng_) * maxx;
                b.ly = u01(rng_) * maxy;
            }

            // Keep Floorplan's parallel W/H arrays consistent so downstream
            // helpers (channels, finalize) see matching dims.
            fp_.W[i] = w;
            fp_.H[i] = h;
        }
    }

private:
    Floorplan& fp_;
    Design&    d_;
    std::mt19937 rng_;
};
