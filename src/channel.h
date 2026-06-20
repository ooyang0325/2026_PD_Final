#pragma once
#include "types.h"
#include <vector>
#include <set>
#include <algorithm>
#include <sstream>

// Channel calculation: vertical-strip based decomposition.
// For each pair of consecutive x-coordinates, find y-intervals not covered by blocks.
// Each uncovered rectangle in a vertical strip becomes one channel.

class ChannelCalculator {
public:
    static std::vector<Channel> compute(const std::vector<Block>& blocks,
                                         double outline_w, double outline_h) {
        std::set<double> xs_set;
        xs_set.insert(0.0);
        xs_set.insert(outline_w);
        for (auto& b : blocks) {
            xs_set.insert(b.lx);
            xs_set.insert(b.lx + b.width);
        }
        std::vector<double> xs(xs_set.begin(), xs_set.end());
        std::sort(xs.begin(), xs.end());

        std::vector<Channel> channels;
        int ch_id = 0;

        for (int k = 0; k + 1 < (int)xs.size(); k++) {
            double x0 = xs[k], x1 = xs[k+1];
            if (x1 - x0 < 1e-9) continue;

            // Find blocks that overlap this vertical strip (x0,x1)
            std::vector<std::pair<double,double>> covered;
            for (auto& b : blocks) {
                double bx0 = b.lx, bx1 = b.lx + b.width;
                if (bx0 < x1 - 1e-9 && bx1 > x0 + 1e-9) {
                    covered.push_back({b.ly, b.ly + b.height});
                }
            }

            // Invert covered y-intervals to produce channel rectangles
            auto free = uncovered_intervals(covered, 0.0, outline_h);
            for (auto& [y0, y1] : free) {
                if (y1 - y0 < 1e-9) continue;
                Channel ch;
                ch.name = "CH" + std::to_string(ch_id++);
                ch.lx = x0;
                ch.ly = y0;
                ch.width = x1 - x0;
                ch.height = y1 - y0;
                channels.push_back(ch);
            }
        }
        return channels;
    }

private:
    static std::vector<std::pair<double,double>>
    uncovered_intervals(std::vector<std::pair<double,double>> covered,
                        double lo, double hi) {
        // Sort and merge covered intervals, then invert into free intervals
        std::sort(covered.begin(), covered.end());
        std::vector<std::pair<double,double>> merged;
        for (auto& iv : covered) {
            double a = std::max(iv.first, lo);
            double b = std::min(iv.second, hi);
            if (b <= a + 1e-9) continue;
            if (!merged.empty() && a <= merged.back().second + 1e-9) {
                merged.back().second = std::max(merged.back().second, b);
            } else {
                merged.push_back({a, b});
            }
        }
        std::vector<std::pair<double,double>> free;
        double cur = lo;
        for (auto& iv : merged) {
            if (iv.first - cur > 1e-9) free.push_back({cur, iv.first});
            cur = std::max(cur, iv.second);
        }
        if (hi - cur > 1e-9) free.push_back({cur, hi});
        return free;
    }
};
