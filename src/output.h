#pragma once
#include "types.h"
#include "channel.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

class OutputWriter {
public:
    // Write the .cfg output file
    static void write(const Design& d, const std::string& out_path) {
        std::ofstream f(out_path);
        if (!f) { std::cerr << "Cannot write " << out_path << "\n"; return; }

        // Outline: current packed width/height
        f << "OUTLINE\n";
        f << fmt(d.outline.cur_width) << " " << fmt(d.outline.cur_height) << "\n";
        f << "END\n\n";

        // Blocks: lower-left coords and dimensions
        f << "BLOCK\n";
        for (auto& b : d.blocks) {
            f << b.name << " "
              << fmt(b.lx) << " " << fmt(b.ly) << " "
              << fmt(b.width) << " " << fmt(b.height) << "\n";
        }
        f << "END\n\n";

        // Channels: vertical-strip rectangles for routing
        f << "CHANNEL\n";
        for (auto& ch : d.channels) {
            f << ch.name << " "
              << fmt(ch.lx) << " " << fmt(ch.ly) << " "
              << fmt(ch.width) << " " << fmt(ch.height) << "\n";
        }
        f << "END\n\n";

        // Paths: sequence of (rect, edge) entries per spec
        f << "PATH\n";
        for (auto& path : d.paths) {
            if (path.segments.empty()) continue;
            f << "PATH " << path.nets;
            for (auto& seg : path.segments) {
                if (seg.edge_in == 0 && seg.edge_out != 0) {
                    // Source block
                    f << " " << seg.rect_name << " " << seg.edge_out;
                } else if (seg.edge_out == 0) {
                    // Sink block
                    f << " " << seg.rect_name << " " << seg.edge_in;
                } else {
                    // Intermediate: enter via edge_in, exit via edge_out
                    f << " " << seg.rect_name << " " << seg.edge_in
                      << " " << seg.rect_name << " " << seg.edge_out;
                }
            }
            f << "\n";
        }
        f << "END\n";

        f.close();
        std::cerr << "[Output] Written to " << out_path << "\n";
    }

    // Print summary to stderr
    static void print_summary(const Design& d) {
        double area = d.outline.cur_width * d.outline.cur_height;
        // Delegate HPWL and cost calculation to external consistent logic or re-call get_guiding_point here
        // For brevity, print overflow details
        std::cerr << "=== Summary ===\n"
                  << "Outline: " << fmt(d.outline.cur_width) << " x " << fmt(d.outline.cur_height) << "\n"
                  << "Area: " << fmt(area) << "\n";

        int n_overflow = 0;
        for (auto& ch : d.channels) {
            if (ch.overflowed()) {
                n_overflow++;
                std::cerr << "  [Overflow] " << ch.name << " Nets X:" << ch.nets_x << "/" << ch.cap_x() 
                          << " Y:" << ch.nets_y << "/" << ch.cap_y() << "\n";
            }
        }
        std::cerr << "Channel overflow: " << n_overflow << " / " << d.channels.size() << "\n";
    }

private:
    static std::string fmt(double v) {
        // Round to 2 decimal places
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << v;
        return ss.str();
    }
};
