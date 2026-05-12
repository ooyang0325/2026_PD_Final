#pragma once
#include "types.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

// ─── CSV field splitter (handles quoted fields) ───────────────────────────────
static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    bool in_quote = false;
    for (char c : line) {
        if (c == '"') {
            in_quote = !in_quote;
        } else if (c == ',' && !in_quote) {
            fields.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    fields.push_back(cur);
    return fields;
}

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Parse a percentage like "20%" -> 0.2, or a decimal "0.2" -> 0.2
static double parse_pct(const std::string& s) {
    std::string t = trim(s);
    if (t.empty()) return 0.0;
    bool is_pct = (!t.empty() && t.back() == '%');
    if (is_pct) t.pop_back();
    double v = 0.0;
    try { v = std::stod(t); } catch (...) {}
    if (is_pct) v /= 100.0;
    return v;
}

class Parser {
public:
    // Detect if file is CSV (ends in .csv) or legacy .in format
    static Design load(const std::string& path) {
        std::string lp = path;
        std::transform(lp.begin(), lp.end(), lp.begin(), ::tolower);
        if (lp.size() >= 4 && lp.substr(lp.size()-4) == ".csv")
            return load_csv(path);
        return load_in(path);
    }

    static Design load_csv(const std::string& path) {
        Design d;
        std::ifstream f(path);
        if (!f) { std::cerr << "Cannot open " << path << "\n"; return d; }

        std::vector<std::vector<std::string>> rows;
        std::string line;
        while (std::getline(f, line)) {
            // strip \r
            if (!line.empty() && line.back() == '\r') line.pop_back();
            rows.push_back(split_csv(line));
        }

        enum Section { NONE, BLOCK, OUTLINE, ALPHA_SEC, CONN_HEADER, CONN_DATA } sec = NONE;
        std::vector<std::string> block_names;
        bool got_conn_col_header = false;

        for (auto& row : rows) {
            if (row.empty()) continue;
            std::string c0 = trim(row[0]);

            // Section detection
            if (c0 == "BLOCK") { sec = BLOCK; continue; }
            if (c0 == "OUTLINE") { sec = OUTLINE; continue; }

            // Alpha row: "alpha" or UTF-8 alpha character
            {
                std::string lc0 = c0;
                std::transform(lc0.begin(), lc0.end(), lc0.begin(), ::tolower);
                // strip utf-8 alpha character (α = 0xCE 0xB1 in utf-8)
                // also accept "a" or strings containing alpha
                bool is_alpha = (lc0 == "a" || lc0 == "alpha" ||
                                 lc0.find('\xce') != std::string::npos || // UTF-8 α
                                 lc0.find("alpha") != std::string::npos);
                if (is_alpha && row.size() > 1) {
                    std::string val = trim(row[1]);
                    if (!val.empty()) {
                        try { d.alpha = std::stod(val); } catch (...) {}
                    }
                    sec = NONE;
                    continue;
                }
            }

            // Connection matrix header detection
            if (c0.find("CONN") != std::string::npos ||
                c0.find("INTERFACE") != std::string::npos ||
                c0.find("MATRIX") != std::string::npos) {
                sec = CONN_HEADER;
                got_conn_col_header = false;
                continue;
            }

            // Skip blank/separator rows
            bool all_blank = true;
            for (auto& f2 : row) if (!trim(f2).empty()) { all_blank = false; break; }
            if (all_blank) continue;

            // ─── BLOCK section ───────────────────────────────────────────────
            if (sec == BLOCK) {
                // Skip header rows and tier label rows
                if (c0.empty() || c0 == "BLOCK" || c0 == "FT CONVERSION" ||
                    c0.substr(0,2) == "<=" || c0.substr(0,1) == ">") continue;
                if (c0.substr(0,3) != "BLK") continue;

                Block b;
                b.name = c0;

                // AREA
                b.area = 0.0;
                if (row.size() > 1 && !trim(row[1]).empty())
                    try { b.area = std::stod(trim(row[1])); } catch (...) {}

                // WIDTH, HEIGHT
                b.width = b.height = 0.0;
                if (row.size() > 2 && !trim(row[2]).empty())
                    try { b.width = std::stod(trim(row[2])); } catch (...) {}
                if (row.size() > 3 && !trim(row[3]).empty())
                    try { b.height = std::stod(trim(row[3])); } catch (...) {}

                // ASPECT RATIO RANGE
                b.min_ar = b.max_ar = 1.0;
                if (row.size() > 4 && !trim(row[4]).empty()) {
                    std::string ar_str = trim(row[4]);
                    auto comma = ar_str.find(',');
                    if (comma != std::string::npos) {
                        try { b.min_ar = std::stod(ar_str.substr(0, comma)); } catch (...) {}
                        try { b.max_ar = std::stod(ar_str.substr(comma+1)); } catch (...) {}
                    } else {
                        try {
                            double v = std::stod(ar_str);
                            b.min_ar = b.max_ar = v;
                        } catch (...) {}
                    }
                }

                // BLOCK TYPE
                b.type = BlockType::SOFT;
                b.has_fixed_wh = false;
                if (row.size() > 5 && !trim(row[5]).empty()) {
                    std::string bt = trim(row[5]);
                    if (bt == "EDGE")        { b.type = BlockType::EDGE;       b.has_fixed_wh = true; }
                    else if (bt == "MACRO")  { b.type = BlockType::HARD_MACRO; b.has_fixed_wh = true; }
                    else                     { b.type = BlockType::SOFT;       b.has_fixed_wh = false; }
                }

                // LOCATION
                if (row.size() > 6 && !trim(row[6]).empty()) {
                    std::string loc_str = trim(row[6]);
                    if (loc_str != "NONE") {
                        std::istringstream ls(loc_str);
                        std::string tok;
                        while (std::getline(ls, tok, ','))
                            if (!trim(tok).empty()) b.locations.push_back(trim(tok));
                    }
                }

                // FT CONVERSION (4 values, possibly percentages)
                double defaults[4] = {0.2, 0.4, 0.8, 1.0};
                for (int k = 0; k < 4; k++) {
                    b.ft.rate[k] = defaults[k];
                    if (row.size() > (size_t)(7+k) && !trim(row[7+k]).empty())
                        b.ft.rate[k] = parse_pct(row[7+k]);
                }

                // Default dimensions for SOFT blocks
                if (b.width == 0.0 && b.height == 0.0 && b.area > 0.0)
                    b.width = b.height = std::sqrt(b.area);
                else if (b.area == 0.0 && b.width > 0.0 && b.height > 0.0)
                    b.area = b.width * b.height;

                block_names.push_back(b.name);
                d.blocks.push_back(b);
                continue;
            }

            // ─── OUTLINE section ─────────────────────────────────────────────
            if (sec == OUTLINE) {
                if (c0 == "MAX") {
                    if (row.size() > 1 && !trim(row[1]).empty())
                        try { d.outline.max_width = std::stod(trim(row[1])); } catch (...) {}
                    if (row.size() > 2 && !trim(row[2]).empty())
                        try { d.outline.max_height = std::stod(trim(row[2])); } catch (...) {}
                    d.outline.cur_width  = d.outline.max_width;
                    d.outline.cur_height = d.outline.max_height;
                }
                continue;
            }

            // ─── CONNECTION section ──────────────────────────────────────────
            if (sec == CONN_HEADER) {
                // First non-empty row after header contains column names
                if (!got_conn_col_header) {
                    // The first field is blank (row-label column), rest are block names
                    // Verify they look like block names
                    bool has_blocks = false;
                    for (size_t i = 1; i < row.size(); i++)
                        if (trim(row[i]).substr(0,3) == "BLK") { has_blocks = true; break; }
                    if (has_blocks) {
                        got_conn_col_header = true;
                        // Trust the order from block_names vector, not from CSV header
                        // (block_names already parsed from BLOCK section)
                        sec = CONN_DATA;
                    }
                    continue;
                }
            }

            if (sec == CONN_DATA) {
                if (c0.substr(0,3) != "BLK") continue;
                std::string from_name = c0;
                int fi = -1;
                for (int i = 0; i < (int)block_names.size(); i++)
                    if (block_names[i] == from_name) { fi = i; break; }
                if (fi < 0) continue;

                for (int j = 0; j < (int)block_names.size(); j++) {
                    int nets = 0;
                    if ((size_t)(j+1) < row.size() && !trim(row[j+1]).empty()) {
                        try { nets = std::stoi(trim(row[j+1])); } catch (...) {}
                    }
                    if (nets <= 0) continue;

                    // Deduplicate: keep max for symmetric entries
                    int from_i = std::min(fi, j);
                    int to_i   = std::max(fi, j);
                    if (from_i == to_i) continue;

                    bool found = false;
                    for (auto& c : d.connections) {
                        if (c.from == from_i && c.to == to_i) {
                            c.nets = std::max(c.nets, nets);
                            found = true; break;
                        }
                    }
                    if (!found) {
                        Connection c;
                        c.from = from_i;
                        c.to   = to_i;
                        c.nets = nets;
                        d.connections.push_back(c);
                    }
                }
                continue;
            }
        }

        return d;
    }

    // ─── Legacy .in format parser (unchanged) ────────────────────────────────
    static Design load_in(const std::string& in_path) {
        Design d;
        std::ifstream f(in_path);
        if (!f) { std::cerr << "Cannot open " << in_path << "\n"; return d; }

        std::string line;
        enum Section { NONE, BLOCKS, CONNS } sec = NONE;
        std::map<std::string, int> name_to_idx;

        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string tag;
            ss >> tag;

            if (tag == "OUTLINE") {
                ss >> d.outline.max_width >> d.outline.max_height;
                d.outline.cur_width = d.outline.max_width;
                d.outline.cur_height = d.outline.max_height;
            } else if (tag == "ALPHA") {
                ss >> d.alpha;
            } else if (tag == "BLOCKS") {
                sec = BLOCKS;
            } else if (tag == "END_BLOCKS") {
                sec = NONE;
            } else if (tag == "CONNECTIONS") {
                sec = CONNS;
            } else if (tag == "END_CONNECTIONS") {
                sec = NONE;
            } else if (sec == BLOCKS && tag.substr(0,3) == "BLK") {
                Block b;
                b.name = tag;
                double area, w, h, min_ar, max_ar;
                std::string btype, loc_str;
                double ft0, ft1, ft2, ft3;
                ss >> area >> w >> h >> min_ar >> max_ar >> btype >> loc_str
                   >> ft0 >> ft1 >> ft2 >> ft3;
                b.area = area;
                b.width = w;
                b.height = h;
                b.min_ar = min_ar;
                b.max_ar = max_ar;
                b.has_fixed_wh = (btype != "SOFT");
                b.ft.rate[0] = ft0;
                b.ft.rate[1] = ft1;
                b.ft.rate[2] = ft2;
                b.ft.rate[3] = ft3;

                if (btype == "SOFT")       b.type = BlockType::SOFT;
                else if (btype == "MACRO") b.type = BlockType::HARD_MACRO;
                else                       b.type = BlockType::EDGE;

                if (loc_str != "NONE") {
                    std::istringstream ls(loc_str);
                    std::string tok;
                    while (std::getline(ls, tok, ','))
                        if (!tok.empty()) b.locations.push_back(tok);
                }

                name_to_idx[b.name] = (int)d.blocks.size();
                d.blocks.push_back(b);

            } else if (sec == CONNS) {
                std::string from_name, to_name;
                int nets;
                ss.str(line); ss.clear();
                ss >> from_name >> to_name >> nets;
                if (name_to_idx.count(from_name) && name_to_idx.count(to_name)) {
                    int fi = name_to_idx[from_name];
                    int ti = name_to_idx[to_name];
                    bool found = false;
                    for (auto& c : d.connections) {
                        if ((c.from == fi && c.to == ti) ||
                            (c.from == ti && c.to == fi)) {
                            c.nets = std::max(c.nets, nets);
                            found = true; break;
                        }
                    }
                    if (!found) {
                        Connection c;
                        c.from = std::min(fi, ti);
                        c.to   = std::max(fi, ti);
                        c.nets = nets;
                        d.connections.push_back(c);
                    }
                }
            }
        }
        return d;
    }
};
