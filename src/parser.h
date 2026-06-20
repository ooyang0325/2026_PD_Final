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

namespace {

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Parse percentage ("20%" → 0.2) or bare decimal ("0.2" → 0.2).
double parse_pct(const std::string& s) {
    std::string t = trim(s);
    if (t.empty()) return 0.0;
    bool is_pct = (!t.empty() && t.back() == '%');
    if (is_pct) t.pop_back();
    double v = 0.0;
    try { v = std::stod(t); } catch (...) {}
    if (is_pct) v /= 100.0;
    return v;
}

} // namespace

class Parser {
public:
    static Design load_csv(const std::string& path) {
        Design d;
        std::ifstream f(path, std::ios::binary);
        if (!f) { std::cerr << "Cannot open " << path << "\n"; return d; }

        // Read entire file into a string.
        std::string content((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

        // Strip UTF-8 BOM if present.
        if (content.size() >= 3 &&
            (unsigned char)content[0] == 0xEF &&
            (unsigned char)content[1] == 0xBB &&
            (unsigned char)content[2] == 0xBF) {
            content = content.substr(3);
        }

        // Robust CSV parse: supports quoted fields (including embedded newlines).
        std::vector<std::vector<std::string>> rows;
        std::vector<std::string> current_row;
        std::string field;
        bool in_quote = false;

        for (size_t i = 0; i < content.size(); i++) {
            char c = content[i];
            if (c == '"') {
                if (in_quote && i + 1 < content.size() && content[i+1] == '"') {
                    field += '"'; // escaped quote ""
                    i++;
                } else {
                    in_quote = !in_quote;
                }
            } else if (c == ',' && !in_quote) {
                current_row.push_back(field);
                field.clear();
            } else if ((c == '\n' || c == '\r') && !in_quote) {
                if (c == '\r' && i + 1 < content.size() && content[i+1] == '\n') {
                    i++; // consume CRLF as one line ending
                }
                current_row.push_back(field);
                rows.push_back(current_row);
                current_row.clear();
                field.clear();
            } else {
                field += c;
            }
        }
        if (!field.empty() || !current_row.empty()) {
            current_row.push_back(field);
            rows.push_back(current_row);
        }

        enum Section { NONE, BLOCK, OUTLINE, CONN_HEADER, CONN_DATA } sec = NONE;
        std::vector<std::string> block_names;
        bool got_conn_col_header = false;
        // The 2026-06-17 testcase format inserts a "PORT EDGE" column between
        // LOCATION (col 6) and FT CONVERSION (cols 7+).  Detect it from the BLOCK
        // header so FT rates are read from the correct columns (backward compatible
        // with the old format, which has no PORT EDGE column).
        bool has_port_edge = false;

        for (auto& row : rows) {
            if (row.empty()) continue;
            std::string c0 = trim(row[0]);

            // Section detection.
            if (c0 == "BLOCK") {
                sec = BLOCK;
                for (auto& cell : row) {
                    std::string t = trim(cell);
                    std::transform(t.begin(), t.end(), t.begin(), ::toupper);
                    if (t.find("PORT EDGE") != std::string::npos) { has_port_edge = true; break; }
                }
                continue;
            }
            if (c0 == "OUTLINE") { sec = OUTLINE; continue; }

            // Detect alpha row (compatible with various UTF-8 representations of α).
            {
                std::string lc0 = c0;
                std::transform(lc0.begin(), lc0.end(), lc0.begin(), ::tolower);
                bool is_alpha = (lc0 == "a" || lc0 == "alpha" ||
                                 lc0.find("α") != std::string::npos ||
                                 lc0.find('\xce') != std::string::npos);
                if (is_alpha && row.size() > 1) {
                    std::string val = trim(row[1]);
                    if (!val.empty()) {
                        try { d.alpha = std::stod(val); } catch (...) {}
                    }
                    sec = NONE;
                    continue;
                }
            }

            if (c0.find("CONN") != std::string::npos ||
                c0.find("INTERFACE") != std::string::npos ||
                c0.find("MATRIX") != std::string::npos) {
                sec = CONN_HEADER;
                got_conn_col_header = false;
                continue;
            }

            // Skip blank rows.
            bool all_blank = true;
            for (auto& f2 : row) if (!trim(f2).empty()) { all_blank = false; break; }
            if (all_blank) continue;

            // ─── BLOCK section ────────────────────────────────────────────────
            if (sec == BLOCK) {
                if (c0.empty() || c0 == "BLOCK" || c0 == "FT CONVERSION" ||
                    c0.substr(0,2) == "<=" || c0.substr(0,1) == ">" || c0 == "AREA") continue;
                if (c0.substr(0,3) != "BLK") continue;

                Block b;
                b.name = c0;

                b.area = 0.0;
                if (row.size() > 1 && !trim(row[1]).empty())
                    try { b.area = std::stod(trim(row[1])); } catch (...) {}

                b.width = b.height = 0.0;
                if (row.size() > 2 && !trim(row[2]).empty())
                    try { b.width = std::stod(trim(row[2])); } catch (...) {}
                if (row.size() > 3 && !trim(row[3]).empty())
                    try { b.height = std::stod(trim(row[3])); } catch (...) {}

                // Aspect ratio: supports "0.5,2" range or a single fixed value.
                b.min_ar = b.max_ar = 1.0;
                if (row.size() > 4 && !trim(row[4]).empty()) {
                    std::string ar_str = trim(row[4]);
                    ar_str.erase(std::remove(ar_str.begin(), ar_str.end(), '"'), ar_str.end());
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

                // Block type from the LOCATION column keyword.
                b.type = BlockType::SOFT;
                b.has_fixed_wh = false;
                if (row.size() > 5 && !trim(row[5]).empty()) {
                    std::string bt = trim(row[5]);
                    if (bt.find("EDGE") != std::string::npos)        { b.type = BlockType::EDGE;       b.has_fixed_wh = true; }
                    else if (bt.find("MACRO") != std::string::npos)  { b.type = BlockType::HARD_MACRO; b.has_fixed_wh = true; }
                    else                                             { b.type = BlockType::SOFT;       b.has_fixed_wh = false; }
                }

                // LOCATION constraint (supports comma-separated codes like "BR,RB").
                if (row.size() > 6 && !trim(row[6]).empty()) {
                    std::string loc_str = trim(row[6]);
                    loc_str.erase(std::remove(loc_str.begin(), loc_str.end(), '"'), loc_str.end());
                    if (loc_str != "NONE" && loc_str != "") {
                        std::istringstream ls(loc_str);
                        std::string tok;
                        while (std::getline(ls, tok, ',')) {
                            std::string t = trim(tok);
                            if (!t.empty()) b.locations.push_back(t);
                        }
                    }
                }

                // PORT EDGE (col 7, new format only): 1=left,2=top,3=right,4=bottom; 0/empty=unrestricted.
                b.port_edge = 0;
                if (has_port_edge && row.size() > 7 && !trim(row[7]).empty()) {
                    try {
                        int pe = std::stoi(trim(row[7]));
                        if (pe >= 1 && pe <= 4) b.port_edge = pe;
                    } catch (...) {}
                }

                // FT conversion rates (columns shift right by 1 when PORT EDGE is present).
                int ftc = has_port_edge ? 8 : 7;
                double defaults[4] = {0.2, 0.4, 0.8, 1.0};
                for (int k = 0; k < 4; k++) {
                    b.ft.rate[k] = defaults[k];
                    if (row.size() > (size_t)(ftc+k) && !trim(row[ftc+k]).empty())
                        b.ft.rate[k] = parse_pct(row[ftc+k]);
                }

                // Derive missing dimension: if only area given, assume square.
                if (b.width == 0.0 && b.height == 0.0 && b.area > 0.0) {
                    b.width = std::sqrt(b.area);
                    b.height = std::sqrt(b.area);
                } else if (b.area == 0.0 && b.width > 0.0 && b.height > 0.0) {
                    b.area = b.width * b.height;
                }

                block_names.push_back(b.name);
                d.blocks.push_back(b);
                continue;
            }

            // ─── OUTLINE section ──────────────────────────────────────────────
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

            // ─── CONNECTION MATRIX section ────────────────────────────────────
            if (sec == CONN_HEADER) {
                if (!got_conn_col_header) {
                    bool has_blocks = false;
                    for (size_t i = 1; i < row.size(); i++)
                        if (trim(row[i]).substr(0,3) == "BLK") { has_blocks = true; break; }
                    if (has_blocks) {
                        got_conn_col_header = true;
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

};
