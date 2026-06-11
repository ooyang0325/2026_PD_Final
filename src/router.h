#pragma once
#include "types.h"
#include "channel.h"
#include "config.h"
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <numeric>

// ─── Directed Face Graph ────────────────────────────────────────────────────
// Nodes: each rect (block or channel) × each edge × {entry, exit}
//   face_enter(r, e) = r*8 + (e-1)      [arriving at rect r from outside via edge e]
//   face_exit (r, e) = r*8 + (e-1) + 4  [leaving rect r to outside via edge e]
//
// Edges:
//   Intra-rect: face_enter(R, e_in) → face_exit(R, e_out)  [interior traversal]
//               only for channels and soft blocks (allow_ft=true)
//   Inter-rect: face_exit(A, e_out) → face_enter(B, e_in)  [boundary crossing, zero cost]
//               added for every pair of touching rects, both directions
//
// This directed model prevents "boundary relay" (chaining multiple adjacency edges
// without interior traversal), which would create invalid PATH entries.

struct FaceId {
    int rect_idx;
    int edge;   // 1=L, 2=T, 3=R, 4=B
};

struct RectInfo {
    std::string name;
    double lx, ly, w, h;
    bool is_channel;
    bool allow_ft;
};

static bool intervals_overlap(double a0, double a1, double b0, double b1) {
    return a0 < b1 - 1e-6 && b0 < a1 - 1e-6;
}

class GlobalRouter {
public:
    std::vector<RectInfo> rects;
    std::vector<double> ch_penalty;
    // Max FT nets each SOFT block can absorb by in-place growth into its
    // surrounding whitespace (optimistic, used only by the relief pass to
    // flag truly hopeless overloads).
    std::vector<double> ft_cap;
    // Blocks whose interior is priced out during the relief re-route.
    std::set<int> ft_relief;
    int n_blocks = 0;

    std::unordered_map<std::string, int> name_to_idx;

    void init(const std::vector<Block>& blocks,
              const std::vector<Channel>& channels,
              double ow, double oh) {
        rects.clear();
        name_to_idx.clear();
        n_blocks = (int)blocks.size();

        for (auto& b : blocks) {
            RectInfo r;
            r.name = b.name;
            r.lx = b.lx; r.ly = b.ly; r.w = b.width; r.h = b.height;
            r.is_channel = false;
            r.allow_ft = (b.type == BlockType::SOFT);
            name_to_idx[b.name] = (int)rects.size();
            rects.push_back(r);
        }
        for (auto& ch : channels) {
            RectInfo r;
            r.name = ch.name;
            r.lx = ch.lx; r.ly = ch.ly; r.w = ch.width; r.h = ch.height;
            r.is_channel = true;
            r.allow_ft = true;
            name_to_idx[ch.name] = (int)rects.size();
            rects.push_back(r);
        }
        ch_penalty.assign(channels.size(), 1.0);
        compute_ft_caps(blocks, ow, oh);
    }

    // Per-block absorbable FT load: how far can the block grow in place
    // (bounded by neighbors and the outline) and how many nets does that area
    // buy under the tiered conversion rates?  f(n) = (n/25)*rate(n)/2 is the
    // side extension a load of n demands; invert per tier, keep the largest
    // feasible candidate.
    void compute_ft_caps(const std::vector<Block>& blocks, double ow, double oh) {
        ft_cap.assign(blocks.size(), 0.0);
        for (int i = 0; i < (int)blocks.size(); i++) {
            const Block& b = blocks[i];
            if (b.type != BlockType::SOFT) continue;
            double L = b.lx, R = ow - (b.lx + b.width);
            double D = b.ly, U = oh - (b.ly + b.height);
            for (int j = 0; j < (int)blocks.size(); j++) {
                if (j == i) continue;
                const Block& o = blocks[j];
                bool yo = o.ly < b.ly + b.height - 1e-9 && o.ly + o.height > b.ly + 1e-9;
                bool xo = o.lx < b.lx + b.width  - 1e-9 && o.lx + o.width  > b.lx + 1e-9;
                if (yo) {
                    if (o.lx + o.width <= b.lx + 1e-9)      L = std::min(L, b.lx - (o.lx + o.width));
                    else if (o.lx >= b.lx + b.width - 1e-9) R = std::min(R, o.lx - (b.lx + b.width));
                }
                if (xo) {
                    if (o.ly + o.height <= b.ly + 1e-9)      D = std::min(D, b.ly - (o.ly + o.height));
                    else if (o.ly >= b.ly + b.height - 1e-9) U = std::min(U, o.ly - (b.ly + b.height));
                }
            }
            double maxW = b.width  + std::max(0.0, L) + std::max(0.0, R);
            double maxH = b.height + std::max(0.0, D) + std::max(0.0, U);
            double ext = std::sqrt(std::max(maxW * maxH, b.area)) - std::sqrt(b.area);

            const double* r = b.ft.rate;
            auto nets_at = [&](double rate) { return rate > 0 ? ext * 50.0 / rate : 1e18; };
            double cap = std::min(nets_at(r[0]), 3000.0);
            double n1 = nets_at(r[1]); if (n1 > 3000.0) cap = std::max(cap, std::min(n1, 6000.0));
            double n2 = nets_at(r[2]); if (n2 > 6000.0) cap = std::max(cap, std::min(n2, 9000.0));
            double n3 = nets_at(r[3]); if (n3 > 9000.0) cap = std::max(cap, n3);
            ft_cap[i] = cap;
        }
    }

    // Build directed face-graph with channel penalties for congestion avoidance.
    std::vector<std::vector<std::pair<int,double>>> build_adj() const {
        int N = (int)rects.size() * 8;
        std::vector<std::vector<std::pair<int,double>>> adj(N);

        // Intra-rect edges: enter(R,e_in) → exit(R,e_out) for traversable rects
        for (int i = 0; i < (int)rects.size(); i++) {
            auto& R = rects[i];
            if (!R.is_channel && !R.allow_ft) continue;
            for (int ein = 1; ein <= 4; ein++) {
                for (int eout = 1; eout <= 4; eout++) {
                    if (ein == eout) continue;
                    auto [inx, iny]   = ec(R, ein);
                    auto [outx, outy] = ec(R, eout);
                    double cost = std::abs(outx-inx) + std::abs(outy-iny);
                    if (R.is_channel) {
                        int ci = i - n_blocks;
                        if (ci >= 0 && ci < (int)ch_penalty.size())
                            cost *= ch_penalty[ci];
                    } else {
                        // Feeding through a SOFT block is expensive: it forces the
                        // module to grow (FT area conversion).  Charge a heavy
                        // multiplier so the router prefers routing AROUND through
                        // channels, and only feeds through as a last resort when
                        // no channel path exists.  Blocks under relief are priced
                        // out (connectivity preserved — no FAIL risk).
                        cost *= cfg::FT_TRAVERSE_PENALTY;
                        if (ft_relief.count(i)) cost *= 1e9;
                    }
                    adj[face_enter(i, ein)].push_back({face_exit(i, eout), cost});
                }
            }
        }

        // Inter-rect boundary edges: exit(A,e) → enter(B,e'), both directions
        for (int i = 0; i < (int)rects.size(); i++) {
            for (int j = i+1; j < (int)rects.size(); j++) {
                add_boundary_edges(i, rects[i], j, rects[j], adj);
            }
        }

        return adj;
    }

    // A routing demand (a connection, or a split fragment of one).
    struct Demand { int from, to; long nets; };

    // One rip-up/reroute episode over a fixed set of demands.  Routes each demand
    // on a single shortest path; the returned paths are aligned to `dems`.
    bool ripup_route(Design& d, const std::vector<Demand>& dems,
                     std::vector<RoutePath>& dem_paths, int max_rr) {
        std::fill(ch_penalty.begin(), ch_penalty.end(), 1.0);
        std::vector<int> order(dems.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b){ return dems[a].nets > dems[b].nets; });

        bool all_ok = false;
        for (int rr = 0; rr <= max_rr; rr++) {
            auto adj = build_adj();
            for (auto& ch : d.channels) { ch.nets_x = 0; ch.nets_y = 0; }
            dem_paths.assign(dems.size(), RoutePath{});
            all_ok = true;

            for (int k : order) {
                auto path = route_one(dems[k].from, dems[k].to, (int)dems[k].nets, adj, d);
                if (path.segments.empty()) all_ok = false;
                else { accum_nets(path, (int)dems[k].nets, d); dem_paths[k] = path; }
            }

            bool overflow = false;
            for (int i = 0; i < (int)d.channels.size(); i++) {
                if (d.channels[i].overflowed()) {
                    overflow = true;
                    ch_penalty[i] = std::min(ch_penalty[i] * 2.0, 1e6);
                }
            }
            if (!overflow && all_ok) break;
        }
        return all_ok;
    }

    static double total_overflow(const Design& d) {
        double o = 0;
        for (auto& ch : d.channels) {
            o += std::max(0.0, ch.nets_x - ch.cap_x());
            o += std::max(0.0, ch.nets_y - ch.cap_y());
        }
        return o;
    }

    // Hybrid net splitting (per the official Q&A: a connection's demand may be
    // split across multiple paths).  Round 0 routes every connection on a single
    // path (identical to the plain rip-up router).  After each round, any demand
    // whose path crosses an over-capacity channel is split in half into two
    // demands and the set is re-routed; the lowest-total-overflow result across
    // all rounds is kept.  Splitting therefore only ever happens for connections
    // involved in an overflow, and can never make the result worse than the
    // single-path solution (round 0 is always a candidate).
    bool route_all(Design& d, int max_rr = 8) {
        std::vector<Demand> dems;
        dems.reserve(d.connections.size());
        for (auto& c : d.connections) dems.push_back({c.from, c.to, (long)c.nets});

        const int MAX_ROUNDS = 6;
        std::vector<RoutePath> best_paths;
        double best_ovf = 1e30;
        bool best_ok = false;

        std::vector<RoutePath> dem_paths;
        for (int round = 0; round < MAX_ROUNDS; round++) {
            bool ok = ripup_route(d, dems, dem_paths, max_rr);
            double ovf = total_overflow(d);
            if (ovf < best_ovf - 1e-6) {
                best_ovf = ovf;
                best_paths.clear();
                for (auto& p : dem_paths) if (!p.segments.empty()) best_paths.push_back(p);
                best_ok = ok;
            }
            if (ovf < 1e-6) break;

            // Split demands whose path crosses an over-capacity channel.
            std::set<std::string> over;
            for (auto& ch : d.channels) if (ch.overflowed()) over.insert(ch.name);

            std::vector<Demand> next;
            next.reserve(dems.size() * 2);
            bool split_any = false;
            for (size_t k = 0; k < dems.size(); k++) {
                bool crosses = false;
                if (dems[k].nets > 1)
                    for (auto& seg : dem_paths[k].segments)
                        if (over.count(seg.rect_name)) { crosses = true; break; }
                if (crosses) {
                    long h = dems[k].nets / 2;
                    next.push_back({dems[k].from, dems[k].to, h});
                    next.push_back({dems[k].from, dems[k].to, dems[k].nets - h});
                    split_any = true;
                } else {
                    next.push_back(dems[k]);
                }
            }
            if (!split_any) break;
            dems.swap(next);
        }

        // Surgical artery relief, then adopt the routing and recompute usage.
        relief_pass(d, best_paths);
        d.paths = best_paths;
        for (auto& ch : d.channels) { ch.nets_x = 0; ch.nets_y = 0; }
        for (auto& p : d.paths) accum_nets(p, p.nets, d);
        return best_ok;
    }

    // ─── Artery relief pass ─────────────────────────────────────────────────
    // One-shot and surgical: blocks whose final FT load exceeds what their
    // surrounding whitespace can absorb (ft_cap) get their interiors priced
    // out, and ONLY the paths crossing them are re-routed.  The result is
    // kept only if the penalty-aligned badness improves.  Unlike in-anneal
    // load feedback (which measurably oscillated and regressed 3-5x), this
    // cannot make the routing worse.
    void relief_pass(Design& d, std::vector<RoutePath>& paths) {
        auto load = soft_loads(paths);
        std::set<std::string> over;
        for (int i = 0; i < n_blocks; i++)
            if (load[i] > ft_cap[i] + 1e-9) { over.insert(rects[i].name); ft_relief.insert(i); }
        if (over.empty()) return;

        auto recompute = [&](std::vector<RoutePath>& ps) {
            for (auto& ch : d.channels) { ch.nets_x = 0; ch.nets_y = 0; }
            for (auto& p : ps) accum_nets(p, p.nets, d);
        };
        recompute(paths);
        double before = routing_badness(d, paths);

        auto adj = build_adj();   // with relief pricing
        ft_relief.clear();

        std::vector<RoutePath> cand = paths;
        for (auto& p : cand) {
            if (p.segments.size() < 3) continue;
            bool crosses = false;
            for (size_t k = 1; k + 1 < p.segments.size(); k++)
                if (over.count(p.segments[k].rect_name)) { crosses = true; break; }
            if (!crosses) continue;
            auto fi = name_to_idx.find(p.segments.front().rect_name);
            auto ti = name_to_idx.find(p.segments.back().rect_name);
            if (fi == name_to_idx.end() || ti == name_to_idx.end()) continue;
            auto np = route_one(fi->second, ti->second, (int)p.nets, adj, d);
            if (!np.segments.empty()) p = np;
        }
        recompute(cand);
        double after = routing_badness(d, cand);
        if (after < before - 1e-6) paths.swap(cand);
        recompute(paths);
    }

    // Per-SOFT-block feedthrough nets implied by a set of paths.
    std::vector<double> soft_loads(const std::vector<RoutePath>& paths) const {
        std::vector<double> load(n_blocks, 0.0);
        for (auto& p : paths) {
            for (size_t k = 1; k + 1 < p.segments.size(); k++) {
                auto it = name_to_idx.find(p.segments[k].rect_name);
                if (it != name_to_idx.end() && it->second < n_blocks &&
                    rects[it->second].allow_ft)
                    load[it->second] += p.nets;
            }
        }
        return load;
    }

    // Penalty-count-aligned routing badness: the evaluator charges one penalty
    // per overflowing channel-direction and one per unsatisfiable soft block,
    // independent of magnitude — count first, excess nets as tie-break.
    double routing_badness(const Design& d, const std::vector<RoutePath>& paths) const {
        int cnt = 0;
        double mag = 0;
        for (auto& ch : d.channels) {
            double ox = ch.nets_x - ch.cap_x(), oy = ch.nets_y - ch.cap_y();
            if (ox > 1e-9) { cnt++; mag += ox; }
            if (oy > 1e-9) { cnt++; mag += oy; }
        }
        auto load = soft_loads(paths);
        for (int i = 0; i < n_blocks; i++)
            if (load[i] > ft_cap[i] + 1e-9) { cnt++; mag += load[i] - ft_cap[i]; }
        return cnt * 1e7 + mag;
    }

    // Dijkstra over face-graph; targets any entry face of destination.
    RoutePath route_one(int from_blk, int to_blk, int nets,
                        const std::vector<std::vector<std::pair<int,double>>>& adj,
                        const Design& d) {
        int NF = (int)rects.size() * 8;
        std::vector<double> dist(NF, 1e18);
        std::vector<int> prev(NF, -1);
        using PD = std::pair<double,int>;
        std::priority_queue<PD, std::vector<PD>, std::greater<PD>> pq;

        // Source: start from all exit-faces of source block
        for (int e = 1; e <= 4; e++) {
            int f = face_exit(from_blk, e);
            dist[f] = 0;
            pq.push({0, f});
        }

        while (!pq.empty()) {
            auto [d_cur, f_cur] = pq.top(); pq.pop();
            if (d_cur > dist[f_cur] + 1e-9) continue;

            // Check if we reached any entry-face of destination
            int ri = rect_of(f_cur);
            bool is_enter = is_entry(f_cur);
            if (is_enter && ri == to_blk && ri != from_blk) {
                return reconstruct(from_blk, to_blk, f_cur, prev, nets);
            }

            for (auto& [fn, cost] : adj[f_cur]) {
                double nd = d_cur + cost;
                if (nd < dist[fn]) {
                    dist[fn] = nd;
                    prev[fn] = f_cur;
                    pq.push({nd, fn});
                }
            }
        }

        RoutePath empty; empty.nets = nets;
        return empty;
    }

private:
    // Directed face encoding: 8 nodes per rect
    static int face_enter(int r, int e) { return r * 8 + (e - 1); }
    static int face_exit (int r, int e) { return r * 8 + (e - 1) + 4; }
    static int rect_of(int f) { return f / 8; }
    static int edge_of(int f) { return (f % 4) + 1; }
    static bool is_entry(int f) { return (f % 8) < 4; }

    static std::pair<double,double> ec(const RectInfo& R, int edge) {
        switch (edge) {
            case 1: return {R.lx,         R.ly + R.h/2};
            case 2: return {R.lx + R.w/2, R.ly + R.h};
            case 3: return {R.lx + R.w,   R.ly + R.h/2};
            case 4: return {R.lx + R.w/2, R.ly};
        }
        return {0, 0};
    }

    void add_boundary_edges(int i, const RectInfo& A,
                            int j, const RectInfo& B,
                            std::vector<std::vector<std::pair<int,double>>>& adj) const {
        // A.right touches B.left
        if (std::abs((A.lx + A.w) - B.lx) < 1e-3 &&
            intervals_overlap(A.ly, A.ly+A.h, B.ly, B.ly+B.h)) {
            adj[face_exit(i, 3)].push_back({face_enter(j, 1), 0.0});
            adj[face_exit(j, 1)].push_back({face_enter(i, 3), 0.0});
        }
        // B.right touches A.left
        if (std::abs((B.lx + B.w) - A.lx) < 1e-3 &&
            intervals_overlap(A.ly, A.ly+A.h, B.ly, B.ly+B.h)) {
            adj[face_exit(j, 3)].push_back({face_enter(i, 1), 0.0});
            adj[face_exit(i, 1)].push_back({face_enter(j, 3), 0.0});
        }
        // A.top touches B.bottom
        if (std::abs((A.ly + A.h) - B.ly) < 1e-3 &&
            intervals_overlap(A.lx, A.lx+A.w, B.lx, B.lx+B.w)) {
            adj[face_exit(i, 2)].push_back({face_enter(j, 4), 0.0});
            adj[face_exit(j, 4)].push_back({face_enter(i, 2), 0.0});
        }
        // B.top touches A.bottom
        if (std::abs((B.ly + B.h) - A.ly) < 1e-3 &&
            intervals_overlap(A.lx, A.lx+A.w, B.lx, B.lx+B.w)) {
            adj[face_exit(j, 2)].push_back({face_enter(i, 4), 0.0});
            adj[face_exit(i, 4)].push_back({face_enter(j, 2), 0.0});
        }
    }

    static double seg_overlap(double a0, double a1, double b0, double b1) {
        return std::max(0.0, std::min(a1,b1) - std::max(a0,b0));
    }

    void accum_nets(const RoutePath& path, int nets, Design& d) const {
        for (auto& seg : path.segments) {
            if (seg.rect_name.size() >= 2 && seg.rect_name.substr(0,2) == "CH") {
                for (auto& ch : d.channels) {
                    if (ch.name == seg.rect_name) {
                        // 判定跨越方向
                        bool has_x = (seg.edge_in == 1 || seg.edge_in == 3 || seg.edge_out == 1 || seg.edge_out == 3);
                        bool has_y = (seg.edge_in == 2 || seg.edge_in == 4 || seg.edge_out == 2 || seg.edge_out == 4);
                        if (has_x) const_cast<Channel&>(ch).nets_x += nets;
                        if (has_y) const_cast<Channel&>(ch).nets_y += nets;
                        break;
                    }
                }
            }
        }
    }

    // Convert face sequence into rect-edge segments for output formatting.
    RoutePath reconstruct(int from_blk, int to_blk, int f_end,
                          const std::vector<int>& prev, int nets) const {
        // Trace back: the chain alternates exit/enter nodes.
        // Pattern: exit(SRC) → enter(R1) → exit(R1) → enter(R2) → ... → enter(DST)
        std::vector<int> faces;
        int f = f_end;
        int guard = 200000;
        while (f != -1 && --guard > 0) {
            faces.push_back(f);
            f = prev[f];
        }
        std::reverse(faces.begin(), faces.end());

        RoutePath path;
        path.nets = nets;
        if (faces.empty()) return path;

        // Group consecutive same-rect faces into segments.
        // Each rect appears as at most two faces: exit (source), enter+exit (intermediate), enter (dest).
        int cur_rect = -1;
        for (int fi : faces) {
            int ri = rect_of(fi);
            int ei = edge_of(fi);
            bool entering = is_entry(fi);

            if (ri != cur_rect) {
                PathSegment seg;
                seg.rect_name = rects[ri].name;
                if (entering) {
                    seg.edge_in  = ei;
                    seg.edge_out = 0; // will be filled by exit face or left 0 for dest
                } else {
                    // exit face (source or intermediate exit)
                    seg.edge_in  = 0; // will be filled by prev enter or left 0 for source
                    seg.edge_out = ei;
                }
                path.segments.push_back(seg);
                cur_rect = ri;
            } else {
                // Second face of same rect
                if (entering) {
                    path.segments.back().edge_in = ei;
                } else {
                    path.segments.back().edge_out = ei;
                }
            }
        }

        // Validate: source should have edge_in=0 (no entry), dest should have edge_out=0
        if (!path.segments.empty()) path.segments.front().edge_in = 0;
        if (!path.segments.empty()) path.segments.back().edge_out = 0;

        return path;
    }
};
