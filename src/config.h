#pragma once
#include <cstdlib>

// Global tunable parameters.  Defaults are the production values; each can be
// overridden via an environment variable for offline parameter sweeps.  Values
// are read once at program start (cfg::load_from_env) and treated as read-only
// thereafter, so they are safe to share across worker threads.
namespace cfg {

inline double HALO               = 30.0; // routing-channel gap around every block
inline double FT_TRAVERSE_PENALTY = 3.0; // router cost multiplier for feedthrough
inline double FTW                = 0.4; // SA feedthrough-minimizing penalty weight
inline double PCONGW             = 1.5; // SA placement-congestion weight (routing-aware)
inline double PCONG_THRESH       = 500; // only high-demand connections drive congestion

// SA feedthrough-concentration weight: penalizes single soft blocks whose
// estimated through-traffic exceeds FT_EST_CAP (the top of the cheapest
// conversion-rate tier).  Targets "artery" blocks the routed FT piles onto.
inline double FTCW               = 3.0;
inline double FT_EST_CAP         = 3000.0;

// SA cut-overflow weight: penalizes the worst vertical/horizontal cut whose
// straddling net demand exceeds the cut's carrying capacity (25 nets per um
// of whitespace + bounded FT credit through bridging soft blocks).  Demand
// above capacity at a cut is unroutable by ANY router — only moving the
// connected blocks to the same side fixes it (the b50u65 "moat" failure).
inline double MOATW              = 10.0;

// Routed-feedback weight: at checkpoints during the anneal the REAL router
// runs on the best layout so far, and soft blocks it overloads beyond their
// whitespace-absorbable capacity become positional repulsion marks in the SA
// cost.  This is the only artery signal not filtered through the line-of-
// sight estimate (which under-counts routed wandering).  0 disables.
inline double RFBW               = 6.0;

// Route-driven legalize loop (post-SA refinement).
inline int    LEG_ITERS          = 128;   // max outer iterations
inline double LEG_TIME           = 30.0; // wall-clock cap (seconds)
inline int    LEG_TRIES          = 8;    // displace mini-SA tries per attempt
inline int    LEG_ENABLE         = 1;    // 0 = skip the loop (use legacy ft_iter)

// Final-pass analytical (force-directed) legalizer.  Breaks the B*-tree and
// outputs continuous coordinates.  Strict rollback if it doesn't improve.
inline int    ANA_ENABLE         = 1;    // 0 = skip the analytical pass
inline int    ANA_ITERS          = 128;   // force-directed iterations
inline double ANA_STEP           = 0.03; // per-iter step size (fraction of outline dim)
inline double ANA_REPEL          = 8.0;  // repulsion-vs-attraction weight

// In-place FT expansion (post-route): grow undersized soft blocks into the
// adjacent whitespace without repacking, then re-route under a strict gate.
inline int    EXP_ENABLE         = 1;    // 0 = skip the in-place expansion pass
inline int    EXP_ITERS          = 8;    // max grow→reroute rounds

// Minimum routing channel (um) that soft-block growth leaves between a block and
// each NEIGHBOR (the die edge needs none).  DEFAULT 0 (off): a uniform reservation
// globally restricts FT growth and measured WORSE on b50u55 (27 vs 13 penalties) —
// it starves growth everywhere to protect channels that may not be the binding
// ones.  Kept as an env knob (FP_CHMIN) for experiments; the demand-aware channel
// widening (gated) is the targeted replacement.
inline double CHMIN              = 0.0;

// Restarts per worker.  The legacy schedule (sa1 = 0.75 x remaining) gives one
// restart regardless of budget; equal slices give every worker independent
// draws — run-to-run penalty variance (7 vs 18 on b50u60) dwarfs the per-draw
// quality loss from a shorter anneal.
inline int    RESTARTS           = 2;

// Router: per-soft-block feedthrough load cap.  Loads at or below this stay in
// the cheapest FT-rate tier; the router charges escalating cost above it.
inline double FT_SOFT_CAP        = 3000.0;

// ── Deterministic evaluation harness ─────────────────────────────────────────
// Wall-clock-terminated SA makes every run non-reproducible (the iteration count
// depends on machine speed), so single runs can't A/B a post-processing change
// against the ±10-penalty SA noise.  With SA_ITERS>0 the SA instead runs EXACTLY
// that many moves (temperature scheduled by move fraction) and the post-SA
// legalize loop is bounded by iterations only, making the whole pipeline a pure
// function of the seed.  WORKERS>0 forces the worker count (use 1 for a single
// deterministic draw).  ANA_PASSES selects how many escalating analytical-
// legalizer attempts to run (1 = the legacy single mild pass).  Production
// leaves SA_ITERS=0 / WORKERS=0 (time-based, multi-worker).
inline long   SA_ITERS           = 0;
inline int    WORKERS            = 0;
// Report instrumentation: when 1, run_once prints a parseable STAGEREPORT line
// (fails / channel-pen / ft-pen / cost) after each post-processing stage, so a
// benchmark harness can capture the per-stage pre/post penalty breakdown in a
// single run.  Off in production.
inline int    REPORT_STAGES      = 0;
// Analytical-legalizer escalation count.  DEFAULT 1 (the proven single mild
// pass): each extra escalating pass re-routes the layout, and routing dominates
// the post-SA tail, so ANA_PASSES=4 multiplied wall time enough to blow the
// budget.  The escalation is gated (can only cut penalties, never regress) and
// kept behind this knob for harness A/B once a fast valid draw is available.
inline int    ANA_PASSES         = 1;

// ── Pre-SA partitioning stage ────────────────────────────────────────────────
// A congestion-aware recursive min-cut bisection runs BEFORE simulated
// annealing and seeds the initial B*-tree, so the SA starts from a layout where
// (a) strongly-connected (high-net) blocks are already adjacent — relieving long
// high-capacity connections — and (b) connectivity load is spread so no single
// region is a dense knot competing for limited channel/feedthrough capacity.
// The bipartition objective is EXPLICITLY dual-aware: it sums cut weight AND a
// local-congestion-overload term (plus an area-balance term to keep the slicing
// inside the fixed outline).  All terms are carried in "nets" units so the
// weights are directly comparable to the cut.
inline int    PART_ENABLE        = 1;     // 0 = skip partitioning (legacy seed)
inline double PART_CONGW         = 0.6;   // local-congestion-overload weight
inline double PART_BALW          = 0.4;   // area-balance weight (fit the outline)
inline double PART_KCAP          = 25.0;  // channel-capacity proxy: nets per unit
                                          // of region linear dimension (sqrt-area)
inline int    PART_FMPASS        = 8;     // FM refinement passes per bipartition

inline void load_from_env() {
    if (const char* e = std::getenv("FP_HALO"))    HALO = std::atof(e);
    if (const char* e = std::getenv("FP_FTP"))     FT_TRAVERSE_PENALTY = std::atof(e);
    if (const char* e = std::getenv("FP_FTW"))     FTW = std::atof(e);
    if (const char* e = std::getenv("FP_PCONGW"))  PCONGW = std::atof(e);
    if (const char* e = std::getenv("FP_PCTHR"))   PCONG_THRESH = std::atof(e);
    if (const char* e = std::getenv("FP_FTCW"))    FTCW = std::atof(e);
    if (const char* e = std::getenv("FP_FTEC"))    FT_EST_CAP = std::atof(e);
    if (const char* e = std::getenv("FP_MOATW"))   MOATW = std::atof(e);
    if (const char* e = std::getenv("FP_RFBW"))    RFBW = std::atof(e);
    if (const char* e = std::getenv("FP_LEG_ITERS")) LEG_ITERS  = std::atoi(e);
    if (const char* e = std::getenv("FP_LEG_TIME"))  LEG_TIME   = std::atof(e);
    if (const char* e = std::getenv("FP_LEG_TRIES")) LEG_TRIES  = std::atoi(e);
    if (const char* e = std::getenv("FP_LEG_ENABLE")) LEG_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_EXP_ENABLE")) EXP_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_EXP_ITERS"))  EXP_ITERS  = std::atoi(e);
    if (const char* e = std::getenv("FP_CHMIN"))      CHMIN      = std::atof(e);
    if (const char* e = std::getenv("FP_SA_ITERS"))   SA_ITERS   = std::atol(e);
    if (const char* e = std::getenv("FP_WORKERS"))    WORKERS    = std::atoi(e);
    if (const char* e = std::getenv("FP_ANA_PASSES")) ANA_PASSES = std::atoi(e);
    if (const char* e = std::getenv("FP_REPORT_STAGES")) REPORT_STAGES = std::atoi(e);
    if (const char* e = std::getenv("FP_RESTARTS"))   RESTARTS   = std::atoi(e);
    if (const char* e = std::getenv("FP_FT_CAP"))     FT_SOFT_CAP = std::atof(e);
    if (const char* e = std::getenv("FP_PART_ENABLE")) PART_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_PART_CONGW"))  PART_CONGW  = std::atof(e);
    if (const char* e = std::getenv("FP_PART_BALW"))   PART_BALW   = std::atof(e);
    if (const char* e = std::getenv("FP_PART_KCAP"))   PART_KCAP   = std::atof(e);
    if (const char* e = std::getenv("FP_PART_FMPASS")) PART_FMPASS = std::atoi(e);
}

} // namespace cfg
