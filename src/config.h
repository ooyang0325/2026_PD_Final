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

// Router: per-soft-block feedthrough load cap.  Loads at or below this stay in
// the cheapest FT-rate tier; the router charges escalating cost above it.
inline double FT_SOFT_CAP        = 3000.0;

inline void load_from_env() {
    if (const char* e = std::getenv("FP_HALO"))    HALO = std::atof(e);
    if (const char* e = std::getenv("FP_FTP"))     FT_TRAVERSE_PENALTY = std::atof(e);
    if (const char* e = std::getenv("FP_FTW"))     FTW = std::atof(e);
    if (const char* e = std::getenv("FP_PCONGW"))  PCONGW = std::atof(e);
    if (const char* e = std::getenv("FP_PCTHR"))   PCONG_THRESH = std::atof(e);
    if (const char* e = std::getenv("FP_FTCW"))    FTCW = std::atof(e);
    if (const char* e = std::getenv("FP_FTEC"))    FT_EST_CAP = std::atof(e);
    if (const char* e = std::getenv("FP_MOATW"))   MOATW = std::atof(e);
    if (const char* e = std::getenv("FP_LEG_ITERS")) LEG_ITERS  = std::atoi(e);
    if (const char* e = std::getenv("FP_LEG_TIME"))  LEG_TIME   = std::atof(e);
    if (const char* e = std::getenv("FP_LEG_TRIES")) LEG_TRIES  = std::atoi(e);
    if (const char* e = std::getenv("FP_LEG_ENABLE")) LEG_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_EXP_ENABLE")) EXP_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_EXP_ITERS"))  EXP_ITERS  = std::atoi(e);
    if (const char* e = std::getenv("FP_FT_CAP"))     FT_SOFT_CAP = std::atof(e);
}

} // namespace cfg
