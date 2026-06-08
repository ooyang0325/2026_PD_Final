#pragma once
#include <cstdlib>

// Global tunable parameters.  Defaults are the production values; each can be
// overridden via an environment variable for offline parameter sweeps.  Values
// are read once at program start (cfg::load_from_env) and treated as read-only
// thereafter, so they are safe to share across worker threads.
namespace cfg {

inline double HALO               = 10.0; // routing-channel gap around every block
inline double FT_TRAVERSE_PENALTY = 3.0; // router cost multiplier for feedthrough
inline double FTW                = 0.4; // SA feedthrough-minimizing penalty weight
inline double PCONGW             = 1.5; // SA placement-congestion weight (routing-aware)
inline double PCONG_THRESH       = 500; // only high-demand connections drive congestion

// Route-driven legalize loop (post-SA refinement).
inline int    LEG_ITERS          = 30;   // max outer iterations
inline double LEG_TIME           = 30.0; // wall-clock cap (seconds)
inline int    LEG_TRIES          = 8;    // displace mini-SA tries per attempt
inline int    LEG_ENABLE         = 1;    // 0 = skip the loop (use legacy ft_iter)

// Final-pass analytical (force-directed) legalizer.  Breaks the B*-tree and
// outputs continuous coordinates.  Strict rollback if it doesn't improve.
inline int    ANA_ENABLE         = 1;    // 0 = skip the analytical pass
inline int    ANA_ITERS          = 25;   // force-directed iterations
inline double ANA_STEP           = 0.03; // per-iter step size (fraction of outline dim)
inline double ANA_REPEL          = 8.0;  // repulsion-vs-attraction weight

inline void load_from_env() {
    if (const char* e = std::getenv("FP_HALO"))    HALO = std::atof(e);
    if (const char* e = std::getenv("FP_FTP"))     FT_TRAVERSE_PENALTY = std::atof(e);
    if (const char* e = std::getenv("FP_FTW"))     FTW = std::atof(e);
    if (const char* e = std::getenv("FP_PCONGW"))  PCONGW = std::atof(e);
    if (const char* e = std::getenv("FP_PCTHR"))   PCONG_THRESH = std::atof(e);
    if (const char* e = std::getenv("FP_LEG_ITERS")) LEG_ITERS  = std::atoi(e);
    if (const char* e = std::getenv("FP_LEG_TIME"))  LEG_TIME   = std::atof(e);
    if (const char* e = std::getenv("FP_LEG_TRIES")) LEG_TRIES  = std::atoi(e);
    if (const char* e = std::getenv("FP_LEG_ENABLE")) LEG_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_ANA_ENABLE")) ANA_ENABLE = std::atoi(e);
    if (const char* e = std::getenv("FP_ANA_ITERS"))  ANA_ITERS  = std::atoi(e);
    if (const char* e = std::getenv("FP_ANA_STEP"))   ANA_STEP   = std::atof(e);
    if (const char* e = std::getenv("FP_ANA_REPEL"))  ANA_REPEL  = std::atof(e);
}

} // namespace cfg
