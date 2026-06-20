#pragma once
#include <cstdlib>
#include <string>

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

// Route-driven legalize loop (post-SA refinement).
inline int    LEG_ITERS          = 128;   // max outer iterations
inline double LEG_TIME           = 30.0; // wall-clock cap (seconds)
inline int    LEG_TRIES          = 8;    // displace mini-SA tries per attempt
inline int    LEG_ENABLE         = 1;    // 0 = skip the loop (use legacy ft_iter)

// Final-pass analytical (force-directed) legalizer.  Breaks the B*-tree and
// outputs continuous coordinates.  Strict rollback if it doesn't improve.
inline int    ANA_ENABLE         = 1;    // 0 = skip the analytical pass
inline int    ANA_ITERS          = 128;  // force-directed iterations
inline double ANA_STEP           = 0.03; // per-iter step size (fraction of outline dim)

// Mathematical-programming (analytical / ePlace-style) engine.
inline std::string ENGINE        = "mp"; // "sa" | "mp" (env: FP_ENGINE)
inline int    MP_GRID            = 128;  // density-bin count per axis
inline double MP_TARGET_OVF      = 0.08; // overflow target (Eq. 37) for termination
inline int    MP_MAX_ITER        = 1500; // Nesterov iteration cap
inline double MP_INIT_LAMBDA     = 8e-5; // λ0 scale (ePlace Eq. 35)
inline double MP_PHI_MIN         = 0.95; // phiCoef clamp lower
inline double MP_PHI_MAX         = 1.05; // phiCoef clamp upper
inline double MP_HALO_SCALE      = 1.0;  // demand-driven halo multiplier
inline int    MP_RB_ITERS        = 30;   // routability outer-loop iterations
                                          // (drives route → SEP-bump → re-legalize
                                          //  feedback; 30 lets the boost matrix
                                          //  saturate for outliers AND newly-
                                          //  discovered overflows after earlier
                                          //  ones clear)
inline double MP_TDENSITY        = 0.0;  // ρt target (0 = auto from utilization)
inline int    MP_RB_DEBUG        = 0;    // 1 = trace routability outer loop to stderr

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
    if (const char* e = std::getenv("FP_ENGINE"))      ENGINE = e;
    if (const char* e = std::getenv("MP_GRID"))        MP_GRID = std::atoi(e);
    if (const char* e = std::getenv("MP_TARGET_OVF"))  MP_TARGET_OVF = std::atof(e);
    if (const char* e = std::getenv("MP_MAX_ITER"))    MP_MAX_ITER = std::atoi(e);
    if (const char* e = std::getenv("MP_INIT_LAMBDA")) MP_INIT_LAMBDA = std::atof(e);
    if (const char* e = std::getenv("MP_PHI_MIN"))     MP_PHI_MIN = std::atof(e);
    if (const char* e = std::getenv("MP_PHI_MAX"))     MP_PHI_MAX = std::atof(e);
    if (const char* e = std::getenv("MP_HALO_SCALE"))  MP_HALO_SCALE = std::atof(e);
    if (const char* e = std::getenv("MP_RB_ITERS"))    MP_RB_ITERS = std::atoi(e);
    if (const char* e = std::getenv("MP_TDENSITY"))    MP_TDENSITY = std::atof(e);
    if (const char* e = std::getenv("MP_RB_DEBUG"))    MP_RB_DEBUG = std::atoi(e);
}

} // namespace cfg
