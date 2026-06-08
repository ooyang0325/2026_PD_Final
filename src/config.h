#pragma once
#include <cstdlib>

// Global tunable parameters.  Defaults are the production values; each can be
// overridden via an environment variable for offline parameter sweeps.  Values
// are read once at program start (cfg::load_from_env) and treated as read-only
// thereafter, so they are safe to share across worker threads.
namespace cfg {

inline double HALO               = 2.0; // routing-channel gap around every block
inline double FT_TRAVERSE_PENALTY = 3.0; // router cost multiplier for feedthrough
inline double FTW                = 0.4; // SA feedthrough-minimizing penalty weight
inline double PCONGW             = 1.5; // SA placement-congestion weight (routing-aware)
inline double PCONG_THRESH       = 500; // only high-demand connections drive congestion

inline void load_from_env() {
    if (const char* e = std::getenv("FP_HALO"))    HALO = std::atof(e);
    if (const char* e = std::getenv("FP_FTP"))     FT_TRAVERSE_PENALTY = std::atof(e);
    if (const char* e = std::getenv("FP_FTW"))     FTW = std::atof(e);
    if (const char* e = std::getenv("FP_PCONGW"))  PCONGW = std::atof(e);
    if (const char* e = std::getenv("FP_PCTHR"))   PCONG_THRESH = std::atof(e);
}

} // namespace cfg
