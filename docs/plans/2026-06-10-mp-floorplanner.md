# Plan: Mathematical-Programming Floorplanner (ePlace-style analytical NLP)

**Branch:** `feature/mp-floorplanner` (created, branched from `feature/floorplanner`)
**Date:** 2026-06-10
**Execution:** phase-by-phase via `/claude-mem:do`; each phase is self-contained and runnable in a fresh context.

## Goal and priority hierarchy (USER DIRECTIVE — overrides everything)

Replace the simulated-annealing global engine with an analytical (nonlinear-programming, ePlace-style) floorplanner. The known weakness of the SA flow is that it **cannot generate suitable channels for nets to route** — channel overflow and soft-block area (feedthrough) overflow penalties dominate.

Optimization priority, strictly lexicographic:

1. **0 FAILs** — outline, no overlap, in-bounds, EDGE locations, routing completeness
2. **0 penalties** — channel overflow (`nets > dimension × 25`) and soft-block undersize (`area < FT target`)
3. **Minimize cost** — `chip_w × chip_h + α·HPWL`

Wirelength is NOT the first target. **Never accept a change that trades penalty count for cost.** Every acceptance test in this plan compares `(fails, penalties, cost)` lexicographically.

## Locked decisions (user-confirmed)

| Decision | Choice |
|---|---|
| MP method | Analytical NLP: WA wirelength + electrostatic density, Nesterov solver (ePlace) |
| Dependencies | Vendoring open-source OK; plan vendors only OpenROAD's FFT/Poisson files (BSD-3) |
| Scope | Replace `SAOptimizer` only, behind `FP_ENGINE=mp\|sa`; keep LegalizeLoop / AnalyticalLegalizer / GlobalRouter / scoring / `run_search` harness |

---

## Phase 0: Documentation Discovery — COMPLETE (consolidated findings)

Two discovery passes were run (codebase map + math/reference research). All formulas below were verified against primary sources; reference copies live **in this repo** under `docs/refs/`.

### 0.1 Reference materials (local copies)

| Path | Content | License |
|---|---|---|
| `docs/refs/openroad-gpl/fft.h`, `fft.cpp` | Complete spectral Poisson solve (`doFFT()`: DCT → coefficients → ψ, ξx, ξy) | BSD-3 |
| `docs/refs/openroad-gpl/fftsg.cpp`, `fftsg2d.cpp` | Ooura FFT/DCT/DST package (`ddct2d`, `ddsct2d`, `ddcst2d`) — plain C, zero deps | BSD-3 |
| `docs/refs/openroad-gpl/nesterovBase.cpp/.h` | WA forces/gradient (lines 1293–1394, 1460–1524), preconditioner (1526–1530, 2572–2578), steplength (2725–2743, 3107–3167), λ init/update (2701–2723, 2975–2986, 3098), density gradient (2582–2601) | BSD-3 |
| `docs/refs/openroad-gpl/nesterovPlace.cpp/.h` | Nesterov loop + backtracking (871–928, 1064–1138), γ schedule (1172–1186), divergence guards (560–640) | BSD-3 |
| `docs/refs/openroad-gpl/Replace.h` | All default parameter values (lines 62–83) | BSD-3 |
| `docs/refs/dreamplace/wa_functional.h` | WA gradient kernel, cleanest formula (lines ~82–187) | BSD-3 |
| `docs/refs/dreamplace/electric_potential.py`, `electric_force.cpp`, `discrete_spectral_transform.py` | DCT pipeline cross-check | BSD-3 |
| `docs/refs/papers/eplace-todaes2015.pdf` | ePlace: Eq. 6 (WA), 14–16 (density objective), 20–24 (Poisson/spectral), 25 (local smoothing), 29 (Lipschitz), 30–33 (preconditioner), 35–36 (λ), 37 (overflow), 38 (γ), Alg. 2 (Nesterov) | — |
| `docs/refs/papers/pef-tcad2023.pdf` | PeF fixed-outline floorplanning: Eq. 12 (model w/ width variables), 16/18 (∂/∂w), Alg. 1–2 (alternating Nesterov + projection), §IV-B (constraint-graph legalization) | — |
| `docs/refs/papers/iccad2023-rotation-fp.pdf` | Eq. 5: wirelength sizing force ∂LSE/∂w via pin-offset chain rule; two-step schedule | — |

### 0.2 Existing-codebase interfaces (the contract the new engine must honor)

| Item | Location |
|---|---|
| Pipeline | `src/main.cpp:83–373` (`run_once`), `384–421` (`run_search`, multi-threaded multi-start) |
| Data model | `src/types.h` — `Block{lx,ly,width,height,area,min_ar,max_ar,type,locations,ft}`, `Connection{from,to,nets}`, `Channel{cap_x=height×25, cap_y=width×25}`, `Design` |
| FT target area | `src/types.h:29–36` — `get_target_area(ft_nets) = (√area + (ft_nets/25)·rate/2)²` |
| SA engine being replaced | `src/sa_optimizer.h` (ctor `(Floorplan&, unsigned seed)`, `run()`) |
| Scoring | `src/main.cpp:48–78` (`compute_final_cost`), `202–221` (`route_and_score`); ground truth `script/evaluator.py` |
| FT collection after routing | `src/main.cpp:102–114` (`collect_ft`) |
| Edge snapping | `src/floorplan.h:125–149` (overlap-safe boundary snap) |
| Channels from coordinates | `src/channel.h` (vertical-strip decomposition — works on raw coords, reusable as-is) |
| Router | `src/router.h:179` `route_all(Design&, int max_rr)` — face-graph Dijkstra + rip-up/reroute |
| Downstream legalizers | `src/legalize_loop.h` (tree-based DISPLACE — needs a no-tree mode), `src/analytical_legalizer.h` (coordinate-based — compatible as-is) |
| Config pattern | `src/config.h` + `load_from_env()` (`main.cpp:431`) |
| Build | `Makefile`: g++ `-std=c++17 -O3 -march=native -Wall -Wextra`, single TU (`src/main.cpp`), header-only style |

### 0.3 Allowed APIs / formulas (cite-only — do NOT invent variants)

- **WA wirelength + gradient**: copy `docs/refs/dreamplace/wa_functional.h` or `nesterovBase.cpp:1460–1524`. Mandatory: bounding-box exponent shift (`exp((x_i−x_max)/γ)`), exponent cutoff −300.
- **γ schedule**: ePlace Eq. 38 / `nesterovPlace.cpp:1172–1186`: overflow>1.0 → 1/γ=0.1·base; <0.1 → 10·base; else `1/pow(10,(τ−0.1)·20/9−1)·base`, `base = 0.25/avg_bin_size`.
- **Poisson solve**: vendor `fft.cpp/fft.h/fftsg.cpp/fftsg2d.cpp` **wholesale**. `doFFT()` already does: `ddct2d` forward → scale → `ψ=ρ/(wx²+wy²)`, `ξx=ψ·wx`, `ξy=ψ·wy`, zero DC → `ddct2d(ψ)`, `ddsct2d(ξx)`, `ddcst2d(ξy)`, including non-square-bin frequency scaling (fft.cpp:54–57).
- **Density gradient per block**: Σ over overlapped bins of `overlapArea × field` (`nesterovBase.cpp:2582–2601`); local smoothing for sub-bin dims (ePlace Eq. 25); **fixed objects' charge scaled by ρt** (ePlace §3.2).
- **Nesterov**: ePlace Alg. 2 as implemented in `nesterovPlace.cpp:1064–1138` + `doBackTracking` (871–928): `a_{k+1}=(1+√(4a_k²+1))/2`, steplength `=‖Δcoord‖/‖Δgrad‖`, accept backtrack when `newStep > 0.95·step`, max 10 rounds, floor 0.01, NaN/Inf → divergence revert.
- **Preconditioner**: `h_i = max(1, #pins_i + λ·q_i)` (ePlace Eq. 30–33; `nesterovBase.cpp:1526–1530, 2572–2578, 2788–2801`).
- **λ init/update**: Eq. 35 (`λ0 = Σ|∇W| / Σ q|ξ| × 8e-5`); update `λ *= phiCoef`, `phiCoef = maxPhi·pow(maxPhi, −ΔHPWL/refHPWL)` clamped to `[0.95, 1.05]` (`nesterovBase.cpp:2975–2986`).
- **Overflow**: τ = Σ_bins max(ρ'_b − ρt, 0)·A_b / Σ movable area (Eq. 37), ρ' counts movable area only.
- **Soft-block sizing forces** (stretch only): PeF Eq. 18 (energy ∂/∂w, midpoint approximation) and ICCAD'23 Eq. 5 (`∂W/∂w_i = ∂W/∂x_pin·0.5 − ∂W/∂y_pin·0.5·A_i/w_i²` for center pins), with AR projection — PeF Alg. 1 alternation.

### 0.4 Anti-patterns (documented pitfalls — guard in every phase)

1. **Sign convention**: OpenROAD stores *descent* directions and updates `x += α·g`. Pick this convention, document it in one header comment, unit-test on a 2-pin net before anything else.
2. **DCT/DST pairing**: ξx uses `ddsct2d` (sin-x/cos-y), ξy uses `ddcst2d`. Swapping produces plausible but rotationally wrong forces.
3. **Never mix normalization conventions**: paper uses `2πu/m`, Ooura/OpenROAD `πu/m` with its own scalings. Copy `fft.cpp` wholesale; do not re-derive coefficients.
4. **γ too small early** → vanishing gradients. Use the overflow-driven schedule; never hardcode a tiny γ.
5. **Whitespace handling**: vanilla ePlace needs fillers OR overflow measured against ρt < 1. This plan uses **ρt = inflated utilization** (no fillers) — and that is a *feature* here: whitespace IS channel capacity. Add fillers only if blocks clump/oscillate (Phase 7 fallback).
6. **Fixed-object density scaling**: scale EDGE (fixed) block charge by ρt, or they over-repel and create dead halos.
7. **Divergence guards**: NaN/Inf steplength → revert to best snapshot, retry with smaller λ0; keep a best-`(fails,penalties,cost)` snapshot at all times.
8. **Grid too coarse** → oscillation. 64–256 bins is the right range for ≤60 blocks; make it a knob.
9. **Invented APIs**: do not add parameters or methods not present in the cited references or this repo. If a formula isn't in §0.3, it doesn't exist.

---

## Architecture (synthesis)

```
run_once (FP_ENGINE=mp path)                        [main.cpp branch]
  Floorplan fp(d)            — reused only for W/H arrays + helpers
  MPOptimizer mp(fp, d, seed)
    ├─ seeded random init inside FULL max outline   (multi-start via run_search workers)
    ├─ ROUTABILITY OUTER LOOP (the core, Phase 5):
    │    1. Nesterov global place (x,y) w/ demand-inflated block charges
    │    2. CG-legalize lightly → channels → GlobalRouter.route_all (cheap at this scale)
    │    3. measure: per-channel overflow, per-soft-block ft_nets
    │    4. resize soft blocks to get_target_area(ft_nets); bump halos near
    │       overflowed channels (×1.2); warm-restart Nesterov
    │    5. stop at 0 penalties or MP_RB_ITERS; keep lexicographic-best snapshot
    └─ writes d.blocks[i].lx/ly + fp.W/H
  CGLegalizer (overlap removal with halo spacing; NO aggressive compaction)
  edge snap → collect_ft → [AnalyticalLegalizer] → [LegalizeLoop, no-tree mode]
  optional final compaction: accept ONLY if (fails,penalties,cost) improves
  route_and_score → best
```

**Routability-first design decisions (from the user directive):**

- **Placement region = the full max outline.** Cost (area) is tertiary; whitespace is channel capacity. Shrinking happens only in the final strictly-validated compaction step.
- **Demand-driven halos**: block i with total incident nets `T_i = Σ_j nets_ij` needs `T_i/25` µm of channel cross-section at its boundary. Per-side halo: `h_i = max(cfg::HALO/2, MP_HALO_SCALE · T_i/(2·k_i·25))`, `k_i` = usable sides (4 interior, 3 EDGE). Adjacent gap ≈ `h_i + h_j`. Density charges and legalizer spacing use inflated dims `w_i + 2h_i`.
- **FT-aware sizing uses the REAL router in the loop** (it costs ~ms at ≤60 blocks), not just the FTAFP estimate — soft undersize penalties are killed by construction because areas are set to the measured FT target before the final pass.
- **Target density ρt = Σ inflated block area / outline area** (clamped to ≤ 0.9): equilibrium = uniform spread with demand-proportional whitespace everywhere.

**New files** (matching flat header-only style): `src/mp_optimizer.h` (façade), `src/mp/wa_wirelength.h`, `src/mp/density.h`, `src/mp/nesterov.h`, `src/mp/cg_legalizer.h`, `src/vendor/{fft.h,fft.cpp,fftsg.cpp,fftsg2d.cpp}`, `tests/test_mp.cpp`.

**New config knobs** (in `src/config.h`, env-driven like existing ones): `FP_ENGINE` ("sa" default), `MP_GRID` (128), `MP_TARGET_OVF` (0.08), `MP_MAX_ITER` (1500), `MP_INIT_LAMBDA` (8e-5), `MP_PHI_MIN/MAX` (0.95/1.05), `MP_HALO_SCALE` (1.0), `MP_RB_ITERS` (4), `MP_TDENSITY` (0 = auto from utilization).

---

## Phase 1: Scaffold, vendoring, engine switch

**Implement:**
1. Copy `docs/refs/openroad-gpl/{fft.h,fft.cpp,fftsg.cpp,fftsg2d.cpp}` → `src/vendor/`, keeping BSD-3 headers. Strip OpenROAD-only includes (`utl/Logger` etc.) from `fft.cpp/h` — replace logging with nothing; verify the four files compile standalone.
2. `Makefile`: add `src/vendor/fft.cpp src/vendor/fftsg.cpp src/vendor/fftsg2d.cpp` as translation units (this repo gains its first multi-TU build — keep flags identical).
3. `src/config.h`: add the knobs table above, following the exact `load_from_env` pattern at `config.h:29–39`.
4. `src/mp_optimizer.h`: stub `MPOptimizer{ MPOptimizer(Floorplan&, Design&, unsigned seed); void run(); }` — for now seeded random scatter of blocks inside the outline writing `d.blocks[i].lx/ly`.
5. `src/main.cpp` `run_once`: branch on `cfg::ENGINE` — `"mp"` skips the two SA passes and tree-finalize, calls the stub, then proceeds to `route_and_score`. **Do not touch the SA path.**

**Docs to read first:** `Makefile`, `src/config.h` (whole file), `src/main.cpp:83–235`, vendored sources after copying.

**Verification:**
- [ ] `make` clean with `-Wall -Wextra` (no new warnings)
- [ ] `FP_ENGINE` unset → bit-identical behavior: run one case with a fixed seed before/after the change, same final score
- [ ] `FP_ENGINE=mp ./bin/solver <case>` runs to completion and writes a parseable `.cfg` (FAILs acceptable at this phase)
- [ ] `grep -rn "utl::" src/vendor/` → empty

**Anti-pattern guards:** keep license headers; no edits inside SA/legalizer code; guard #9.

## Phase 2: Density engine + spectral Poisson

**Implement** `src/mp/density.h`: bin grid (`MP_GRID`² over the outline), per-bin movable/fixed density accumulation with area-overlap, local smoothing (ePlace Eq. 25) for sub-bin dims, EDGE-block charge scaled by ρt, overflow τ (Eq. 37 vs ρt), Poisson solve via vendored `doFFT` pattern, per-block density descent direction by bin-overlap accumulation (`nesterovBase.cpp:2582–2601`). Density uses **inflated dims** (`w+2h_i`) — implement the halo computation here (formula in Architecture) with `MP_HALO_SCALE`.

**Implement** `tests/test_mp.cpp` + `make test` target (plain asserts, no framework — repo has none).

**Docs to read first:** `docs/refs/openroad-gpl/fft.cpp` (all 185 lines), `nesterovBase.cpp:2582–2601`, ePlace pp. 14–17 (Eq. 20–25, 37), §0.4 guards 2/3/5/6.

**Verification:**
- [ ] FD test: ∂N/∂x_i matches central finite difference within 1e-2 relative on a 5-block toy (this validates DCT/DST pairing — guard #2)
- [ ] Two overlapping blocks → equal-and-opposite separating forces
- [ ] Σψ ≈ 0 (DC removed); uniform density → ‖ξ‖ ≈ 0
- [ ] τ = 0 when blocks spread uniformly at ρt

## Phase 3: WA wirelength

**Implement** `src/mp/wa_wirelength.h`: nets from `Design.connections` (2-pin, pin = block center `lx+w/2`, net weight = `nets`), WA cost + per-block gradient copied from `docs/refs/dreamplace/wa_functional.h` (exp shift + cutoff), γ from the overflow schedule. Store **descent** directions (guard #1) — one comment block documenting the convention.

**Verification:**
- [ ] FD gradient check within 1e-3 relative (random 10-block instance, several γ)
- [ ] 2-pin net: gradient pulls pins together; with γ → 0.01·span, WA → exact weighted HPWL within 1%
- [ ] Weighted: doubling `nets` doubles gradient

## Phase 4: Nesterov global placement (x, y)

**Implement** `src/mp/nesterov.h` + wire into `MPOptimizer::run`: copy the OpenROAD loop skeleton (a_k, coeff, backtracking, steplength, accept rule, floor, divergence revert + λ0-shrink retry, synthetic prev-point init with ×10 perturbation retry), Jacobi preconditioner, λ0 init (Eq. 35 × `MP_INIT_LAMBDA`), λ·phiCoef update, coordinate clamping to outline. EDGE blocks fixed at nearest allowed `locations` entry (flips come later); HARD macros movable with fixed dims. Init: seeded scatter (different per `run_search` worker). Terminate at `MP_TARGET_OVF`, `MP_MAX_ITER`, or stall.

**Docs to read first:** `nesterovPlace.cpp:871–928, 1064–1186`, `nesterovBase.cpp:2701–2743, 2975–3167`, ePlace Alg. 2 + Eq. 29–36.

**Verification:**
- [ ] On `50_block` and `50b60u` inputs: τ descends to ≤ 0.10, no NaN, < 2000 iters, < ~2 s per start
- [ ] Log per-iter `(iter, τ, HPWL, λ, step)` behind a debug env flag; trajectory monotone-ish in τ
- [ ] Different seeds → different layouts (multi-start works), all reaching target τ

## Phase 5: Routability outer loop (CORE — attacks the user's stated pain point)

**Implement** in `MPOptimizer::run` the outer loop from the Architecture diagram: light CG-legalize (Phase 6 header, can be built first in skeletal form) → `channel.h` channels → `GlobalRouter::route_all` → measure per-channel overflow + per-block `ft_nets` (copy `main.cpp:102–114`) → resize soft blocks to `get_target_area(ft_nets)` preserving current AR (clamp to `[min_ar, max_ar]`) → bump halos of blocks adjacent to overflowed channels ×1.2 → warm-restart Nesterov (positions kept, λ re-initialized ×0.1). Keep the lexicographic-best `(fails, penalties, cost)` snapshot across iterations; return it.

**Stretch (only after the loop works):** in-loop continuous width variable per PeF Alg. 1 / ICCAD'23 Eq. 5 with AR projection.

**Verification:**
- [ ] On the worst-penalty case (b50u60 class, SA baseline ≈ 21+ penalties): penalties strictly decrease across outer iterations
- [ ] Soft-undersize penalties = 0 after the loop (areas set to measured targets by construction)
- [ ] Snapshot logic: forced-bad final iteration still returns the earlier best

**Anti-pattern guards:** never shrink a halo below its demand bound; never accept penalty regression for cost (user directive); guard #7 snapshots.

## Phase 6: Legalization bridge + full pipeline integration

**Implement:**
1. `src/mp/cg_legalizer.h`: constraint-graph legalization per PeF §IV-B — build HCG/VCG from relative positions, longest-path placement using **halo-inflated dims** as spacing, push minimally to clear overlaps inside the outline; do NOT compact whitespace away. Infeasible → report, caller falls back to best snapshot.
2. Coordinate-based EDGE snap (pattern: `floorplan.h:125–149`).
3. `run_once` MP path completed: MP → CG-legalize → edge snap → `collect_ft` → AnalyticalLegalizer (already coordinate-based) → LegalizeLoop with a new `tree_moves=false` flag (DISPLACE skipped or replaced by small coordinate nudges; ROTATE/EXPAND unchanged) → `route_and_score`.
4. Final optional compaction: longest-path compaction toward one corner with halos preserved, accepted only on lexicographic improvement verified by `route_and_score`.

**Verification:**
- [ ] `python3 script/evaluator.py` on ALL inputs in `input/`: **0 FAILs** with `FP_ENGINE=mp`
- [ ] Output `.cfg` matches `output.h` format (OUTLINE/BLOCK/CHANNEL/PATH sections)
- [ ] `FP_ENGINE=sa` still bit-identical to Phase-1 baseline
- [ ] Compaction step: zero cases where penalties increased

## Phase 7: Benchmark & tune

Run `script/bench.sh` for both engines across all cases (same wall-clock budget). Produce an A/B table: fails / penalties / cost / runtime per case, vs the SA baselines in `bench_baseline_v5.txt` and `bench_halo30_*.txt`. Tune: `MP_GRID` ∈ {64,128}, `MP_TARGET_OVF` ∈ {0.05,0.08,0.12}, `MP_HALO_SCALE` ∈ {0.5,1.0,1.5}, ρt clamp. Fallback if blocks clump/oscillate: add ePlace fillers (`nesterovBase.cpp:2060–2228` pattern).

**Acceptance gates (ordered, lexicographic):**
- [ ] **G1:** 0 FAILs on every case (hard requirement)
- [ ] **G2:** total penalties ≤ SA baseline on every case; target 0 on the key cases (b50u60)
- [ ] **G3:** cost reported side-by-side (win not required for plan completion — document the gap)

## Phase 8: Final verification (skill-mandated)

- [ ] `make test` green (all FD checks from Phases 2–3 still pass)
- [ ] Full evaluator suite re-run, results recorded in a `docs/` bench note
- [ ] Anti-pattern greps: `grep -rn "2 \* M_PI\|2\*M_PI" src/mp src/vendor` (convention leakage — must be empty in our code), `grep -rn "rand()" src/mp` (must use seeded `mt19937`), `grep -n ddsct2d src/mp src/vendor/fft.cpp` confirms ξx/ξy pairing untouched
- [ ] Implementation-vs-citation audit: every formula in `src/mp/*` carries a comment naming its source (ePlace Eq. N / file:line from §0.3); no uncited math
- [ ] `FP_ENGINE=sa` regression: still identical to pre-branch behavior
- [ ] Update `README.md` with the engine switch and knobs; commit history follows `feat:`/`docs:` convention
