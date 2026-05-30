# ICCAD 2026 Problem E Optimization Tracking

This document tracks the routing/fail-reduction optimization roadmap for the current codebase.

## Goal

Primary objective:
- Reduce `FAIL` count to 0 (especially `Routing open` and `Edge constraint violation`).

Secondary objective:
- Reduce penalties (`channel overflow`, `feedthrough overflow`) and total cost under fixed runtime budget.

## Baseline (Current Status)

Based on recent evaluator output:
- `Total Fails: 13`
- `Total Penalties: 2`
- Key fail types:
  - `Edge constraint violation`
  - `Routing open` (multiple pairs with routed nets = 0)

Known root-cause summary:
- Routing graph becomes disconnected (multiple traversable components), causing unreachable source/sink pairs.
- Edge constraints are treated mainly by penalty and can still end up violated.

---

## Metrics and Evaluation Protocol

Use a fixed regression set (including current failing testcase) and track the following per run:

- `Fail count` (highest priority)
- `Routing open` pair count
- `Edge constraint` fail count
- `Penalty count` (channel overflow + feedthrough overflow)
- `Total cost`
- Runtime (same time limit setting for fair comparison)

Run command template:

```bash
./bin/solver <input.csv> <output.cfg> [time_limit_sec]
python3 script/evaluator.py <input.csv> <output.cfg>
```

Acceptance priority:
1. Fail count down toward 0
2. Penalty reduction
3. Cost reduction
4. Runtime stability

---

## Phase Plan

## Phase 1 - Connectivity-Driven SA (Highest Priority)

Objective:
- Prevent disconnected routing states during floorplanning, so `Routing open` fails are avoided early.

Planned changes:
- Add fast connectivity proxy (`P_disc`) in SA cost:
  - Build channels from current placement.
  - Build simplified traversable adjacency graph.
  - Compute connected components.
  - Penalize high-net pairs that are unreachable.
- Use top-k critical nets only for speed (`k = 20~40`, start from 30).
- Add phase-based weight schedule:
  - Early SA: low `beta`
  - Late SA: higher `beta`

Target formula:

```text
C = Area + alpha*HPWL + beta(t)*P_disc + gamma*P_edge + delta*P_outline
```

Expected impact:
- Significant drop in `Routing open` fails.

Status checklist:
- [ ] Define data flow for connectivity proxy in floorplan/SA path.
- [ ] Implement traversable graph builder for current placement.
- [ ] Implement top-k critical-net unreachable penalty.
- [ ] Add `beta(t)` schedule in SA.
- [ ] Add debug logs for connectivity-related terms.
- [ ] Run regression and record before/after.

---

## Phase 2 - Edge Constraint Hardening

Objective:
- Reduce/eradicate edge placement violations.

Planned changes:
- Add dedicated edge-aware moves (snap/slide toward legal boundary options).
- Reject obviously illegal edge states earlier (semi-hard constraint).
- Add post-SA edge repair pass before final routing.

Expected impact:
- Edge-related fail count trends to 0.

Status checklist:
- [ ] Implement edge-specific move operator.
- [ ] Add legality gate/reject policy for edge placement.
- [ ] Add lightweight edge-repair pass.
- [ ] Validate on edge-heavy testcases.

---

## Phase 3 - Early Congestion Estimation (RUDY-like)

Objective:
- Estimate and suppress congestion before full global routing.

Planned changes:
- Add coarse grid bins (`16x16` or `20x20`).
- Demand estimation from net bounding boxes (RUDY-like).
- Capacity approximation from available free space / channels.
- Add overflow-based congestion penalty term:

```text
P_cong = sum_bins( max(0, dx/cx - 1)^q + max(0, dy/cy - 1)^q ), q in [1.5, 2]
```

Updated objective:

```text
C = Area + alpha*HPWL + beta*P_disc + eta*P_cong + gamma*P_edge
```

Expected impact:
- Lower channel overflow and better downstream routability.

Status checklist:
- [ ] Implement coarse grid model.
- [ ] Implement RUDY-like demand accumulation.
- [ ] Implement capacity approximation.
- [ ] Add `eta(t)` schedule and parameter tuning.
- [ ] Validate fail/penalty/cost/runtime tradeoff.

---

## Phase 4 - Router-Side Robustness

Objective:
- Reduce empty-path outcomes and improve critical-net success rate.

Planned changes:
- Add routeability precheck for each connection pair.
- Annotate disconnected blockers for debugging/feedback.
- Enhance critical-net ordering and reroute strategy for hard nets.

Expected impact:
- Fewer empty paths in final output.

Status checklist:
- [ ] Add pre-route reachability diagnostics.
- [ ] Improve route ordering heuristics for critical nets.
- [ ] Add focused reroute retries for failed pairs.
- [ ] Verify evaluator `Routing open` reduction.

---

## Phase 5 (Optional) - Analytical Seeding (Hybrid)

Objective:
- Improve initial placements for SA via lightweight analytical initialization.

Rationale:
- Problem size is small (`n <= 50`), analytical seed generation is practical and cheap.
- Better seeds can reduce SA time to feasibility.

Planned approach:
- Generate 3-5 analytical initial placements (quadratic/force-directed style with boundary soft constraints).
- Run existing SA+routing pipeline per seed.
- Keep best final solution.

Expected impact:
- Faster convergence to feasible solutions; potentially better final cost.

Status checklist:
- [ ] Prototype analytical initializer.
- [ ] Integrate multi-seed launch path.
- [ ] Compare against pure-random/multi-start baseline.

---

## Suggested Initial Hyperparameters

- Critical net count: `k = 30`
- Connectivity weight schedule: `beta: 0 -> 5e5`
- Congestion weight schedule: `eta: 0 -> 5e3` (enabled later in SA)
- Congestion exponent: `q = 1.8`
- Edge legality tolerance: `1e-3`

Note:
- These are initial guesses and should be tuned with regression metrics, not fixed assumptions.

---

## Implementation Order (Recommended)

1. Phase 1 (connectivity proxy)  
2. Phase 2 (edge hardening)  
3. Phase 3 (RUDY-like congestion)  
4. Phase 4 (router robustness)  
5. Phase 5 optional (analytical seed)

Reason:
- Current largest blocker is fail feasibility (`Routing open`), so routeability-first optimization is highest ROI.

---

## Progress Log

Use this section to append dated updates while implementing.

### 2026-05-28
- Created optimization tracking document.
- Confirmed primary blocker: routing-graph disconnection causing `Routing open`.
- Agreed to prioritize fail elimination before cost refinement.
- Ran short baseline benchmark on 5 generated small cases (`n=12,14,16,18,20`, `time_limit=20s`):
  - Fails: `1,2,8,5,13`
  - Penalties: `0,0,1,3,0`
- Implemented first trial of connectivity-aware SA cost:
  - Added `connectivity_penalty()` in `floorplan.h`.
  - Added time-ramped connectivity weight in `sa_optimizer.h`.
- Re-ran same benchmark after change:
  - Fails: `1,4,12,10,12`
  - Penalties: `0,0,0,3,1`
- Result summary:
  - Improvement was not consistent; several cases regressed in `Routing open`.
  - Current conclusion: first connectivity penalty design/weighting needs retuning or redesign before continuing.
- Disabled connectivity weight (`conn_weight_max=0`) to continue with later phases.
- Implemented edge/congestion trial:
  - Edge strengthening:
    - Increased edge penalty coefficient.
    - Edge penalty now takes best legal location option per EDGE block (instead of fixed active option only).
  - Congestion proxy:
    - Added RUDY-like coarse-bin congestion penalty in SA cost.
    - Added time-ramped congestion weight schedule in SA.
- Re-ran same short benchmark (`20s`, same 5 seeds):
  - Fails: `1,4,5,4,12`
  - Penalties: `0,0,2,5,0`
  - Routing-open counts: `0,4,5,3,9`
  - Edge fails: `1,0,0,1,3`
- Result summary:
  - Mid-size cases (`n=16,18`) improved fail count and routing-open count.
  - Largest case (`n=20`) still high fail count and edge fails.
  - Congestion proxy helps some routability but currently increases penalties on some instances; weight tuning is still needed.
- Added CLI switches for fast mode toggling without code edits:
  - `--time <sec>`
  - `--enable-congestion | --disable-congestion`
  - `--cong-weight <v> --cong-bins <n>`
  - `--enable-connectivity | --disable-connectivity`
  - `--conn-weight <v> --conn-topk <k>`
  - `--edge-penalty <v>`
  - `--edge-best-location | --edge-active-location`
- Time-limit comparison (same 5-case short suite, congestion on, connectivity off):
  - `20s`: fails = `0,0,13,5,7`
  - `40s`: fails = `1,0,9,6,16`
  - Observation: longer runtime does not monotonically improve quality in current stochastic multi-start setting; some cases improve while others regress.
- Implemented final repair pass with "hard-at-edge" push:
  - Added optional post-route repair iterations that nudge `HARD_MACRO` blocks outward.
  - Each repair iteration re-runs channelization + routing.
  - Accept criterion: improve `routing open` count; tie-break with fewer channel overflows.
  - Added CLI toggles: `--enable-repair/--disable-repair`, `--repair-iters`, `--repair-step-ratio`.
- Quick A/B on 5 short cases (`20s`, congestion on, connectivity off):
  - Repair OFF fails: `0,3,9,5,3`
  - Repair ON fails:  `0,4,5,6,3`
  - Observation: helps some harder routing-open cases (notably `n=16`), but can regress some mid cases (`n=14,18`).
  - Conclusion: repair pass is useful but still stochastic; should be used as optional mode and tuned by testcase profile.

### 2026-05-29 — Sprint A + B (minimal, library-aligned)

**Sprint A — Observability + gated repair**
- `route_all(..., &failed_conn)` records connection indices that fail to route.
- `log_failed_pairs()` prints `routing_open=N/M` and up to 8 block-pair names.
- `[Route] after main flow` log after final routing in each worker.
- Repair (geo + SA) runs only when `routing_open >= repair_trigger_open` (default `1`).
- `count_routing_open()` now uses `failed_conn.size()` (consistent with router).

**Sprint B — Repair-SA (hard-at-center via sequence pair)**
- `hard_center_penalty()` in SA cost (ramps with temperature fraction in `sa_optimizer.h`).
- `run_repair_sa()`: short SA pass (`repair_sa_time`, default `12s`) with `hard_center_weight_max` (default `3e5`), connectivity off, half congestion if enabled.
- Pipeline when repair triggered: **Repair-SA → reroute → log** (geo push removed 2026-05-29).
- CLI:
  - `--enable-repair` / `--disable-repair` (Repair-SA; aliases `--enable-repair-sa` / `--disable-repair-sa`)
  - `--repair-sa-time <sec>`
  - `--hard-center-weight <value>`
  - `--repair-trigger-open <int>`

Recommended command:

```bash
./bin/solver in.csv out.cfg --time 30 \
  --disable-connectivity \
  --enable-congestion --cong-weight 3500 --cong-bins 10 \
  --edge-best-location --edge-penalty 2000000 \
  --enable-repair --repair-sa-time 12 --hard-center-weight 300000 \
  --repair-trigger-open 1
```

Smoke test: `testcase_official.csv` @ 25s → evaluator **0 FAIL** (repair skipped when `routing_open=0`).

### 2026-05-29 — Remove geo repair (Repair-Geo)

- Removed `run_final_repair()` (direct HARD_MACRO nudge toward outline).
- `--enable-repair` now controls **Repair-SA only**; removed `--repair-iters`, `--repair-step-ratio`.

