#!/usr/bin/env python3
"""Benchmark Tier-1 cost opts vs prior P0+P1 baseline (seeds 801-803)."""
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "bin" / "solver"
EVAL = ROOT / "script" / "evaluator.py"
OUT = ROOT / "benchcases" / "tier1_run"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    ROOT / "test1.csv",
    ROOT / "testcase_auto.csv",
    ROOT / "testcase_official.csv",
    ROOT / "benchcases" / "case_n12_s101.csv",
    ROOT / "benchcases" / "case_n16_s303.csv",
    ROOT / "benchcases" / "case_n20_s505.csv",
]
SEEDS = [801, 802, 803]
TIME = 25.0

# P0+P1 baseline medians (bench_p0p1.py)
BASELINE = {
    "test1": {"fails": 3.0, "cost": 1828220701.0},
    "testcase_auto": {"fails": 0.0, "cost": 118762698.8},
    "testcase_official": {"fails": 0.0, "cost": 3875915.6},
    "case_n12_s101": {"fails": 0.0, "cost": 8120704.1},
    "case_n16_s303": {"fails": 0.0, "cost": 19250442.2},
    "case_n20_s505": {"fails": 0.0, "cost": 55770567.9},
}


def run_solver(csv_path: Path, cfg_path: Path, seed: int) -> None:
    cmd = [
        str(SOLVER), str(csv_path), str(cfg_path), str(TIME),
        "--enable-congestion",
        "--cong-weight", "3500", "--cong-bins", "10",
        "--edge-best-location", "--edge-penalty", "2000000",
        "--enable-repair", "--repair-sa-time", "12",
        "--repair-trigger-open", "1",
        "--hard-center-weight", "100000", "--repair-hard-center-scale", "0.35",
        "--seed-base", str(seed),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_eval(text: str) -> dict:
    m = {"fails": None, "penalties": None, "cost": None, "area": None, "hpwl": None,
         "routing_open": 0, "outline_fail": 0}
    for line in text.splitlines():
        if "Total Fails" in line:
            m["fails"] = int(re.search(r":\s*(\d+)", line).group(1))
        elif "Total Penalties" in line:
            m["penalties"] = int(re.search(r":\s*(\d+)", line).group(1))
        elif "Area Cost" in line:
            m["area"] = float(re.search(r":\s*([\d.]+)", line).group(1))
        elif "HPWL" in line:
            m["hpwl"] = float(re.search(r":\s*([\d.]+)", line).group(1))
        elif "Total Cost" in line:
            m["cost"] = float(re.search(r":\s*([\d.]+)", line).group(1))
        elif "Routing open" in line:
            m["routing_open"] += 1
        elif "Outline constraint" in line:
            m["outline_fail"] += 1
    return m


def main():
    rows = []
    for case in CASES:
        name = case.stem
        for seed in SEEDS:
            cfg = OUT / f"{name}_s{seed}.cfg"
            print(f"Running {name} seed={seed}...", flush=True)
            run_solver(case, cfg, seed)
            p = subprocess.run(
                [sys.executable, str(EVAL), str(case), str(cfg)],
                cwd=ROOT, capture_output=True, text=True,
            )
            ev = parse_eval(p.stdout + p.stderr)
            rows.append((name, seed, ev))
            print(f"  fail={ev['fails']} cost={ev['cost']:.1f} pen={ev['penalties']}")

    print("\n=== Tier1 vs P0+P1 baseline (median) ===")
    by = {}
    for name, _s, ev in rows:
        by.setdefault(name, []).append(ev)

    cost_deltas = []
    print(f"{'case':<22} {'fail':>5} {'bfail':>5} {'cost':>14} {'bcost':>14} {'Δcost':>8} note")
    print("-" * 78)
    for name in sorted(by):
        evs = by[name]
        mf = statistics.median([e["fails"] for e in evs])
        mc = statistics.median([e["cost"] for e in evs if e["cost"] is not None])
        base = BASELINE.get(name, {})
        bf = base.get("fails", float("nan"))
        bc = base.get("cost", float("nan"))
        delta = (mc - bc) / bc * 100 if bc else float("nan")
        note = []
        if mf < bf: note.append("fail↓")
        elif mf > bf: note.append("fail↑")
        if mf == 0 and bf == 0 and delta < -1: note.append("cost↓")
        elif mf == 0 and bf == 0 and delta > 1: note.append("cost↑")
        if mf == 0 and bf == 0 and bc:
            cost_deltas.append(delta)
        print(f"{name:<22} {mf:>5.0f} {bf:>5.0f} {mc:>14.1f} {bc:>14.1f} {delta:>+7.1f}% {','.join(note)}")

    if cost_deltas:
        print(f"\nMean cost delta (both PASS): {statistics.mean(cost_deltas):+.2f}%")


if __name__ == "__main__":
    main()
