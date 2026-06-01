#!/usr/bin/env python3
"""Benchmark P0+P1 vs previous slack-only baseline (seeds 801-803)."""
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "bin" / "solver"
EVAL = ROOT / "script" / "evaluator.py"
OUT = ROOT / "benchcases" / "p0p1_run"
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

# Median from prior slack_on run (bench_slack_ab.py, before P0/P1)
BASELINE = {
    "test1": {"fails": 3.0, "cost": 1942626486.6, "outline": 1.0},
    "testcase_auto": {"fails": 0.0, "cost": 71213305.1, "outline": 0.0},
    "testcase_official": {"fails": 0.0, "cost": 3877019.8, "outline": 0.0},
    "case_n12_s101": {"fails": 0.0, "cost": 8139904.9, "outline": 0.0},
    "case_n16_s303": {"fails": 0.0, "cost": 21741295.1, "outline": 0.0},
    "case_n20_s505": {"fails": 0.0, "cost": 60882017.1, "outline": 0.0},
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
    m = {"fails": None, "penalties": None, "cost": None,
         "routing_open": 0, "outline_fail": 0}
    for line in text.splitlines():
        if "Total Fails" in line:
            m["fails"] = int(re.search(r":\s*(\d+)", line).group(1))
        elif "Total Penalties" in line:
            m["penalties"] = int(re.search(r":\s*(\d+)", line).group(1))
        elif "Total Cost" in line:
            m["cost"] = float(re.search(r":\s*([\d.]+)", line).group(1))
        elif "Routing open" in line:
            m["routing_open"] += 1
        elif "Outline constraint" in line:
            m["outline_fail"] += 1
    return m


def main():
    rows = []
    print(f"{'case':<22} {'seed':>5} {'fail':>4} {'open':>4} {'outl':>4} {'pen':>4} {'cost':>14}")
    print("-" * 70)
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
            print(f"{name:<22} {seed:>5} {ev['fails']:>4} {ev['routing_open']:>4} "
                  f"{ev['outline_fail']:>4} {ev['penalties']:>4} {ev['cost']:>14.1f}")

    print("\n=== P0+P1 vs prior slack-only baseline (median, seeds 801-803) ===")
    by = {}
    for name, _seed, ev in rows:
        by.setdefault(name, []).append(ev)

    hdr = f"{'case':<22} {'fail':>5} {'bfail':>5} {'cost':>14} {'bcost':>14} {'Δcost':>8} note"
    print(hdr)
    print("-" * len(hdr))
    deltas = []
    for name in sorted(by):
        evs = by[name]
        mf = statistics.median([e["fails"] for e in evs])
        mc = statistics.median([e["cost"] for e in evs if e["cost"] is not None])
        mo = statistics.median([e["outline_fail"] for e in evs])
        base = BASELINE.get(name, {})
        bf = base.get("fails", float("nan"))
        bc = base.get("cost", float("nan"))
        delta = (mc - bc) / bc * 100 if bc else float("nan")
        if mf == 0 and bf == 0 and bc:
            deltas.append(delta)
        note = []
        if mf < bf:
            note.append("fail↓")
        elif mf > bf:
            note.append("fail↑")
        if mf == 0 and bf == 0 and delta < -1:
            note.append("cost↓")
        elif mf == 0 and bf == 0 and delta > 1:
            note.append("cost↑")
        print(f"{name:<22} {mf:>5.0f} {bf:>5.0f} {mc:>14.1f} {bc:>14.1f} {delta:>+7.1f}% {','.join(note)}")

    if deltas:
        print(f"\nMean cost delta (both PASS): {statistics.mean(deltas):+.2f}%")


if __name__ == "__main__":
    main()
