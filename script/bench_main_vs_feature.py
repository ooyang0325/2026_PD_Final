#!/usr/bin/env python3
"""Compare main (raw/ws) vs feature branch on shared small cases."""
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVAL = ROOT / "script" / "evaluator.py"
OUT = ROOT / "benchcases" / "compare" / "runs"
OUT.mkdir(parents=True, exist_ok=True)

SOLVERS = {
    "main_raw": ROOT.parent / "2026_PD_Final_compare_main" / "bin" / "solver_main_raw",
    "main_ws": ROOT.parent / "2026_PD_Final_compare_main" / "bin" / "solver_main_ws",
    "feature": ROOT / "bin" / "solver_feature",
}

CASES = [
    ROOT / "testcase_auto.csv",
    ROOT / "testcase_official.csv",
    ROOT / "benchcases" / "compare" / "case_n12.csv",
    ROOT / "benchcases" / "compare" / "case_n20.csv",
    ROOT / "benchcases" / "compare" / "case_n28.csv",
]
SEEDS = [801, 802]
TIME = 20.0


def run_solver(solver: Path, csv_path: Path, cfg_path: Path, mode: str) -> None:
    cmd = [str(solver), str(csv_path), str(cfg_path), str(TIME)]
    if mode == "feature":
        cmd += [
            "--enable-congestion",
            "--cong-weight", "3500", "--cong-bins", "10",
            "--edge-best-location", "--edge-penalty", "2000000",
            "--enable-repair", "--repair-sa-time", "12",
            "--repair-trigger-open", "1",
            "--hard-center-weight", "100000",
            "--repair-hard-center-scale", "0.35",
            "--seed-base", "801",
        ]
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_eval(text: str) -> dict:
    m = {"fails": None, "penalties": None, "cost": None, "routing_open": 0, "outline_fail": 0}
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
    for name, path in SOLVERS.items():
        if not path.exists():
            print(f"Missing solver: {path}", file=sys.stderr)
            sys.exit(1)

    rows = []
    print(f"{'case':<20} {'mode':<10} {'seed':>5} {'fail':>4} {'open':>4} {'pen':>4} {'cost':>14}")
    print("-" * 72)
    for case in CASES:
        cname = case.stem
        for mode in ["main_raw", "main_ws", "feature"]:
            for seed in SEEDS:
                cfg = OUT / f"{cname}_{mode}_s{seed}.cfg"
                print(f"Running {cname} {mode} seed={seed}...", flush=True)
                run_solver(SOLVERS[mode], case, cfg, mode)
                p = subprocess.run(
                    [sys.executable, str(EVAL), str(case), str(cfg)],
                    cwd=ROOT, capture_output=True, text=True,
                )
                ev = parse_eval(p.stdout + p.stderr)
                rows.append((cname, mode, seed, ev))
                print(f"{cname:<20} {mode:<10} {seed:>5} {ev['fails']:>4} {ev['routing_open']:>4} "
                      f"{ev['penalties']:>4} {ev['cost']:>14.1f}")

    print("\n=== Median summary ===")
    by = {}
    for cname, mode, _seed, ev in rows:
        by.setdefault(cname, {}).setdefault(mode, []).append(ev)

    modes = ["main_raw", "main_ws", "feature"]
    hdr = f"{'case':<20}" + "".join(f"{m:>16}" for m in modes) + f"{'best':>10}"
    print(hdr)
    print("-" * len(hdr))
    wins = {m: 0 for m in modes}
    for cname in sorted(by):
        parts = []
        best_mode = None
        best_key = None
        for m in modes:
            evs = by[cname].get(m, [])
            if not evs:
                parts.append(f"{'n/a':>16}")
                continue
            mf = statistics.median([e["fails"] for e in evs])
            mc = statistics.median([e["cost"] for e in evs if e["cost"] is not None])
            mp = statistics.median([e["penalties"] for e in evs])
            parts.append(f"{mf:>3.0f}/{mp:>2.0f}/{mc/1e6:>5.1f}M")
            key = (mf, mp, mc if mf == 0 else 1e30)
            if best_key is None or key < best_key:
                best_key = key
                best_mode = m
        if best_mode:
            wins[best_mode] += 1
        print(f"{cname:<20}" + "".join(parts) + f"{best_mode or '':>10}")

    print(f"\nLexicographic wins (fail, pen, cost): {wins}")
    print("Format per cell: fail/pen/cost(M)")


if __name__ == "__main__":
    main()
