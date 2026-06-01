#!/usr/bin/env python3
"""Quick multi-case benchmark with evaluator metrics."""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "bin" / "solver"
EVAL = ROOT / "script" / "evaluator.py"
OUT = ROOT / "benchcases" / "slack_run"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    ROOT / "test1.csv",
    ROOT / "testcase_auto.csv",
    ROOT / "testcase_official.csv",
    ROOT / "benchcases" / "case_n12_s101.csv",
    ROOT / "benchcases" / "case_n16_s303.csv",
    ROOT / "benchcases" / "case_n20_s505.csv",
]
SEEDS = [701, 702, 703]
TIME = 25.0


def run_solver(csv_path: Path, cfg_path: Path, seed: int) -> None:
    cmd = [
        str(SOLVER),
        str(csv_path),
        str(cfg_path),
        str(TIME),
        "--enable-congestion",
        "--cong-weight", "3500",
        "--cong-bins", "10",
        "--edge-best-location",
        "--edge-penalty", "2000000",
        "--enable-repair",
        "--repair-sa-time", "12",
        "--repair-trigger-open", "1",
        "--hard-center-weight", "100000",
        "--repair-hard-center-scale", "0.35",
        "--seed-base", str(seed),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_eval(text: str) -> dict:
    m = {"fails": None, "penalties": None, "cost": None, "routing_open": 0,
         "edge_fail": 0, "outline_fail": 0, "overlap_fail": 0}
    for line in text.splitlines():
        if "Total Fails" in line:
            m["fails"] = int(re.search(r":\s*(\d+)", line).group(1))
        elif "Total Penalties" in line:
            m["penalties"] = int(re.search(r":\s*(\d+)", line).group(1))
        elif "Total Cost" in line:
            m["cost"] = float(re.search(r":\s*([\d.]+)", line).group(1))
        elif "Routing open" in line:
            m["routing_open"] += 1
        elif "Edge constraint" in line:
            m["edge_fail"] += 1
        elif "Outline constraint" in line:
            m["outline_fail"] += 1
        elif "Overlap" in line and "constraint" in line.lower():
            m["overlap_fail"] += 1
    return m


def main():
    rows = []
    for case in CASES:
        if not case.exists():
            print(f"SKIP missing {case}", file=sys.stderr)
            continue
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
            print(f"  fails={ev['fails']} open={ev['routing_open']} outline={ev['outline_fail']} "
                  f"edge={ev['edge_fail']} pen={ev['penalties']} cost={ev['cost']}")

    print("\n=== Summary (median over seeds) ===")
    import statistics
    by_case = {}
    for name, seed, ev in rows:
        by_case.setdefault(name, []).append(ev)
    hdr = f"{'case':<22} {'fails':>6} {'open':>5} {'outline':>8} {'edge':>5} {'pen':>6} {'cost':>14}"
    print(hdr)
    print("-" * len(hdr))
    for name in sorted(by_case):
        evs = by_case[name]
        med = lambda k: statistics.median([e[k] for e in evs if e[k] is not None])
        print(f"{name:<22} {med('fails'):>6.1f} {med('routing_open'):>5.1f} "
              f"{med('outline_fail'):>8.1f} {med('edge_fail'):>5.1f} "
              f"{med('penalties'):>6.1f} {med('cost'):>14.1f}")


if __name__ == "__main__":
    main()
