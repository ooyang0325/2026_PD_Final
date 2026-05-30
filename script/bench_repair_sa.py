#!/usr/bin/env python3
"""A/B benchmark: Repair-SA (hard-center) on vs off."""
import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "bin" / "solver"
EVAL = ROOT / "script" / "evaluator.py"
OUT_DIR = ROOT / "benchcases" / "ab_repair_sa"


def run_solver(csv_path: Path, cfg_path: Path, time_sec: float, repair_sa: bool) -> None:
    cmd = [
        str(SOLVER),
        str(csv_path),
        str(cfg_path),
        "--time",
        str(time_sec),
        "--disable-connectivity",
        "--enable-congestion",
        "--cong-weight",
        "3500",
        "--cong-bins",
        "10",
        "--edge-best-location",
        "--edge-penalty",
        "2000000",
        "--repair-trigger-open",
        "1",
        "--repair-sa-time",
        "12",
        "--hard-center-weight",
        "300000",
    ]
    if repair_sa:
        cmd.append("--enable-repair")
    else:
        cmd.append("--disable-repair")

    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_evaluator(csv_path: Path, cfg_path: Path) -> dict:
    p = subprocess.run(
        [sys.executable, str(EVAL), str(csv_path), str(cfg_path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    text = p.stdout + p.stderr
    m = {
        "fails": None,
        "penalties": None,
        "cost": None,
        "routing_open": 0,
        "edge_fail": 0,
    }
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
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", type=float, default=25.0)
    ap.add_argument("--cases", nargs="*", help="CSV paths (default: bench + auto)")
    args = ap.parse_args()

    if args.cases:
        cases = [Path(c) for c in args.cases]
    else:
        cases = sorted((ROOT / "benchcases").glob("case_n*.csv"))
        cases.append(ROOT / "testcase_auto.csv")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    print(f"{'case':<28} {'mode':<12} {'fails':>5} {'r_open':>6} {'edge':>5} {'pen':>4} {'cost':>14}")
    print("-" * 80)

    for csv_path in cases:
        name = csv_path.stem
        for mode, sa_on in [("sa_on", True), ("sa_off", False)]:
            cfg = OUT_DIR / f"{name}_{mode}.cfg"
            try:
                run_solver(csv_path, cfg, args.time, sa_on)
                m = parse_evaluator(csv_path, cfg)
            except Exception as e:
                print(f"{name:<28} {mode:<12} ERROR {e}")
                continue
            rows.append((name, mode, m))
            print(
                f"{name:<28} {mode:<12} {m['fails']:>5} {m['routing_open']:>6} "
                f"{m['edge_fail']:>5} {m['penalties']:>4} {m['cost']:>14.2f}"
            )

    print("\n=== Pairwise (sa_on vs sa_off) ===")
    by_name = {}
    for name, mode, m in rows:
        by_name.setdefault(name, {})[mode] = m
    wins_on = wins_off = ties = 0
    for name, d in sorted(by_name.items()):
        if "sa_on" not in d or "sa_off" not in d:
            continue
        a, b = d["sa_on"], d["sa_off"]
        if a["fails"] < b["fails"]:
            wins_on += 1
            verdict = "SA better (fewer fails)"
        elif a["fails"] > b["fails"]:
            wins_off += 1
            verdict = "SA worse"
        elif (a["routing_open"], a["penalties"], a["cost"] or 1e30) < (
            b["routing_open"],
            b["penalties"],
            b["cost"] or 1e30,
        ):
            wins_on += 1
            verdict = "SA better (tie fails, better secondary)"
        elif (a["routing_open"], a["penalties"], a["cost"] or 1e30) > (
            b["routing_open"],
            b["penalties"],
            b["cost"] or 1e30,
        ):
            wins_off += 1
            verdict = "SA worse (tie fails)"
        else:
            ties += 1
            verdict = "tie"
        print(
            f"  {name}: fails {a['fails']} vs {b['fails']}, "
            f"r_open {a['routing_open']} vs {b['routing_open']} -> {verdict}"
        )
    print(f"\nSummary: SA_on wins {wins_on}, SA_off wins {wins_off}, ties {ties} (of {len(by_name)} cases)")


if __name__ == "__main__":
    main()
