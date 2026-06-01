#!/usr/bin/env python3
"""A/B benchmark: Repair-SA on vs off with multi-seed medians."""
import argparse
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "bin" / "solver"
EVAL = ROOT / "script" / "evaluator.py"
OUT_DIR = ROOT / "benchcases" / "ab_repair_sa"


def run_solver(csv_path: Path, cfg_path: Path, time_sec: float, repair_sa: bool, seed: int) -> None:
    cmd = [
        str(SOLVER),
        str(csv_path),
        str(cfg_path),
        "--time",
        str(time_sec),
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
        "--seed-base",
        str(seed),
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
        "outline_fail": 0,
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
        elif "Outline constraint" in line:
            m["outline_fail"] += 1
    return m


def median_metrics(metrics: list[dict]) -> dict:
    out = {}
    for k in ["fails", "routing_open", "edge_fail", "outline_fail", "penalties", "cost"]:
        vals = [m[k] for m in metrics if m[k] is not None]
        out[k] = statistics.median(vals) if vals else None
    return out


def case_mode_better(a: dict, b: dict) -> bool:
    return (
        a["fails"],
        a["routing_open"],
        a["penalties"],
        a["cost"] if a["cost"] is not None else 1e30,
    ) < (
        b["fails"],
        b["routing_open"],
        b["penalties"],
        b["cost"] if b["cost"] is not None else 1e30,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", type=float, default=25.0)
    ap.add_argument("--cases", nargs="*", help="CSV paths (default: bench + auto)")
    ap.add_argument("--seeds", nargs="*", type=int, default=[1001, 2002, 3003], help="Seed bases for repeated runs")
    ap.add_argument("--penalty-tolerance", type=float, default=1.0, help="Allowed median penalty increase for gate pass")
    args = ap.parse_args()

    if args.cases:
        cases = [Path(c) for c in args.cases]
    else:
        cases = sorted((ROOT / "benchcases").glob("case_n*.csv"))
        cases.append(ROOT / "testcase_auto.csv")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str, int, dict]] = []

    print(f"{'case':<22} {'mode':<8} {'seed':>6} {'fails':>5} {'r_open':>6} {'edge':>5} {'outl':>5} {'pen':>4} {'cost':>14}")
    print("-" * 98)

    for csv_path in cases:
        name = csv_path.stem
        for mode, sa_on in [("sa_on", True), ("sa_off", False)]:
            for seed in args.seeds:
                cfg = OUT_DIR / f"{name}_{mode}_s{seed}.cfg"
                try:
                    run_solver(csv_path, cfg, args.time, sa_on, seed)
                    m = parse_evaluator(csv_path, cfg)
                except Exception as e:
                    print(f"{name:<22} {mode:<8} {seed:>6} ERROR {e}")
                    continue
                rows.append((name, mode, seed, m))
                print(
                    f"{name:<22} {mode:<8} {seed:>6} {m['fails']:>5} {m['routing_open']:>6} "
                    f"{m['edge_fail']:>5} {m['outline_fail']:>5} {m['penalties']:>4} {m['cost']:>14.2f}"
                )

    print("\n=== Median Summary by Case ===")
    by_case_mode: dict[str, dict[str, list[dict]]] = {}
    for name, mode, _seed, m in rows:
        by_case_mode.setdefault(name, {}).setdefault(mode, []).append(m)

    wins_on = wins_off = ties = 0
    gate_pass = gate_total = 0
    for name, d in sorted(by_case_mode.items()):
        if "sa_on" not in d or "sa_off" not in d:
            continue
        med_on = median_metrics(d["sa_on"])
        med_off = median_metrics(d["sa_off"])

        if case_mode_better(med_on, med_off):
            wins_on += 1
            verdict = "SA_on better"
        elif case_mode_better(med_off, med_on):
            wins_off += 1
            verdict = "SA_off better"
        else:
            ties += 1
            verdict = "tie"

        gate_total += 1
        gate_ok = (
            med_on["routing_open"] <= med_off["routing_open"]
            and med_on["fails"] <= med_off["fails"]
            and (med_on["penalties"] - med_off["penalties"]) <= args.penalty_tolerance
        )
        if gate_ok:
            gate_pass += 1

        print(
            f"{name}: median fails {med_on['fails']} vs {med_off['fails']}, "
            f"open {med_on['routing_open']} vs {med_off['routing_open']}, "
            f"pen {med_on['penalties']} vs {med_off['penalties']} -> {verdict}; "
            f"gate={'PASS' if gate_ok else 'FAIL'}"
        )

    print(
        f"\nSummary: SA_on wins {wins_on}, SA_off wins {wins_off}, ties {ties} "
        f"(cases={len(by_case_mode)}). Gate pass {gate_pass}/{gate_total} "
        f"with penalty_tolerance={args.penalty_tolerance}."
    )


if __name__ == "__main__":
    main()
