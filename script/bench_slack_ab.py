#!/usr/bin/env python3
"""A/B: slack on (current defaults) vs slack off (pre-slack penalty profile)."""
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "bin" / "solver"
EVAL = ROOT / "script" / "evaluator.py"
OUT = ROOT / "benchcases" / "slack_ab"
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


def run_solver(csv_path: Path, cfg_path: Path, seed: int, slack_on: bool) -> None:
    cmd = [
        str(SOLVER), str(csv_path), str(cfg_path), str(TIME),
        "--enable-congestion",
        "--cong-weight", "3500", "--cong-bins", "10",
        "--edge-best-location", "--edge-penalty", "2000000",
        "--enable-repair", "--repair-sa-time", "12",
        "--repair-trigger-open", "1",
        "--seed-base", str(seed),
    ]
    if slack_on:
        cmd += ["--hard-center-weight", "100000", "--repair-hard-center-scale", "0.35"]
    else:
        cmd += [
            "--sa1-slack-weight", "0", "--sa2-slack-weight", "0",
            "--repair-slack-weight", "0",
            "--hard-center-weight", "300000", "--repair-hard-center-scale", "0.70",
        ]
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_eval(text: str) -> dict:
    m = {"fails": None, "penalties": None, "cost": None,
         "routing_open": 0, "edge_fail": 0, "outline_fail": 0}
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


def med(vals):
    return statistics.median(vals) if vals else None


def main():
    rows = []
    print(f"{'case':<22} {'mode':<10} {'seed':>5} {'fail':>4} {'open':>4} {'outl':>4} {'pen':>4} {'cost':>14}")
    print("-" * 78)
    for case in CASES:
        if not case.exists():
            continue
        name = case.stem
        for mode, slack_on in [("slack_on", True), ("slack_off", False)]:
            for seed in SEEDS:
                cfg = OUT / f"{name}_{mode}_s{seed}.cfg"
                run_solver(case, cfg, seed, slack_on)
                p = subprocess.run(
                    [sys.executable, str(EVAL), str(case), str(cfg)],
                    cwd=ROOT, capture_output=True, text=True,
                )
                ev = parse_eval(p.stdout + p.stderr)
                rows.append((name, mode, seed, ev))
                print(f"{name:<22} {mode:<10} {seed:>5} {ev['fails']:>4} {ev['routing_open']:>4} "
                      f"{ev['outline_fail']:>4} {ev['penalties']:>4} {ev['cost']:>14.1f}")

    print("\n=== Median by case ===")
    by = {}
    for name, mode, _seed, ev in rows:
        by.setdefault(name, {}).setdefault(mode, []).append(ev)

    cost_wins_on = cost_wins_off = cost_tie = 0
    pass_both = []
    print(f"{'case':<22} {'fail_on':>7} {'fail_off':>8} {'cost_on':>14} {'cost_off':>14} {'delta%':>8} verdict")
    print("-" * 90)
    for name in sorted(by):
        on = by[name].get("slack_on", [])
        off = by[name].get("slack_off", [])
        mf_on, mf_off = med([e["fails"] for e in on]), med([e["fails"] for e in off])
        mc_on = med([e["cost"] for e in on if e["cost"] is not None])
        mc_off = med([e["cost"] for e in off if e["cost"] is not None])
        mp_on = med([e["penalties"] for e in on])
        mp_off = med([e["penalties"] for e in off])

        # Cost verdict only meaningful when both pass (fails=0)
        if mf_on == 0 and mf_off == 0 and mc_on and mc_off:
            delta = (mc_on - mc_off) / mc_off * 100.0
            if mc_on < mc_off * 0.995:
                v = "slack_on cheaper"
                cost_wins_on += 1
            elif mc_off < mc_on * 0.995:
                v = "slack_off cheaper"
                cost_wins_off += 1
            else:
                v = "~tie"
                cost_tie += 1
            pass_both.append((name, mc_on, mc_off, delta, mp_on, mp_off))
        elif mf_on < mf_off:
            v = "slack_on fewer fails"
        elif mf_off < mf_on:
            v = "slack_off fewer fails"
        else:
            delta = ((mc_on - mc_off) / mc_off * 100.0) if mc_on and mc_off else float("nan")
            v = "both fail (cost N/A)" if mf_on > 0 else "~tie"
            if mf_on == 0 and mc_on and mc_off:
                if mc_on < mc_off * 0.995:
                    cost_wins_on += 1
                    v = "slack_on cheaper"
                elif mc_off < mc_on * 0.995:
                    cost_wins_off += 1
                    v = "slack_off cheaper"
        delta_s = ""
        if mf_on == 0 and mf_off == 0 and mc_on and mc_off:
            delta_s = f"{(mc_on-mc_off)/mc_off*100:+.1f}%"
        print(f"{name:<22} {mf_on:>7.0f} {mf_off:>8.0f} {mc_on:>14.1f} {mc_off:>14.1f} {delta_s:>8} {v}")

    if pass_both:
        avg_delta = statistics.mean([d for _, _, _, d, _, _ in pass_both])
        print(f"\nAmong {len(pass_both)} cases with fails=0 in both modes:")
        print(f"  cost: slack_on wins {cost_wins_on}, slack_off wins {cost_wins_off}, tie {cost_tie}")
        print(f"  mean median cost delta (on vs off): {avg_delta:+.2f}%")


if __name__ == "__main__":
    main()
