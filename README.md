# ICCAD 2026 Problem E Floorplanner

C++ floorplanner and global router for ICCAD 2026 Problem E. The solver reads a CSV, runs simulated annealing with feedthrough-aware sizing and routing, and writes a .cfg output.

## Build

```bash
make
```

Binary is created at `bin/solver`.

## Run

```bash
./bin/solver <input.csv> <output.cfg> [time_limit_sec]
```

- `input.csv`: official CSV format(probably, they only provide .xlsx file)
- `output.cfg`: result file in problem-specified format
- `time_limit_sec`: optional, default `7000`

## Output

The solver writes:
- `OUTLINE` (current packed width/height)
- `BLOCK` placements
- `CHANNEL` rectangles
- `PATH` routing segments

## Scripts

- `script/generator.py`: generate random CSV testcases for local testing.
  ```bash
  python3 script/generator.py -n <# of blocks> -u <utilization> -a <alpha> -s <seed> -o <output filename>
  ```

- `script/evaluator.py`: check a CSV/CFG pair and report violations and score.
  ```bash
  python3 script/evaluator.py <input.csv> <output.cfg>
  ```

## Notes

- Parallel multi-start search is enabled; the solver runs multiple independent searches and keeps the best result.
- CPU worker count is derived from `hardware_concurrency()` as `max(1, cpu_count - 2)`.
