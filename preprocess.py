#!/usr/bin/env python3
"""Convert ICCAD Problem E xlsx to clean .in format for C++ solver."""
import sys, openpyxl, math

def convert(xlsx_path, out_path):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))

    blocks = []
    outline_w = outline_h = 0.0
    alpha = 1.0
    conn_matrix = []
    block_names = []

    sec = None
    for row in rows:
        if not any(v is not None for v in row):
            continue
        c0 = str(row[0]).strip() if row[0] is not None else ''

        if c0 == 'BLOCK':
            sec = 'BLOCK'
            continue
        elif c0 == 'OUTLINE':
            sec = 'OUTLINE'
            continue
        elif c0.replace(' ', '') in ('α', 'alpha', 'a'):
            try: alpha = float(row[1])
            except: pass
            continue
        elif 'CONN' in c0 or 'MATRIX' in c0 or 'INTERFACE' in c0:
            sec = 'CONN'
            continue

        if sec == 'BLOCK':
            if c0 in ('', 'MAX', 'FT CONVERSION') or c0.startswith('<=') or c0.startswith('>'):
                continue
            if not c0.startswith('BLK'):
                continue
            name = c0
            area = float(row[1]) if row[1] is not None else 0.0
            width = float(row[2]) if row[2] is not None else 0.0
            height = float(row[3]) if row[3] is not None else 0.0
            ar_raw = str(row[4]).strip() if row[4] is not None else '1'
            if ',' in ar_raw:
                parts = ar_raw.split(',')
                min_ar, max_ar = float(parts[0]), float(parts[1])
            else:
                min_ar = max_ar = float(ar_raw)
            btype = str(row[5]).strip() if row[5] is not None else 'SOFT'
            loc_raw = str(row[6]).strip() if row[6] is not None else ''
            locs = [l.strip() for l in loc_raw.split(',') if l.strip()] if loc_raw else []
            ft = [float(row[7] or 0.2), float(row[8] or 0.4),
                  float(row[9] or 0.8), float(row[10] or 1.0)]
            # Normalize: if values > 1, divide by 100
            ft = [v/100.0 if v > 1.0 else v for v in ft]
            if width == 0 and height == 0 and area > 0:
                width = height = math.sqrt(area)
            elif area == 0 and width > 0 and height > 0:
                area = width * height
            blocks.append((name, area, width, height, min_ar, max_ar, btype, locs, ft))
            block_names.append(name)

        elif sec == 'OUTLINE':
            if c0 == 'MAX':
                try: outline_w = float(row[1])
                except: pass
                try: outline_h = float(row[2])
                except: pass

        elif sec == 'CONN':
            if row[0] is None or not str(row[0]).startswith('BLK'):
                continue
            from_name = str(row[0]).strip()
            row_data = []
            for i, bn in enumerate(block_names):
                v = row[i+1] if (i+1) < len(row) else None
                row_data.append(int(v) if v is not None else 0)
            conn_matrix.append((from_name, row_data))

    # Write .in file
    with open(out_path, 'w') as f:
        f.write(f"OUTLINE {outline_w:.6f} {outline_h:.6f}\n")
        f.write(f"ALPHA {alpha}\n")
        f.write(f"BLOCKS {len(blocks)}\n")
        for (name, area, w, h, min_ar, max_ar, btype, locs, ft) in blocks:
            loc_str = ','.join(locs) if locs else 'NONE'
            f.write(f"{name} {area:.6f} {w:.6f} {h:.6f} "
                    f"{min_ar} {max_ar} {btype} {loc_str} "
                    f"{ft[0]} {ft[1]} {ft[2]} {ft[3]}\n")
        f.write("END_BLOCKS\n")

        # Write connection matrix as flat pairs (from, to, nets) — symmetric, store all
        f.write("CONNECTIONS\n")
        for from_name, row_data in conn_matrix:
            for to_idx, nets in enumerate(row_data):
                if nets > 0:
                    f.write(f"{from_name} {block_names[to_idx]} {nets}\n")
        f.write("END_CONNECTIONS\n")

    print(f"Written: {out_path}")
    print(f"Outline: {outline_w:.1f} x {outline_h:.1f}, alpha={alpha}")
    print(f"Blocks: {len(blocks)}")
    total_conn = sum(1 for _, rd in conn_matrix for v in rd if v > 0)
    print(f"Connections (entries): {total_conn}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: preprocess.py <input.xlsx> <output.in>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
