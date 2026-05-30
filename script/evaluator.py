import sys
import math
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Evaluator:
    def __init__(self, csv_file, cfg_file):
        self.csv_file = csv_file
        self.cfg_file = cfg_file
        self.blocks_info = {}
        self.alpha = 1.0
        self.max_outline = (0, 0)
        self.conn_matrix = {}
        self.conn_headers =[]
        
        self.out_outline = (0, 0)
        self.blocks = {}
        self.channels = {}
        self.paths =[]
        
        self.fails = 0
        self.penalties = 0
        
        self.parse_csv()
        self.parse_cfg()

    def parse_csv(self):
        try:
            with open(self.csv_file, 'r', encoding='utf-8-sig') as f:
                reader = csv.reader(f)
                stage = ""
                for cols in reader:
                    if not cols or all(c.strip() == '' for c in cols): continue
                    
                    c0 = cols[0].strip()
                    if c0 == "BLOCK": stage = "BLOCK"; continue
                    elif c0 == "OUTLINE": stage = "OUTLINE"; continue
                    elif "α" in c0 or "alpha" in c0.lower(): 
                        self.alpha = float(cols[1].strip())
                        continue
                    elif "CONN MATRIX" in c0: stage = "CONN"; continue
                    
                    if stage == "BLOCK" and c0.startswith("BLK"):
                        name = c0
                        area = float(cols[1].strip()) if cols[1].strip() else 0.0
                        btype = cols[5].strip()
                        loc = cols[6].strip() if len(cols) > 6 else ""
                        ft_rates = [float(x.strip().strip('%'))/100 for x in cols[7:11]] if len(cols) >= 11 and cols[7].strip() else [0.2, 0.4, 0.8, 1.0]
                            
                        self.blocks_info[name] = {"area": area, "type": btype, "loc": loc, "ft_rates": ft_rates}
                    
                    elif stage == "OUTLINE" and c0 == "MAX":
                        self.max_outline = (float(cols[1].strip()), float(cols[2].strip()))
                        
                    elif stage == "CONN":
                        if c0 == "" and len(cols) > 1 and "BLK" in cols[1]:
                            self.conn_headers = [c.strip() for c in cols[1:] if c.strip()]
                        elif c0.startswith("BLK"):
                            name = c0
                            for i, val in enumerate(cols[1:len(self.conn_headers)+1]):
                                val = val.strip()
                                if val and float(val) > 0:
                                    target = self.conn_headers[i]
                                    # 建立無向連接
                                    pair = tuple(sorted([name, target]))
                                    self.conn_matrix[pair] = max(self.conn_matrix.get(pair, 0), float(val))
        except Exception as e:
            print(f"[FAIL] Format failed! 無法正確讀取 CSV: {e}")
            self.fails += 1

    def parse_cfg(self):
        try:
            with open(self.cfg_file, 'r', encoding='utf-8') as f:
                lines = [l.strip() for l in f.readlines() if l.strip()]
                
            stage = ""
            for line in lines:
                if line == "END": stage = ""; continue
                if line.startswith("OUTLINE"): stage = "OUTLINE"; continue
                elif line.startswith("BLOCK"): stage = "BLOCK"; continue
                elif line.startswith("CHANNEL"): stage = "CHANNEL"; continue
                elif line.startswith("PATH"): stage = "PATH"; 
                
                if stage == "OUTLINE" and not line.startswith("OUTLINE"):
                    cols = line.split()
                    self.out_outline = (float(cols[0]), float(cols[1]))
                elif stage == "BLOCK" and not line.startswith("BLOCK"):
                    cols = line.split()
                    self.blocks[cols[0]] = {"lx": float(cols[1]), "ly": float(cols[2]), "w": float(cols[3]), "h": float(cols[4])}
                elif stage == "CHANNEL" and not line.startswith("CHANNEL"):
                    cols = line.split()
                    self.channels[cols[0]] = {"lx": float(cols[1]), "ly": float(cols[2]), "w": float(cols[3]), "h": float(cols[4]), "nets_x": 0, "nets_y": 0}
                elif stage == "PATH" and not line == "PATH":
                    cols = line.split()
                    nets = float(cols[1])
                    
                    # 穩健解析 PATH 節點與 Edge (Src -> Node1 -> Node2 -> Tgt)
                    route_info =[]
                    idx = 2
                    src_name = cols[idx]
                    src_out = cols[idx+1]
                    route_info.append({"name": src_name, "in": None, "out": src_out})
                    idx += 2
                    
                    while idx < len(cols) - 2:
                        node_name1 = cols[idx]
                        edge_in = cols[idx+1]
                        node_name2 = cols[idx+2]
                        edge_out = cols[idx+3]
                        route_info.append({"name": node_name1, "in": edge_in, "out": edge_out})
                        idx += 4
                        
                    tgt_name = cols[idx]
                    tgt_in = cols[idx+1]
                    route_info.append({"name": tgt_name, "in": tgt_in, "out": None})
                    
                    self.paths.append({"nets": nets, "route_info": route_info})
        except Exception as e:
            print(f"[FAIL] Format failed! 無法正確解析 CFG: {e}")
            self.fails += 1

    def get_guiding_point(self, r1_name, e1, r2_name, e2):
        rect1 = self.blocks.get(r1_name, self.channels.get(r1_name))
        rect2 = self.blocks.get(r2_name, self.channels.get(r2_name))
        if not rect1 or not rect2: return (0, 0)
        
        if e1 == '3' and e2 == '1': 
            x = rect1['lx'] + rect1['w']
            y = (max(rect1['ly'], rect2['ly']) + min(rect1['ly'] + rect1['h'], rect2['ly'] + rect2['h'])) / 2.0
            return (x, y)
        elif e1 == '1' and e2 == '3': 
            x = rect1['lx']
            y = (max(rect1['ly'], rect2['ly']) + min(rect1['ly'] + rect1['h'], rect2['ly'] + rect2['h'])) / 2.0
            return (x, y)
        elif e1 == '2' and e2 == '4': 
            y = rect1['ly'] + rect1['h']
            x = (max(rect1['lx'], rect2['lx']) + min(rect1['lx'] + rect1['w'], rect2['lx'] + rect2['w'])) / 2.0
            return (x, y)
        elif e1 == '4' and e2 == '2': 
            y = rect1['ly']
            x = (max(rect1['lx'], rect2['lx']) + min(rect1['lx'] + rect1['w'], rect2['lx'] + rect2['w'])) / 2.0
            return (x, y)
        return (0, 0)

    def check_edge_location(self, name, b):
        loc_str = self.blocks_info[name]["loc"]
        if not loc_str or loc_str == "NONE": return True
        
        options = [o.strip() for o in loc_str.split(',')]
        for opt in options:
            match = True
            if 'T' in opt and abs(b['ly'] + b['h'] - self.out_outline[1]) > 1e-3: match = False
            if 'B' in opt and abs(b['ly'] - 0.0) > 1e-3: match = False
            if 'L' in opt and abs(b['lx'] - 0.0) > 1e-3: match = False
            if 'R' in opt and abs(b['lx'] + b['w'] - self.out_outline[0]) > 1e-3: match = False
            if match: return True
        return False

    def evaluate(self):
        print("========== ICCAD 2026 Problem E 評測 ==========")
        
        # 1. Check outline constraint violation
        if self.out_outline[0] > self.max_outline[0] + 1e-3 or self.out_outline[1] > self.max_outline[1] + 1e-3:
            print(f"[FAIL] Outline constraint violation! 配置超出了 MAX: {self.max_outline}")
            self.fails += 1
            
        for name, b in self.blocks.items():
            if b['lx'] < -1e-3 or b['ly'] < -1e-3 or b['lx']+b['w'] > self.out_outline[0]+1e-3 or b['ly']+b['h'] > self.out_outline[1]+1e-3:
                print(f"[FAIL] Outline constraint violation! Block {name} 跑出邊界。")
                self.fails += 1

        # 2. Check block overlap
        b_items = list(self.blocks.items())
        for i in range(len(b_items)):
            n1, b1 = b_items[i]
            for j in range(i+1, len(b_items)):
                n2, b2 = b_items[j]
                x_over = max(0, min(b1['lx']+b1['w'], b2['lx']+b2['w']) - max(b1['lx'], b2['lx']))
                y_over = max(0, min(b1['ly']+b1['h'], b2['ly']+b2['h']) - max(b1['ly'], b2['ly']))
                if x_over > 1e-3 and y_over > 1e-3:
                    print(f"[FAIL] Block overlap! {n1} 與 {n2} 發生重疊，重疊面積: {x_over*y_over:.2f}")
                    self.fails += 1

        # 3. Check edge block constraints
        for name, b in self.blocks.items():
            if self.blocks_info[name]["type"] == "EDGE":
                if not self.check_edge_location(name, b):
                    print(f"[FAIL] Edge constraint violation! {name} 未能放置在要求的邊界 {self.blocks_info[name]['loc']}")
                    self.fails += 1

        # 4. Parse PATH, check routing open, and compute exact HPWL
        total_hpwl = 0
        ft_block_nets = {b: 0 for b in self.blocks if self.blocks_info[b]["type"] == "SOFT"}
        routed_matrix = {}
        
        for p in self.paths:
            nets = p["nets"]
            r_info = p["route_info"]
            
            src = r_info[0]["name"]
            tgt = r_info[-1]["name"]
            pair = tuple(sorted([src, tgt]))
            routed_matrix[pair] = routed_matrix.get(pair, 0) + nets
            
            guiding_points =[]
            for i in range(len(r_info) - 1):
                r1 = r_info[i]["name"]
                e1 = r_info[i]["out"]
                r2 = r_info[i+1]["name"]
                e2 = r_info[i+1]["in"]
                
                # Compute guiding point
                gp = self.get_guiding_point(r1, e1, r2, e2)
                guiding_points.append(gp)
                
                # If the intermediate node is a channel, determine and accumulate directional flow
                if i > 0 and r1.startswith("CH"):
                    in_dir = r_info[i]["in"]
                    out_dir = r_info[i]["out"]
                    has_x = in_dir in ['1', '3'] or out_dir in ['1', '3']
                    has_y = in_dir in ['2', '4'] or out_dir in ['2', '4']
                    if has_x: self.channels[r1]["nets_x"] += nets
                    if has_y: self.channels[r1]["nets_y"] += nets
                    
                # If the intermediate node is a soft block, accumulate FT
                if i > 0 and r1.startswith("BLK") and self.blocks_info[r1]["type"] == "SOFT":
                    ft_block_nets[r1] += nets
            
            # Compute the exact total Manhattan distance across segments
            length = 0
            for k in range(len(guiding_points) - 1):
                pt1 = guiding_points[k]
                pt2 = guiding_points[k+1]
                length += abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1])
            total_hpwl += length * nets

        # 5. Check routing open
        for pair, req_nets in self.conn_matrix.items():
            act_nets = routed_matrix.get(pair, 0)
            if act_nets < req_nets:
                print(f"[FAIL] Routing open! Net {pair[0]}-{pair[1]} 需求 {req_nets}，實際只繞了 {act_nets}")
                self.fails += 1

        # 6. Check channel overflow (Penalty)
        for ch, data in self.channels.items():
            cap_x = data["h"] * 25.0  # 水平穿越看高度
            cap_y = data["w"] * 25.0  # 垂直穿越看寬度
            if data["nets_x"] > cap_x + 1e-3:
                print(f"[PENALTY] Channel overflow! {ch} 水平流量 {data['nets_x']} > 容量 {cap_x:.1f}")
                self.penalties += 1
            if data["nets_y"] > cap_y + 1e-3:
                print(f"[PENALTY] Channel overflow! {ch} 垂直流量 {data['nets_y']} > 容量 {cap_y:.1f}")
                self.penalties += 1

        # 7. Check feedthrough overflow (Penalty)
        for blk, ftnets in ft_block_nets.items():
            rates = self.blocks_info[blk]["ft_rates"]
            if ftnets <= 3000: rate = rates[0]
            elif ftnets <= 6000: rate = rates[1]
            elif ftnets <= 9000: rate = rates[2]
            else: rate = rates[3]
            
            base_side = math.sqrt(self.blocks_info[blk]["area"])
            expected_side = base_side + (ftnets / 25.0) * rate / 2.0
            expected_area = expected_side ** 2
            actual_area = self.blocks[blk]["w"] * self.blocks[blk]["h"]
            
            if actual_area < expected_area - 1e-3:
                print(f"[PENALTY] Feedthrough overflow! Soft Block {blk} 面積不足。需要: {expected_area:.2f}, 實際: {actual_area:.2f}")
                self.penalties += 1

        cost = (self.out_outline[0] * self.out_outline[1]) + self.alpha * total_hpwl
        print("\n=== 評測總結 ===")
        print(f"Outline (WxH) : {self.out_outline[0]:.2f} x {self.out_outline[1]:.2f}")
        print(f"Area Cost     : {self.out_outline[0] * self.out_outline[1]:.2f}")
        print(f"HPWL (精確)   : {total_hpwl:.2f}")
        print(f"Total Cost    : {cost:.2f} (alpha={self.alpha})")
        print("-------------------")
        print(f"Total Fails   : {self.fails}")
        print(f"Total Penalties: {self.penalties}")
        if self.fails == 0:
            print(">> 恭喜！無任何 FAIL 違規。")
        else:
            print(">> 注意！存在 FAIL 違規，此解答將不予計分。")

    def plot(self):
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlim(-100, self.max_outline[0] * 1.1)
        ax.set_ylim(-100, self.max_outline[1] * 1.1)
        
        # 繪製 MAX Outline
        ax.add_patch(patches.Rectangle((0, 0), self.max_outline[0], self.max_outline[1], 
                                       fill=False, edgecolor='red', linewidth=2, linestyle='-.', label='Max Outline'))
        # 繪製 Actual Outline
        ax.add_patch(patches.Rectangle((0, 0), self.out_outline[0], self.out_outline[1], 
                                       fill=False, edgecolor='black', linewidth=2, linestyle='--', label='Actual Outline'))
        
        for name, b in self.blocks.items():
            btype = self.blocks_info[name]["type"]
            color = 'lightcoral' if btype in ['EDGE', 'MACRO'] else 'lightblue'
            ax.add_patch(patches.Rectangle((b['lx'], b['ly']), b['w'], b['h'], facecolor=color, edgecolor='black', alpha=0.7))
            ax.text(b['lx'] + b['w']/2, b['ly'] + b['h']/2, name, ha='center', va='center', fontsize=9, weight='bold')

        for name, c in self.channels.items():
            ax.add_patch(patches.Rectangle((c['lx'], c['ly']), c['w'], c['h'], facecolor='lightgreen', edgecolor='green', alpha=0.3, linestyle=':'))
            ax.text(c['lx'] + c['w']/2, c['ly'] + c['h']/2, name, ha='center', va='center', fontsize=7, color='darkgreen')

        for p in self.paths:
            r_info = p["route_info"]
            gps =[]
            for i in range(len(r_info)-1):
                gp = self.get_guiding_point(r_info[i]["name"], r_info[i]["out"], r_info[i+1]["name"], r_info[i+1]["in"])
                if gp != (0,0): gps.append(gp)
            if gps:
                xs, ys = zip(*gps)
                ax.plot(xs, ys, marker='o', markersize=4, linestyle='-', linewidth=1.5, alpha=0.8)

        plt.title("ICCAD 2026 Early Floorplanning Visualization")
        plt.xlabel("X (um)")
        plt.ylabel("Y (um)")
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.show()
        # save the plot to a file
        plt.savefig('floorplan.png')
        print(f"Plot saved to floorplan.png")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python evaluator.py <input.csv> <output.cfg>")
        sys.exit(1)
        
    evaluator = Evaluator(sys.argv[1], sys.argv[2])
    evaluator.evaluate()
    evaluator.plot()