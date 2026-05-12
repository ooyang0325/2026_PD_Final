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
        
        self.parse_csv()
        self.parse_cfg()

    def parse_csv(self):
        # 使用 Python 內建的 csv 模組，自動處理引號內的逗號 "TL,LT"
        with open(self.csv_file, 'r', encoding='utf-8-sig') as f:
            reader = csv.reader(f)
            stage = ""
            for cols in reader:
                if not cols: continue
                # 略過全空的行 (例如只有一堆逗號 ,,,,,,)
                if all(c.strip() == '' for c in cols): continue
                
                c0 = cols[0].strip()
                if c0 == "BLOCK": 
                    stage = "BLOCK"
                    continue
                elif c0 == "OUTLINE": 
                    stage = "OUTLINE"
                    continue
                elif "α" in c0 or "alpha" in c0.lower(): 
                    self.alpha = float(cols[1].strip())
                    continue
                elif "CONN MATRIX" in c0: 
                    stage = "CONN"
                    continue
                
                # 讀取 BLOCK
                if stage == "BLOCK" and c0.startswith("BLK"):
                    name = c0
                    area = float(cols[1].strip()) if cols[1].strip() else 0.0
                    ar = cols[4].strip()
                    btype = cols[5].strip()
                    loc = cols[6].strip() if len(cols) > 6 else ""
                    # 處理 FT Rates，若為空則給預設值
                    if len(cols) >= 11 and cols[7].strip():
                        ft_rates = [float(x.strip().strip('%'))/100 for x in cols[7:11]]
                    else:
                        ft_rates = [0.2, 0.4, 0.8, 1.0]
                        
                    self.blocks_info[name] = {
                        "area": area, "ar": ar, "type": btype, 
                        "loc": loc, "ft_rates": ft_rates
                    }
                
                # 讀取 OUTLINE
                elif stage == "OUTLINE" and c0 == "MAX":
                    self.max_outline = (float(cols[1].strip()), float(cols[2].strip()))
                    
                # 讀取 CONN MATRIX
                elif stage == "CONN":
                    if c0 == "" and len(cols) > 1 and "BLK" in cols[1]:
                        # 紀錄 Header 對應的目標名稱
                        self.conn_headers = [c.strip() for c in cols[1:] if c.strip()]
                    elif c0.startswith("BLK"):
                        name = c0
                        # 依據 Header 長度映射 Net 數量
                        for i, val in enumerate(cols[1:len(self.conn_headers)+1]):
                            val = val.strip()
                            if val and float(val) > 0:
                                target = self.conn_headers[i]
                                self.conn_matrix[(name, target)] = float(val)

    def parse_cfg(self):
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
                self.channels[cols[0]] = {"lx": float(cols[1]), "ly": float(cols[2]), "w": float(cols[3]), "h": float(cols[4]), "nets": 0}
            elif stage == "PATH" and not line == "PATH":
                cols = line.split()
                nets = float(cols[1])
                route = cols[2:]
                self.paths.append({"nets": nets, "route": route})

    def get_guiding_point(self, r1, e1, r2, e2):
        # 計算相鄰兩個區塊的交界線段中心點 (Guiding Point)
        rect1 = self.blocks.get(r1, self.channels.get(r1))
        rect2 = self.blocks.get(r2, self.channels.get(r2))
        
        if not rect1 or not rect2: return (0, 0)
        
        if e1 == '3' and e2 == '1': # R1 Right, R2 Left
            x = rect1['lx'] + rect1['w']
            y_min = max(rect1['ly'], rect2['ly'])
            y_max = min(rect1['ly'] + rect1['h'], rect2['ly'] + rect2['h'])
            return (x, (y_min + y_max) / 2.0)
        elif e1 == '1' and e2 == '3': # R1 Left, R2 Right
            x = rect1['lx']
            y_min = max(rect1['ly'], rect2['ly'])
            y_max = min(rect1['ly'] + rect1['h'], rect2['ly'] + rect2['h'])
            return (x, (y_min + y_max) / 2.0)
        elif e1 == '2' and e2 == '4': # R1 Top, R2 Bottom
            y = rect1['ly'] + rect1['h']
            x_min = max(rect1['lx'], rect2['lx'])
            x_max = min(rect1['lx'] + rect1['w'], rect2['lx'] + rect2['w'])
            return ((x_min + x_max) / 2.0, y)
        elif e1 == '4' and e2 == '2': # R1 Bottom, R2 Top
            y = rect1['ly']
            x_min = max(rect1['lx'], rect2['lx'])
            x_max = min(rect1['lx'] + rect1['w'], rect2['lx'] + rect2['w'])
            return ((x_min + x_max) / 2.0, y)
        return (0, 0) 

    def evaluate(self):
        print("=== 評測開始 ===")
        errors = 0
        
        # 1. 檢查 Outline 限制
        if self.out_outline[0] > self.max_outline[0] or self.out_outline[1] > self.max_outline[1]:
            print(f"[FAIL] Outline 超出限制! 限制: {self.max_outline}, 實際: {self.out_outline}")
            errors += 1
            
        # 2. 計算精確 HPWL (利用 Guiding Points)
        total_hpwl = 0
        ft_block_nets = {b: 0 for b in self.blocks if self.blocks_info[b]["type"] == "SOFT"}
        
        for p in self.paths:
            nets = p["nets"]
            route = p["route"]
            guiding_points =[]
            
            # 遍歷 PATH 解析 (Block/CH, edge_in, Block/CH, edge_out...)
            for i in range(0, len(route)-2, 2):
                r1, e1 = route[i], route[i+1]
                r2, e2 = route[i+2], route[i+3]
                if r1 == r2: continue # 同一個 block 內部穿越
                
                gp = self.get_guiding_point(r1, e1, r2, e2)
                guiding_points.append(gp)
                
                # 累計 Channel / Soft Block 負載
                if r1.startswith("CH"): self.channels[r1]["nets"] += nets
                if r1.startswith("BLK") and self.blocks_info[r1]["type"] == "SOFT":
                    if i > 0: ft_block_nets[r1] += nets # 中間穿越
            
            if guiding_points:
                min_x = min(pt[0] for pt in guiding_points)
                max_x = max(pt[0] for pt in guiding_points)
                min_y = min(pt[1] for pt in guiding_points)
                max_y = max(pt[1] for pt in guiding_points)
                hpwl = (max_x - min_x) + (max_y - min_y)
                total_hpwl += hpwl * nets

        cost = (self.out_outline[0] * self.out_outline[1]) + self.alpha * total_hpwl
        print(f"[INFO] 總面積: {self.out_outline[0] * self.out_outline[1]:.2f}")
        print(f"[INFO] 總 HPWL: {total_hpwl:.2f}")
        print(f"[INFO] 最終 COST (Area + alpha*HPWL): {cost:.2f}")

        # 3. 檢查 Channel 擁擠度 (Capacity)
        for ch, data in self.channels.items():
            cap = max(data["w"], data["h"]) * 25.0 # 取長邊作為跨越方向的寬度
            if data["nets"] > cap:
                print(f"[WARNING] 通道 {ch} Overflow! 容量: {cap:.1f}, 實際: {data['nets']}")
                errors += 1

        # 4. 檢查 Soft Block 面積膨脹合法性
        for blk, ftnets in ft_block_nets.items():
            rates = self.blocks_info[blk]["ft_rates"]
            # 依據題意區間：<=3000, 3000~6000, 6000~9000, >9000
            if ftnets <= 3000: rate = rates[0]
            elif ftnets <= 6000: rate = rates[1]
            elif ftnets <= 9000: rate = rates[2]
            else: rate = rates[3]
            
            base_side = math.sqrt(self.blocks_info[blk]["area"])
            expected_side = base_side + (ftnets / 25.0) * rate / 2.0
            expected_area = expected_side ** 2
            
            actual_area = self.blocks[blk]["w"] * self.blocks[blk]["h"]
            if actual_area < expected_area - 1.0: # 容許 1.0 的浮點誤差
                print(f"[FAIL] Soft Block {blk} 面積不足以容納 FT! 需要: {expected_area:.2f}, 實際: {actual_area:.2f}")
                errors += 1

        print(f"=== 評測結束，共發現 {errors} 個潛在違規 ===")

    def plot(self):
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlim(0, self.max_outline[0] * 1.1)
        ax.set_ylim(0, self.max_outline[1] * 1.1)
        
        # 繪製 Outline
        ax.add_patch(patches.Rectangle((0, 0), self.out_outline[0], self.out_outline[1], 
                                       fill=False, edgecolor='black', linewidth=2, linestyle='--'))
        
        # 繪製 Blocks
        for name, b in self.blocks.items():
            btype = self.blocks_info[name]["type"]
            color = 'lightcoral' if btype in ['EDGE', 'MACRO'] else 'lightblue'
            ax.add_patch(patches.Rectangle((b['lx'], b['ly']), b['w'], b['h'], 
                                           facecolor=color, edgecolor='black', alpha=0.7))
            ax.text(b['lx'] + b['w']/2, b['ly'] + b['h']/2, name, ha='center', va='center', fontsize=9, weight='bold')

        # 繪製 Channels
        for name, c in self.channels.items():
            ax.add_patch(patches.Rectangle((c['lx'], c['ly']), c['w'], c['h'], 
                                           facecolor='lightgreen', edgecolor='green', alpha=0.3, linestyle=':'))
            ax.text(c['lx'] + c['w']/2, c['ly'] + c['h']/2, name, ha='center', va='center', fontsize=7, color='darkgreen')

        # 繪製 Paths (飛線連線)
        for p in self.paths:
            route = p["route"]
            gps =[]
            for i in range(0, len(route)-2, 2):
                if route[i] != route[i+2]:
                    gps.append(self.get_guiding_point(route[i], route[i+1], route[i+2], route[i+3]))
            if gps:
                xs, ys = zip(*gps)
                ax.plot(xs, ys, marker='o', markersize=4, linestyle='-', linewidth=1.5, alpha=0.8)

        plt.title("Early Floorplanning with Global Route Visualization")
        plt.xlabel("X (um)")
        plt.ylabel("Y (um)")
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python evaluator.py <input.csv> <output.cfg>")
        sys.exit(1)
        
    evaluator = Evaluator(sys.argv[1], sys.argv[2])
    evaluator.evaluate()
    evaluator.plot()