import csv
import random
import math
import argparse

class ICCADTestcaseGenerator:
    def __init__(self, num_blocks, utilization=0.6, alpha=0.5, seed=None):
        self.num_blocks = num_blocks
        self.utilization = utilization
        self.alpha = alpha
        if seed is not None:
            random.seed(seed)
            
        self.blocks =[]
        self.outline_w = 0
        self.outline_h = 0
        self.conn_matrix =[]

    def generate(self):
        self._generate_blocks()
        self._calculate_outline()
        self._generate_connections()

    def _generate_blocks(self):
        # 決定區塊類型的比例：大約 20% EDGE, 30% MACRO, 50% SOFT
        num_edge = max(1, int(self.num_blocks * 0.2))
        num_macro = max(1, int(self.num_blocks * 0.3))
        num_soft = self.num_blocks - num_edge - num_macro

        types = ['EDGE'] * num_edge + ['MACRO'] * num_macro + ['SOFT'] * num_soft
        random.shuffle(types)

        # Edge Block 可用的邊界限制 (參考題目圖 4 的 12 個區域)
        edge_locations = ["TL,LT", "TM", "TR,RT", "LM", "RM", "BL,LB", "BM", "BR,RB", "T", "B", "L", "R"]

        for i, btype in enumerate(types):
            name = f"BLK{i+1:02d}"
            
            # 隨機生成基礎面積 (100,000 ~ 1,000,000)
            area = round(random.uniform(100000.0, 1000000.0), 1)
            
            if btype in ['EDGE', 'MACRO']:
                # 長寬比盡量接近 1:1 ~ 1:2
                ar = random.uniform(0.5, 2.0)
                width = round(math.sqrt(area * ar), 1)
                height = round(area / width, 1)
                # 重新校準面積確保 W * H = Area
                area = round(width * height, 1)
                
                loc = random.choice(edge_locations) if btype == 'EDGE' else ""
                ar_range = "1"
            else: # SOFT
                width = ""
                height = ""
                loc = ""
                ar_range = "0.5,2"

            self.blocks.append({
                "name": name,
                "area": area,
                "w": width,
                "h": height,
                "ar_range": ar_range,
                "type": btype,
                "loc": loc,
                "ft_rates": ["20%", "40%", "80%", "100%"]
            })

    def _calculate_outline(self):
        total_area = sum(b["area"] for b in self.blocks)
        target_outline_area = total_area / self.utilization
        
        # 隨機產生一個 Outline 的長寬比 (0.8 ~ 1.25)
        outline_ar = random.uniform(0.8, 1.25)
        self.outline_w = round(math.sqrt(target_outline_area * outline_ar), 2)
        self.outline_h = round(target_outline_area / self.outline_w, 2)

    def _generate_connections(self):
        n = self.num_blocks
        self.conn_matrix = [[0] * n for _ in range(n)]
        
        # 連線密度：約 30% 的 Block 之間有連線
        sparsity = 0.3
        
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < sparsity:
                    # 隨機生成 Net 數量 (100 ~ 3000)
                    # 刻意產生一些容易觸發 FT Overflow (> 1000) 的測資
                    nets = random.choice([100, 200, 300, 500, 800, 1200, 1500, 2000])
                    self.conn_matrix[i][j] = nets
                    self.conn_matrix[j][i] = nets

    def export_csv(self, filename):
        with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            
            # --- 1. BLOCK Section ---
            writer.writerow(["BLOCK", "AREA", "WIDTH", "HEIGHT", "ASPECT RATIO RANGE", "EDGE\nHARD MACRO\nSOFT", "LOCATION", "FT CONVERSION", "", "", ""])
            writer.writerow(["", "", "", "", "", "", "", "<=1000", ">3000, <=6000", ">6000, <= 9000", ">9000"])
            
            for b in self.blocks:
                row = [
                    b["name"], b["area"], b["w"], b["h"], b["ar_range"], 
                    b["type"], b["loc"]
                ] + b["ft_rates"]
                writer.writerow(row)
                
            writer.writerow([""] * 11)
            writer.writerow([""] * 11)
            
            # --- 2. OUTLINE Section ---
            writer.writerow(["OUTLINE", "WIDTH", "HEIGHT"] + [""] * 8)
            writer.writerow(["MAX", self.outline_w, self.outline_h] + [""] * 8)
            writer.writerow([""] * 11)
            
            # --- 3. ALPHA Section ---
            writer.writerow(["α", self.alpha] + [""] * 9)
            writer.writerow([""] * 11)
            
            # --- 4. CONN MATRIX Section ---
            writer.writerow(["CONN MATRIX"] + [""] * 10)
            writer.writerow(["ON CHIP INTERFACE CONNECTION (1-1) MATRIX"] + [""] * 10)
            
            header = [""] + [b["name"] for b in self.blocks]
            writer.writerow(header)
            
            for i in range(self.num_blocks):
                row = [self.blocks[i]["name"]] + self.conn_matrix[i]
                writer.writerow(row)

        print(f"[Success] 已成功生成測資: {filename}")
        print(f"  - Block 數量: {self.num_blocks}")
        print(f"  - Outline: {self.outline_w} x {self.outline_h} (Utilization: {self.utilization*100:.1f}%)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ICCAD 2026 Problem E Testcase Generator")
    parser.add_argument("-n", "--num_blocks", type=int, default=15, help="Number of blocks (Default: 15)")
    parser.add_argument("-u", "--utilization", type=float, default=0.6, help="Target utilization ratio (Default: 0.6)")
    parser.add_argument("-a", "--alpha", type=float, default=0.5, help="Alpha weight for HPWL (Default: 0.5)")
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument("-o", "--output", type=str, default="testcase_auto.csv", help="Output CSV filename")
    
    args = parser.parse_args()
    
    if args.num_blocks > 50:
        print("[Warning] 大會規定 Block 數量 < 50，您輸入的數量偏高，可能會增加 SA 難度！")
        
    generator = ICCADTestcaseGenerator(
        num_blocks=args.num_blocks,
        utilization=args.utilization,
        alpha=args.alpha,
        seed=args.seed
    )
    
    generator.generate()
    generator.export_csv(args.output)