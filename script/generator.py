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

    # 12 個合法邊界區域 (參考題目圖 4)。注意：每個「角落」(TL/TR/BL/BR) 在幾何上
    # 只能容納「一個」Block (兩個 Block 不可能同時佔據同一個角落而不重疊)，因此每個
    # 區域至多指派給一個 Edge Block，絕不重複，避免產生無解的測資。
    EDGE_LOCATIONS = ["TL,LT", "TR,RT", "BL,LB", "BR,RB",   # 4 個角落 (各限一個)
                      "TM", "BM", "LM", "RM",               # 4 個邊中點
                      "T", "B", "L", "R"]                   # 4 個自由邊

    @staticmethod
    def _edge_axes(loc):
        """取出 location 所約束的邊 (T/B/L/R)。評測器只看 T/B/L/R，M 不予理會。"""
        if not loc:
            return (False, False, False, False)
        opt = loc.split(',')[0]  # 兩個選項等價，取第一個即可
        return ('T' in opt, 'B' in opt, 'L' in opt, 'R' in opt)

    def _generate_blocks(self):
        # 決定區塊類型的比例：大約 20% EDGE, 30% MACRO, 50% SOFT
        num_edge = max(1, int(self.num_blocks * 0.2))
        num_macro = max(1, int(self.num_blocks * 0.3))

        # Edge Block 數量不可超過可用的「不重複」合法區域數，否則必然有兩個 Block
        # 被指派到同一個角落 -> 幾何無解。超出時將多餘的 Edge 改為 MACRO。
        if num_edge > len(self.EDGE_LOCATIONS):
            num_macro += num_edge - len(self.EDGE_LOCATIONS)
            num_edge = len(self.EDGE_LOCATIONS)
        num_soft = self.num_blocks - num_edge - num_macro

        types = ['EDGE'] * num_edge + ['MACRO'] * num_macro + ['SOFT'] * num_soft
        random.shuffle(types)

        # 隨機挑選「不重複」的邊界區域給 Edge Block（每個角落至多一個）。
        available = self.EDGE_LOCATIONS.copy()
        random.shuffle(available)

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

                # 每個 Edge Block 取一個唯一區域（pop 確保不重複）
                loc = available.pop() if btype == 'EDGE' else ""
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

        # ── 邊界可行性保證 ────────────────────────────────────────────────────
        # 即使每個角落唯一，若 Outline 太小，同一條邊上的 Edge Block 仍可能放不下
        # （沿邊並排的總長度超過邊長）而必然重疊/違規。這裡確保每條邊都裝得下其
        # 上的 Edge Block：上下邊看寬度和、左右邊看高度和；不足則放大 Outline。
        top = bot = left = right = 0.0
        for b in self.blocks:
            if b["type"] != "EDGE":
                continue
            has_t, has_b, has_l, has_r = self._edge_axes(b["loc"])
            w, h = float(b["w"]), float(b["h"])
            if has_t: top += w
            if has_b: bot += w
            if has_l: left += h
            if has_r: right += h
        # 留 1% 餘裕避免邊界剛好相切
        self.outline_w = round(max(self.outline_w, top * 1.01, bot * 1.01), 2)
        self.outline_h = round(max(self.outline_h, left * 1.01, right * 1.01), 2)

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