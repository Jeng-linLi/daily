"""掃描線 + 線段樹：求 n 個軸對齊矩形的聯集面積（Union Area）

題意：
    給定 n 個**軸對齊**（axis-aligned）的矩形，第 i 個矩形的左下角與右上角
    分別是 (x1[i], y1[i]) 與 (x2[i], y2[i])，其中 x1[i] < x2[i]、y1[i] < y2[i]。
    求這 n 個矩形**聯集**（union）的面積，也就是把重疊部分只算一次的總面積。
    限制：1 <= n <= 100000，|座標| <= 10^9。

思路：
    ### 1. 為什麼不能直接加總
    矩形彼此可能重疊，直接把面積相加會重複計算。
    樸素做法是「格點化」：把所有 x、y 座標收集起來排序，
    把平面切成 O(n^2) 個小格子再逐格判斷 —— 格子數太多，O(n^2) 不可行。

    ### 2. 掃描線（sweep line）
    讓一條垂直線從左往右掃過 x 軸。
    對每個矩形 (x1, y1, x2, y2) 產生兩個事件：
        (x1, +1, [y1, y2))  —— 進入矩形，把 y 區間 [y1, y2) 的覆蓋次數 +1
        (x2, -1, [y1, y2))  —— 離開矩形，把 y 區間 [y1, y2) 的覆蓋次數 -1
    相鄰兩個事件的 x 之間，覆蓋情況完全不變 —— 這就是一個「條帶」（slab）。
    若條帶寬為 dx、此時 y 方向的聯集長度為 L，則該條帶的面積 = dx * L。
    把所有條帶的面積加起來就是答案。
    所以問題變成：**動態維護一組 y 區間的聯集長度**。

    ### 3. 座標壓縮 + 線段樹
    把所有出現過的 y 座標排序去重，得到 ys[0..m-1]。
    線段樹的葉子不是「點」而是「基本區間」[ys[i], ys[i+1])，共 m-1 個。
    節點維護兩個量：
        cnt   —— 完整覆蓋此區間的「矩形數」（不下傳的覆蓋次數）
        length—— 此區間內目前被覆蓋的總長度
    拉格朗日式的更新規則（經典的「覆蓋次數」寫法）：
        若 cnt > 0：length = 區間的實際長度（整段被蓋住）
        否則若為葉子：length = 0
        否則：length = 左子 + 右子
    根節點的 length 就是當前整個 y 軸的聯集長度。

    ### 4. 為什麼不需要 push-down
    這裡刻意不做 lazy 下傳：cnt 記錄的是「有多少個矩形完整覆蓋了這個節點」，
    而 deeper 的 cnt 記錄的是「在父節點未被完整覆蓋的前提下，
    額外覆蓋此區間的矩形數」。當 cnt 由正變 0 時，
    子節點的 length 早已包含了所有下層更新，直接相加即可。
    這是矩形聯集面積的標準寫法。

    ### 5. 精度與整數性
    所有座標都是整數，條帶寬度與聯集長度也都是整數，
    所以面積一定是整數 —— 直接用整數運算即可，不需要浮點數。

複雜度：
    座標壓縮 + 排序    O(n log n)
    2n 次線段樹更新    O(n log n)
    總計              O(n log n)，空間 O(n)

應用場景（Applications）：
    - **地圖 / GIS 覆蓋分析**：多個行政區、地塊、保護區的聯集面積，
      是選址分析（site selection）的基本運算
    - **資料中心機櫃佔地**：把機櫃視為矩形，求機房內被機櫃佔用的總面積
    - **建築物輪廓**：多棟建物平面圖的聯集，估算基地使用率（容積率、覆蓋率）
    - **影像處理 / UI 渲染**：多個有透明重疊的圖層合成後的實際覆蓋像素數
    - **碰撞與覆蓋檢測**：機械臂多臂在平面上的作業區域聯集
    - **廣告位出租**：一塊版面被多個廣告區塊佔用時的可售面積
"""

import sys
from typing import List, Optional, Sequence, Tuple

NMAX = 100000
COORD_MAX = 10 ** 9


class SegmentTree:
    """y 軸上的「覆蓋次數 + 聯集長度」線段樹。

    葉子 i 代表基本區間 [ys[i], ys[i+1])，共 n = len(ys) - 1 個。
    """

    __slots__ = ("ys", "n", "cnt", "length")

    def __init__(self, ys: Sequence[int]) -> None:
        self.ys = ys
        self.n = len(ys) - 1
        size = 4 * (self.n + 2)
        self.cnt = [0] * size
        self.length = [0] * size

    def update(self, ql: int, qr: int, delta: int) -> None:
        """把基本區間 [ql, qr) 的覆蓋次數加上 delta。"""
        if self.n <= 0 or ql >= qr:
            return
        ys = self.ys
        cnt = self.cnt
        length = self.length

        def rec(node: int, l: int, r: int) -> None:
            if ql <= l and r <= qr:
                cnt[node] += delta
            else:
                mid = (l + r) >> 1
                if ql < mid:
                    rec(node << 1, l, mid)
                if qr > mid:
                    rec((node << 1) | 1, mid, r)
            if cnt[node] > 0:
                length[node] = ys[r] - ys[l]
            elif r - l == 1:
                length[node] = 0
            else:
                length[node] = length[node << 1] + length[(node << 1) | 1]

        rec(1, 0, self.n)

    def covered_length(self) -> int:
        """目前整個 y 軸的聯集長度。"""
        return self.length[1]


def union_area(rects: Sequence[Tuple[int, int, int, int]]) -> int:
    """掃描線 + 線段樹求聯集面積。"""
    if not rects:
        return 0
    ys_set = set()
    events: List[Tuple[int, int, int, int]] = []
    for x1, y1, x2, y2 in rects:
        ys_set.add(y1)
        ys_set.add(y2)
        events.append((x1, 1, y1, y2))
        events.append((x2, -1, y1, y2))
    ys = sorted(ys_set)
    y_index = {y: i for i, y in enumerate(ys)}
    events.sort()

    tree = SegmentTree(ys)
    area = 0
    prev_x = events[0][0]
    i = 0
    total = len(events)
    while i < total:
        x = events[i][0]
        # 條帶 [prev_x, x) 的面積
        area += (x - prev_x) * tree.covered_length()
        # 同一個 x 的所有事件一起處理
        while i < total and events[i][0] == x:
            _, delta, y1, y2 = events[i]
            tree.update(y_index[y1], y_index[y2], delta)
            i += 1
        prev_x = x
    return area


def brute_area(rects: Sequence[Tuple[int, int, int, int]]) -> int:
    """參考實作：座標壓縮後逐格累加（O(n^2)，只給小測試用）。"""
    if not rects:
        return 0
    xs = sorted({r[0] for r in rects} | {r[2] for r in rects})
    ys = sorted({r[1] for r in rects} | {r[3] for r in rects})
    area = 0
    for i in range(len(xs) - 1):
        # 取區間內的一個整數代表點（下取整）
        xm = xs[i] + (xs[i + 1] - xs[i]) // 2
        for j in range(len(ys) - 1):
            ym = ys[j] + (ys[j + 1] - ys[j]) // 2
            for x1, y1, x2, y2 in rects:
                if x1 <= xm < x2 and y1 <= ym < y2:
                    area += (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j])
                    break
    return area


def run_tests() -> None:
    # ---- 固定用例 ----
    # 1) 單一矩形
    assert union_area([(0, 0, 4, 3)]) == 12
    # 2) 兩個不重疊的矩形
    assert union_area([(0, 0, 2, 2), (5, 0, 7, 2)]) == 4 + 4
    # 3) 兩個部分重疊的矩形：[0,2)x[0,2] 與 [1,3)x[0,2]
    #    聯集 = [0,3) x [0,2] = 3*2 = 6
    assert union_area([(0, 0, 2, 2), (1, 0, 3, 2)]) == 6
    # 4) 完全包含：外面 [0,4)x[0,4]，裡面 [1,3)x[1,3]
    assert union_area([(0, 0, 4, 4), (1, 1, 3, 3)]) == 16
    # 5) 完全相同的兩個矩形（重疊只算一次）
    assert union_area([(0, 0, 2, 2), (0, 0, 2, 2)]) == 4
    # 6) 十字形：[0,4)x[1,3] 與 [1,3)x[0,4]
    #    聯集 = 4*2 + 2*4 - 2*2 = 8 + 8 - 4 = 12
    assert union_area([(0, 1, 4, 3), (1, 0, 3, 4)]) == 12
    # 7) 負座標
    assert union_area([(-2, -2, 2, 2)]) == 16
    assert union_area([(-2, -2, 0, 0), (0, 0, 2, 2)]) == 8
    # 8) 共享邊界不算重疊
    assert union_area([(0, 0, 2, 2), (2, 0, 4, 2)]) == 8
    # 9) 三個矩形錯開
    assert union_area([(0, 0, 2, 2), (1, 1, 3, 3), (2, 2, 4, 4)]) == 4 + 4 + 4 - 1 - 1 - 1 + 1

    # ---- 與逐格參考實作對拍 ----
    import random

    rng = random.Random(20261030)
    for _ in range(200):
        n = rng.randint(1, 6)
        rects = []
        for _ in range(n):
            x1 = rng.randint(-8, 8)
            x2 = x1 + rng.randint(1, 6)
            y1 = rng.randint(-8, 8)
            y2 = y1 + rng.randint(1, 6)
            rects.append((x1, y1, x2, y2))
        assert union_area(rects) == brute_area(rects), rects

    # 稍大規模（只跟暴力比，成本可接受）
    for _ in range(20):
        n = rng.randint(5, 12)
        rects = []
        for _ in range(n):
            x1 = rng.randint(-20, 20)
            x2 = x1 + rng.randint(1, 10)
            y1 = rng.randint(-20, 20)
            y2 = y1 + rng.randint(1, 10)
            rects.append((x1, y1, x2, y2))
        assert union_area(rects) == brute_area(rects), rects

    # 純網格（很多矩形對齊到同一組座標）：答案就是外框面積
    rects = [(i, j, i + 1, j + 1) for i in range(10) for j in range(10)]
    assert union_area(rects) == 100

    # 一個大矩形 + 許多小矩形：面積必須等於大矩形面積
    rects = [(0, 0, 100, 100)] + [(i * 3, i * 3, i * 3 + 2, i * 3 + 2) for i in range(30)]
    assert union_area(rects) == 100 * 100


# ---------------------------------------------------------------- IO 模式

def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；非純整數或超出 ±10^18 回 None。"""
    if not tok:
        return None
    i = 0
    if tok[0] in "+-":
        i = 1
    if i == len(tok):
        return None
    for ch in tok[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(tok)
    if val > 10 ** 18:
        return None
    return val


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0
    out: List[str] = []

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    n = _parse_int(nxt() or "")
    if n is None or not (1 <= n <= NMAX):
        bad()
        return
    rects: List[Tuple[int, int, int, int]] = []
    for _ in range(n):
        vals = []
        for _ in range(4):
            v = _parse_int(nxt() or "")
            if v is None or abs(v) > COORD_MAX:
                bad()
                return
            vals.append(v)
        x1, y1, x2, y2 = vals
        if x1 >= x2 or y1 >= y2:
            bad()
            return
        rects.append((x1, y1, x2, y2))
    area = union_area(rects)
    out.append(str(area))
    for line in out:
        emit(line + "\n")


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
