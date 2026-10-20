"""K-D Tree：k 維空間的二分搜尋樹，把「最近鄰」從 O(n) 降到期望 O(log n)

題意：
    給平面上 n 個點，建一棵 **2-d tree**，支援五種空間查詢，並與 O(n) 暴力掃描
    逐項對拍，確保每一次剪枝都沒剪掉真正的答案。

    操作碼（每行 op 及其參數）：
        1  qx qy                 最近鄰：輸出「距離平方 點編號」（平手取編號最小）
        2  qx qy k               k 近鄰：輸出前 k 個距離平方（升序；不足補 -1）
        3  x1 y1 x2 y2           矩形範圍計數（邊界含在內）
        4  qx qy r               半徑計數：距離平方 ≤ r² 的點數
        5  qx qy x1 y1 x2 y2    帶矩形約束的最近鄰（只考慮落在矩形內的點）

思路：
    ### 1. 為什麼二分空間也能「二分搜尋」
    BST 是對**一維**鍵值二分；k-d tree 把它推廣成輪流對**各個維度**二分。
    本題在 2D：每次建節點時，看當前這群點在哪個維度上**散得最開**（spread 最大），
    就沿那個維度排序取中位數切開——這樣切出來的「格子」最接近正方形，
    剪枝效果最好（固定輪流切 x/y 也可以，但遇到狹長分佈會退化）。

    每個節點額外存**這棵子樹所有點的包圍盒**（bounding box, 4 個數字）。
    有了包圍盒，就能算出「查詢點到這整棵子樹的**最小可能距離**」，
    這正是剪枝的依據。

    ### 2. 最近鄰的剪枝（本題核心）
    遞迴時先走「離查詢點比較近」的那一邊，拿到一個候選答案 best 之後，
    再決定另一邊要不要走：

        若 minDist(查詢點, 另一邊的包圍盒) > best  →  整棵子樹不可能更近，剪掉

    注意必須是 `>` 才剪（不是 `>=`）：平手時我們要取**編號最小**的點，
    所以距離相等的子樹仍必須進去看一眼。這是兩語言要一致的第一個細節。

    ### 3. k 近鄰
    維護一個大小 k 的有序緩衝區，第 k 名就是當前的「門檻」。
    剪枝條件一樣用 `>`：包圍盒最小距離等於門檻時仍要進去，
    因為可能存在「距離一樣但編號更小」的點，會把目前的第 k 名擠掉。

    ### 4. 兩語言一致性的三個細節
    - 切分維度平手時固定選 x（`spread_x >= spread_y` → axis 0）。
    - 排序鍵必須是 `(該維度座標, 點編號)` 這種**全序**，不能只比座標，
      否則座標相同的點順序在兩邊可能不同，中位數就換人，整棵樹的形狀跟著不同。
    - 矩形查詢先把 x1>x2、y1>y2 的輸入正規化（min/max），
      半徑 r 先夾到非負，避免兩個語言對「退化矩形」有不同解讀。

    ### 5. 複雜度：是「期望」不是「最壞」
    k-d tree 的查詢在點均勻分佈時是 O(log n)；但遇到點全落在圓周上、
    或查詢點離所有點都很近的最壞情況，仍會退化成 O(n)。
    要嚴格的最壞對數複雜度得用 R-tree / BVH / Voronoi，
    不過實務上 k-d tree 因為常數極小、又好寫，仍是第一選擇。

應用場景：
    - 地理資訊：外送 / 打車配對「最近的司機」、附近餐廳、電子圍欄
    - 機器學習：kNN 分類與迴歸、DBSCAN / 密度聚類的鄰域查詢加速
    - 電腦圖學：光線追蹤的加速結構（k-d tree 是 BVH 的老前輩）、碰撞檢測
    - 電腦視覺：SIFT / ORB 特徵點的最近鄰匹配（FLANN 底層就有 k-d tree）
    - 天文學：星表交叉比對（cross-match），把兩份上百萬筆的星表對起來
    - 物理模擬：SPH / 分子動力學的鄰居清單（neighbor list）

複雜度（n 個點、2 維）：
    建樹                 O(n log² n)（每層排序；用 nth_element 可降到 O(n log n)）
    最近鄰查詢           平均 O(log n)，最壞 O(n)
    k 近鄰查詢           平均 O(k log n)
    矩形 / 半徑計數       平均 O(√n + 輸出數)，最壞 O(n)
    空間                 O(n)
"""

import random
import re
import sys
from typing import List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 400          # IO 模式下的點數上限
QMAX = 200          # IO 模式下的查詢數上限
KMAX = 8            # k 近鄰的 k 上限
CMAX = 10 ** 6      # 座標絕對值上限
RMAX = 3 * 10 ** 6  # 半徑上限

INF = 1 << 62


# ---------------------------------------------------------------- k-d tree


class KDNode:
    """一個 k-d tree 節點：代表的點編號 + 左右子樹 + 子樹包圍盒。"""

    __slots__ = ("idx", "left", "right", "minx", "maxx", "miny", "maxy")

    def __init__(self, idx: int) -> None:
        self.idx = idx
        self.left: Optional["KDNode"] = None
        self.right: Optional["KDNode"] = None
        self.minx = 0
        self.maxx = 0
        self.miny = 0
        self.maxy = 0


class KDTree:
    def __init__(self, pts: List[Tuple[int, int]]) -> None:
        self.pts = pts
        self.root = self._build(list(range(len(pts)))) if pts else None

    def _build(self, ids: List[int]) -> Optional[KDNode]:
        if not ids:
            return None
        minx = min(self.pts[i][0] for i in ids)
        maxx = max(self.pts[i][0] for i in ids)
        miny = min(self.pts[i][1] for i in ids)
        maxy = max(self.pts[i][1] for i in ids)
        # 散得最開的維度優先切；平手時選 x（axis 0）
        axis = 0 if (maxx - minx) >= (maxy - miny) else 1
        # 排序鍵必須是全序（座標, 編號），否則兩語言的中位數可能不同
        ids.sort(key=lambda i: (self.pts[i][axis], i))
        mid = len(ids) // 2
        node = KDNode(ids[mid])
        node.left = self._build(ids[:mid])
        node.right = self._build(ids[mid + 1:])
        node.minx, node.maxx, node.miny, node.maxy = minx, maxx, miny, maxy
        return node

    # -- 包圍盒到點的最小距離平方 --

    @staticmethod
    def _box_dist(node: Optional[KDNode], qx: int, qy: int) -> int:
        if node is None:
            return INF
        dx = 0
        if qx < node.minx:
            dx = node.minx - qx
        elif qx > node.maxx:
            dx = qx - node.maxx
        dy = 0
        if qy < node.miny:
            dy = node.miny - qy
        elif qy > node.maxy:
            dy = qy - node.maxy
        return dx * dx + dy * dy

    @staticmethod
    def _box_hit(node: Optional[KDNode], x1: int, y1: int, x2: int, y2: int) -> bool:
        """包圍盒是否與矩形 [x1,x2]×[y1,y2] 有交集。"""
        if node is None:
            return False
        return not (node.maxx < x1 or node.minx > x2 or node.maxy < y1 or node.miny > y2)

    def _dist2(self, i: int, qx: int, qy: int) -> int:
        dx = self.pts[i][0] - qx
        dy = self.pts[i][1] - qy
        return dx * dx + dy * dy

    # -- 最近鄰 --

    def nearest(self, qx: int, qy: int) -> Tuple[int, int]:
        self._bd = INF
        self._bi = -1
        self._nearest(self.root, qx, qy)
        if self._bi < 0:
            return (-1, -1)
        return (self._bd, self._bi)

    def _nearest(self, node: Optional[KDNode], qx: int, qy: int) -> None:
        if node is None:
            return
        d = self._dist2(node.idx, qx, qy)
        if d < self._bd or (d == self._bd and node.idx < self._bi):
            self._bd = d
            self._bi = node.idx
        dl = self._box_dist(node.left, qx, qy)
        dr = self._box_dist(node.right, qx, qy)
        first, second = (node.left, node.right) if dl <= dr else (node.right, node.left)
        fd, sd = (dl, dr) if dl <= dr else (dr, dl)
        # 用 `<=` 才剪（不是 `<`）：平手時還要進去找編號更小的點
        if fd <= self._bd:
            self._nearest(first, qx, qy)
        if sd <= self._bd:
            self._nearest(second, qx, qy)

    # -- k 近鄰 --

    def knn(self, qx: int, qy: int, k: int) -> List[int]:
        self._buf: List[Tuple[int, int]] = []
        self._k = k
        self._knn(self.root, qx, qy)
        self._buf.sort()
        return [d for d, _ in self._buf]

    def _knn(self, node: Optional[KDNode], qx: int, qy: int) -> None:
        if node is None:
            return
        d = self._dist2(node.idx, qx, qy)
        self._buf.append((d, node.idx))
        self._buf.sort()
        if len(self._buf) > self._k:
            self._buf = self._buf[:self._k]
        thresh = self._buf[-1][0] if len(self._buf) >= self._k else INF
        dl = self._box_dist(node.left, qx, qy)
        dr = self._box_dist(node.right, qx, qy)
        first, second = (node.left, node.right) if dl <= dr else (node.right, node.left)
        fd, sd = (dl, dr) if dl <= dr else (dr, dl)
        # 等於門檻時仍要進去：可能有「距離一樣但編號更小」的點擠掉第 k 名
        if fd <= thresh:
            self._knn(first, qx, qy)
            thresh = self._buf[-1][0] if len(self._buf) >= self._k else INF
        if sd <= thresh:
            self._knn(second, qx, qy)

    # -- 矩形範圍計數 --

    def range_count(self, x1: int, y1: int, x2: int, y2: int) -> int:
        if x1 > x2:
            x1, x2 = x2, x1
        if y1 > y2:
            y1, y2 = y2, y1
        self._rc = (x1, y1, x2, y2)
        self._cnt = 0
        self._range(self.root)
        return self._cnt

    def _range(self, node: Optional[KDNode]) -> None:
        if node is None:
            return
        x1, y1, x2, y2 = self._rc
        if not self._box_hit(node, x1, y1, x2, y2):
            return
        px, py = self.pts[node.idx]
        if x1 <= px <= x2 and y1 <= py <= y2:
            self._cnt += 1
        self._range(node.left)
        self._range(node.right)

    # -- 半徑計數 --

    def radius_count(self, qx: int, qy: int, r: int) -> int:
        if r < 0:
            r = 0
        self._rq = (qx, qy)
        self._r2 = r * r
        self._cnt = 0
        self._radius(self.root)
        return self._cnt

    def _radius(self, node: Optional[KDNode]) -> None:
        if node is None:
            return
        qx, qy = self._rq
        if self._box_dist(node, qx, qy) > self._r2:
            return
        if self._dist2(node.idx, qx, qy) <= self._r2:
            self._cnt += 1
        self._radius(node.left)
        self._radius(node.right)

    # -- 帶矩形約束的最近鄰 --

    def nearest_in_rect(self, qx: int, qy: int, x1: int, y1: int, x2: int, y2: int) -> Tuple[int, int]:
        if x1 > x2:
            x1, x2 = x2, x1
        if y1 > y2:
            y1, y2 = y2, y1
        self._rc2 = (x1, y1, x2, y2)
        self._bd = INF
        self._bi = -1
        self._nearest_rect(self.root, qx, qy)
        if self._bi < 0:
            return (-1, -1)
        return (self._bd, self._bi)

    def _nearest_rect(self, node: Optional[KDNode], qx: int, qy: int) -> None:
        if node is None:
            return
        x1, y1, x2, y2 = self._rc2
        if not self._box_hit(node, x1, y1, x2, y2):
            return
        px, py = self.pts[node.idx]
        if x1 <= px <= x2 and y1 <= py <= y2:
            d = self._dist2(node.idx, qx, qy)
            if d < self._bd or (d == self._bd and node.idx < self._bi):
                self._bd = d
                self._bi = node.idx
        dl = self._box_dist(node.left, qx, qy)
        dr = self._box_dist(node.right, qx, qy)
        first, second = (node.left, node.right) if dl <= dr else (node.right, node.left)
        fd, sd = (dl, dr) if dl <= dr else (dr, dl)
        if fd <= self._bd:
            self._nearest_rect(first, qx, qy)
        if sd <= self._bd:
            self._nearest_rect(second, qx, qy)


# ---------------------------------------------------------------- 暴力對拍


def bf_nearest(pts, qx, qy):
    best = (INF, -1)
    for i, (x, y) in enumerate(pts):
        d = (x - qx) ** 2 + (y - qy) ** 2
        cand = (d, i)
        if cand < best:
            best = cand
    if best[1] < 0:
        return (-1, -1)
    return best


def bf_knn(pts, qx, qy, k):
    ds = sorted(((x - qx) ** 2 + (y - qy) ** 2, i) for i, (x, y) in enumerate(pts))
    return [d for d, _ in ds[:k]]


def bf_range(pts, x1, y1, x2, y2):
    if x1 > x2:
        x1, x2 = x2, x1
    if y1 > y2:
        y1, y2 = y2, y1
    return sum(1 for (x, y) in pts if x1 <= x <= x2 and y1 <= y <= y2)


def bf_radius(pts, qx, qy, r):
    if r < 0:
        r = 0
    r2 = r * r
    return sum(1 for (x, y) in pts if (x - qx) ** 2 + (y - qy) ** 2 <= r2)


def bf_nearest_rect(pts, qx, qy, x1, y1, x2, y2):
    if x1 > x2:
        x1, x2 = x2, x1
    if y1 > y2:
        y1, y2 = y2, y1
    best = (INF, -1)
    for i, (x, y) in enumerate(pts):
        if x1 <= x <= x2 and y1 <= y <= y2:
            cand = ((x - qx) ** 2 + (y - qy) ** 2, i)
            if cand < best:
                best = cand
    if best[1] < 0:
        return (-1, -1)
    return best


# ---------------------------------------------------------------- 主流程


def run(pts: List[Tuple[int, int]], qs: List[List[int]]) -> str:
    tree = KDTree(pts)
    out: List[str] = []
    ok = 1

    for q in qs:
        op = q[0]
        if op == 1:
            _, qx, qy = q
            got = tree.nearest(qx, qy)
            want = bf_nearest(pts, qx, qy)
            out.append(f"{got[0]} {got[1]}")
            if got != want:
                ok = 0
        elif op == 2:
            _, qx, qy, k = q
            got = tree.knn(qx, qy, k)
            want = bf_knn(pts, qx, qy, k)
            got2 = got + [-1] * (k - len(got))
            want2 = want + [-1] * (k - len(want))
            out.append(" ".join(str(v) for v in got2))
            if got2 != want2:
                ok = 0
        elif op == 3:
            _, x1, y1, x2, y2 = q
            got = tree.range_count(x1, y1, x2, y2)
            want = bf_range(pts, x1, y1, x2, y2)
            out.append(str(got))
            if got != want:
                ok = 0
        elif op == 4:
            _, qx, qy, r = q
            got = tree.radius_count(qx, qy, r)
            want = bf_radius(pts, qx, qy, r)
            out.append(str(got))
            if got != want:
                ok = 0
        elif op == 5:
            _, qx, qy, x1, y1, x2, y2 = q
            got = tree.nearest_in_rect(qx, qy, x1, y1, x2, y2)
            want = bf_nearest_rect(pts, qx, qy, x1, y1, x2, y2)
            out.append(f"{got[0]} {got[1]}")
            if got != want:
                ok = 0

    return "\n".join(out) + ("\n" if out else "") + str(ok) + "\n"


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例：五個點 ----
    pts = [(0, 0), (3, 4), (10, 10), (-5, 2), (1, 1)]
    out = run(pts, [[1, 0, 0], [1, 3, 4], [2, 0, 0, 3], [3, -10, -10, 2, 2],
                    [4, 0, 0, 5], [5, 0, 0, -10, -10, 2, 2]])
    lines = out.rstrip("\n").split("\n")
    assert lines[0] == "0 0", lines
    assert lines[1] == "0 1", lines
    assert lines[2] == "0 2 25", lines
    assert lines[3] == "3", lines          # (0,0) (1,1) (-5,2)
    assert lines[4] == "3", lines          # (0,0) (1,1) (3,4)
    assert lines[5] == "0 0", lines
    assert lines[6] == "1", lines

    # ---- 空點集 ----
    out = run([], [[1, 0, 0], [2, 0, 0, 2], [3, 0, 0, 5, 5], [4, 0, 0, 1], [5, 0, 0, 0, 0, 9, 9]])
    lines = out.rstrip("\n").split("\n")
    assert lines[0] == "-1 -1", lines
    assert lines[1] == "-1 -1", lines
    assert lines[2] == "0", lines
    assert lines[3] == "0", lines
    assert lines[4] == "-1 -1", lines
    assert lines[5] == "1", lines

    # ---- 全部點重疊（剪枝的退化情況）----
    pts2 = [(2, 2)] * 7
    out = run(pts2, [[1, 2, 2], [1, 0, 0], [2, 2, 2, 4], [3, 2, 2, 2, 2]])
    lines = out.rstrip("\n").split("\n")
    assert lines[0] == "0 0", lines        # 平手取編號最小
    assert lines[1] == "8 0", lines
    assert lines[2] == "0 0 0 0", lines
    assert lines[3] == "7", lines
    assert lines[4] == "1", lines

    # ---- 隨機對拍（一般分佈）----
    rng = random.Random(20261020)
    for t in range(300):
        n = rng.randint(0, 60)
        pts3 = [(rng.randint(-50, 50), rng.randint(-50, 50)) for _ in range(n)]
        qs = []
        for _ in range(rng.randint(1, 12)):
            op = rng.randint(1, 5)
            if op == 1:
                qs.append([1, rng.randint(-60, 60), rng.randint(-60, 60)])
            elif op == 2:
                qs.append([2, rng.randint(-60, 60), rng.randint(-60, 60), rng.randint(1, 6)])
            elif op == 3:
                qs.append([3, rng.randint(-60, 60), rng.randint(-60, 60),
                           rng.randint(-60, 60), rng.randint(-60, 60)])
            elif op == 4:
                qs.append([4, rng.randint(-60, 60), rng.randint(-60, 60), rng.randint(0, 40)])
            else:
                qs.append([5, rng.randint(-60, 60), rng.randint(-60, 60),
                           rng.randint(-60, 60), rng.randint(-60, 60),
                           rng.randint(-60, 60), rng.randint(-60, 60)])
        out = run(pts3, qs)
        assert out.rstrip("\n").split("\n")[-1] == "1", (t, pts3, qs, out)

    # ---- 隨機對拍（退化分佈：全部在同一條水平線 / 圓周上）----
    for t in range(150):
        n = rng.randint(1, 50)
        style = t % 3
        if style == 0:
            pts4 = [(rng.randint(-40, 40), 0) for _ in range(n)]
        elif style == 1:
            pts4 = [(rng.randint(-40, 40), rng.randint(-40, 40)) for _ in range(n)]
            pts4 = [(x, 3 * x + 1) for (x, _) in pts4]
        else:
            import math
            pts4 = []
            for i in range(n):
                a = 2 * math.pi * i / max(1, n)
                pts4.append((int(1000 * math.cos(a)), int(1000 * math.sin(a))))
        qs = []
        for _ in range(8):
            op = rng.randint(1, 5)
            if op == 1:
                qs.append([1, rng.randint(-1100, 1100), rng.randint(-1100, 1100)])
            elif op == 2:
                qs.append([2, rng.randint(-1100, 1100), rng.randint(-1100, 1100), rng.randint(1, 5)])
            elif op == 3:
                qs.append([3, rng.randint(-1100, 1100), rng.randint(-1100, 1100),
                           rng.randint(-1100, 1100), rng.randint(-1100, 1100)])
            elif op == 4:
                qs.append([4, rng.randint(-1100, 1100), rng.randint(-1100, 1100), rng.randint(0, 600)])
            else:
                qs.append([5, rng.randint(-1100, 1100), rng.randint(-1100, 1100),
                           rng.randint(-1100, 1100), rng.randint(-1100, 1100),
                           rng.randint(-1100, 1100), rng.randint(-1100, 1100)])
        out = run(pts4, qs)
        assert out.rstrip("\n").split("\n")[-1] == "1", (t, style, out)

    # ---- 大規模：500 個點、200 次查詢 ----
    rng2 = random.Random(777)
    pts5 = [(rng2.randint(-10 ** 6, 10 ** 6), rng2.randint(-10 ** 6, 10 ** 6)) for _ in range(500)]
    qs = [[1, rng2.randint(-10 ** 6, 10 ** 6), rng2.randint(-10 ** 6, 10 ** 6)] for _ in range(100)]
    qs += [[4, rng2.randint(-10 ** 6, 10 ** 6), rng2.randint(-10 ** 6, 10 ** 6),
            rng2.randint(0, 3 * 10 ** 5)] for _ in range(100)]
    out = run(pts5, qs)
    assert out.rstrip("\n").split("\n")[-1] == "1", out


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出，避免 Windows 把 \\n 換成 \\r\\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    toks = raw.split()
    idx = 0

    def next_int(default: int) -> int:
        nonlocal idx
        if idx >= len(toks):
            return default
        t = toks[idx]
        idx += 1
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-(10 ** 18), min(10 ** 18, v))

    def cl(v: int) -> int:
        return max(-CMAX, min(CMAX, v))

    n = next_int(0)
    if n < 0:
        n = 0
    if n > NMAX:
        n = NMAX
    pts: List[Tuple[int, int]] = []
    for _ in range(n):
        x = cl(next_int(0))
        y = cl(next_int(0))
        pts.append((x, y))
    q = next_int(0)
    if q < 0:
        q = 0
    if q > QMAX:
        q = QMAX
    qs: List[List[int]] = []
    for _ in range(q):
        op = next_int(0)
        if op == 1:
            qs.append([1, cl(next_int(0)), cl(next_int(0))])
        elif op == 2:
            # 順序是 `2 qx qy k`，座標先讀，最後才是 k
            qx, qy = cl(next_int(0)), cl(next_int(0))
            k = next_int(1)
            if k < 1:
                k = 1
            if k > KMAX:
                k = KMAX
            qs.append([2, qx, qy, k])
        elif op == 3:
            qs.append([3, cl(next_int(0)), cl(next_int(0)), cl(next_int(0)), cl(next_int(0))])
        elif op == 4:
            # 順序是 `4 qx qy r`
            qx, qy = cl(next_int(0)), cl(next_int(0))
            r = next_int(0)
            if r < 0:
                r = 0
            if r > RMAX:
                r = RMAX
            qs.append([4, qx, qy, r])
        elif op == 5:
            qs.append([5, cl(next_int(0)), cl(next_int(0)),
                       cl(next_int(0)), cl(next_int(0)), cl(next_int(0)), cl(next_int(0))])
        else:
            continue
    emit(run(pts, qs))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
