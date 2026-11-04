"""最小圓覆蓋（Smallest Enclosing Circle / Welzl Move-to-Front）

題意：
    給定平面上 n 個整數點（|x|, |y| <= 10^4），求覆蓋所有點的
    半徑最小的圓，輸出圓心座標與半徑（四捨五入到小數點後 2 位，
    四捨五入採「遠離零」規則；-0.00 正規化為 0.00）。
    輸入：第一行 n（1 <= n <= 5000），接下來 n 行每行 x y。
    輸出：一行「cx cy r」，以單一空格分隔，各保留 2 位小數。

思路：
    ### 1. 三個幾何事實（最壞情形的結構）
    最小覆蓋圓唯一，且由至多 3 個點決定：
        * 由 1 個點決定：只有一個點（半徑 0）
        * 由 2 個點決定：以兩點為直徑的圓
        * 由 3 個點決定：三點的外接圓（三角形為銳角或直角）

    ### 2. Welzl 增量算法（move-to-front 隨機化版本）
    把點隨機打亂後逐個加入：
        * 若當前圓已包含新點，什麼都不做；
        * 否則新點 p 一定在**最終**圓的邊界上 —— 帶著「p 在邊界上」
          這個約束，再對前面的點做一遍同樣的增量（一層遞迴）；
          若又遇到不被包含的點 q，則 p、q 都在邊界上，
          帶著雙約束再掃一遍前面的點（此時每個不被包含的點 s
          依叉積分側，只考慮能包住 p、q、s 的外接圓中半徑最小者）。
    期望時間 O(n)（move-to-front 實現；樸素遞迴版是期望 O(n) 但
    有遞迴深度問題，本題用迭代三層循環實現，行為等價）。

    ### 3. 數值細節
    * 洗牌用自製 LCG（線性同餘產生器），Python 與 C++ 產生
      **完全相同**的點序，保證兩語言輸出逐位元組一致。
    * 外接圓計算先把座標平移到點集範圍中心再算，提升精度。
    * 共線三點沒有外接圓（叉積行列式為 0），回傳 None 由上層跳過。

    ### 4. 對拍
    內建測試除了固定用例（單點、兩點、直角三角形外接圓、正方形、
    共線、重複點）外，還用 O(n^4) 暴力：枚舉所有點對（直徑圓）與
    點三元組（外接圓），取「包含全部點的半徑最小者」，與 Welzl 結果比對。

複雜度：
    期望 O(n)（隨機化增量），最壞 O(n^3)（幾乎不會發生）
    空間 O(n)

應用場景（Applications）：
    - **地理圍欄 / 配送半徑**：給定一批送貨點，求最小服務半徑與調度中心
    - **蜂窩網路基站選址**：最小圓覆蓋 = 覆蓋所有用戶的最小基站覆蓋半徑
    - **製造業公差檢測**：一批孔位 / 銷釘的最小包絡圓 = 定位誤差上界
    - **電腦圖學**：包圍球（bounding sphere）預計算，加速視錐剔除與碰撞檢測
"""

import math
import sys
from typing import List, Optional, Sequence, Tuple

NMAX = 5000
XYMAX = 10 ** 4
EPS = 1e-14

Point = Tuple[int, int]
Circle = Tuple[float, float, float]  # (cx, cy, r)


def is_in_circle(c: Optional[Circle], x: float, y: float) -> bool:
    """點是否在圓內（含邊界；容差 eps 相對於半徑）。"""
    if c is None:
        return False
    dx = x - c[0]
    dy = y - c[1]
    return dx * dx + dy * dy <= c[2] * c[2] * (1.0 + EPS)


def make_diameter(ax: float, ay: float, bx: float, by: float) -> Circle:
    """以 a、b 為直徑的圓。"""
    cx = (ax + bx) / 2.0
    cy = (ay + by) / 2.0
    dx = ax - bx
    dy = ay - by
    r = math.sqrt(dx * dx + dy * dy) / 2.0
    return (cx, cy, r)


def make_circumcircle(ax: float, ay: float, bx: float, by: float,
                      cx: float, cy: float) -> Optional[Circle]:
    """三點外接圓；三點共線回傳 None。"""
    # 平移到座標範圍中心，提升數值精度
    ox = (min(ax, bx, cx) + max(ax, bx, cx)) / 2.0
    oy = (min(ay, by, cy) + max(ay, by, cy)) / 2.0
    ax2 = ax - ox
    ay2 = ay - oy
    bx2 = bx - ox
    by2 = by - oy
    cx2 = cx - ox
    cy2 = cy - oy
    d = (ax2 * (by2 - cy2) + bx2 * (cy2 - ay2) + cx2 * (ay2 - by2)) * 2.0
    if d == 0.0:
        return None
    aa = ax2 * ax2 + ay2 * ay2
    bb = bx2 * bx2 + by2 * by2
    cc = cx2 * cx2 + cy2 * cy2
    x = ox + (aa * (by2 - cy2) + bb * (cy2 - ay2) + cc * (ay2 - by2)) / d
    y = oy + (aa * (cx2 - bx2) + bb * (ax2 - cx2) + cc * (bx2 - ax2)) / d
    r1 = math.sqrt((x - ax) * (x - ax) + (y - ay) * (y - ay))
    r2 = math.sqrt((x - bx) * (x - bx) + (y - by) * (y - by))
    r3 = math.sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy))
    r = max(r1, r2, r3)
    return (x, y, r)


def make_circle_two_points(points: Sequence[Point],
                           px: float, py: float, qx: float, qy: float) -> Circle:
    """求覆蓋 points 且邊界經過 p、q 的最小圓。"""
    circ = make_diameter(px, py, qx, qy)
    left: Optional[Circle] = None
    right: Optional[Circle] = None
    left_side = 0.0   # left 圓心相對弦 p→q 的叉積（>0，越遠越大）
    right_side = 0.0  # right 圓心相對弦 p→q 的叉積（<0，越遠越小）
    for (sx, sy) in points:
        if is_in_circle(circ, float(sx), float(sy)):
            continue
        cross = (qx - px) * (sy - py) - (qy - py) * (sx - px)
        c = make_circumcircle(px, py, qx, qy, float(sx), float(sy))
        if c is None:
            continue
        # 候選圓心相對弦 p→q 的叉積（帶符號距離 × |pq|）
        c_side = (qx - px) * (c[1] - py) - (qy - py) * (c[0] - px)
        # 同側候選中取「圓心離弦最遠」者，才能包住同側所有點
        if cross > 0.0 and (left is None or c_side > left_side):
            left = c
            left_side = c_side
        elif cross < 0.0 and (right is None or c_side < right_side):
            right = c
            right_side = c_side
    if left is None and right is None:
        return circ
    elif left is None:
        return right
    elif right is None:
        return left
    else:
        return left if left[2] <= right[2] else right


def make_circle_one_point(points: Sequence[Point],
                          px: float, py: float) -> Circle:
    """求覆蓋 points 且邊界經過 p 的最小圓。"""
    c: Circle = (float(px), float(py), 0.0)
    for i, (sx, sy) in enumerate(points):
        if not is_in_circle(c, float(sx), float(sy)):
            if c[2] == 0.0:
                c = make_diameter(float(px), float(py), float(sx), float(sy))
            else:
                c = make_circle_two_points(points[:i + 1],
                                           float(px), float(py), float(sx), float(sy))
    return c


def lcg_shuffle(points: Sequence[Point], seed: int) -> List[Point]:
    """自製 LCG 的 Fisher-Yates 洗牌（兩語言序列一致）。"""
    a = list(points)
    state = seed & 0x7FFFFFFF
    for i in range(len(a) - 1, 0, -1):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        j = state % (i + 1)
        a[i], a[j] = a[j], a[i]
    return a


def smallest_enclosing_circle(points: Sequence[Point]) -> Circle:
    """Welzl move-to-front：回傳 (cx, cy, r)。points 非空。"""
    shuffled = lcg_shuffle(points, 20261103)
    c: Optional[Circle] = None
    for i, (px, py) in enumerate(shuffled):
        if c is None or not is_in_circle(c, float(px), float(py)):
            c = make_circle_one_point(shuffled[:i + 1], float(px), float(py))
    assert c is not None
    return c


# ---- 暴力對拍（O(n^4)，僅小規模）----

def brute_force(points: Sequence[Point]) -> float:
    """最小覆蓋半徑：枚舉所有點對（直徑圓）與點三元組（外接圓）。"""
    n = len(points)
    best = float("inf")

    def covers(cx: float, cy: float, r: float) -> bool:
        tol = max(r, 1.0) * 1e-7
        for (x, y) in points:
            dx = x - cx
            dy = y - cy
            if dx * dx + dy * dy > r * r + tol:
                return False
        return True

    for i in range(n):
        for j in range(i, n):
            ax, ay = points[i]
            bx, by = points[j]
            cx = (ax + bx) / 2.0
            cy = (ay + by) / 2.0
            dx = ax - bx
            dy = ay - by
            r = math.sqrt(dx * dx + dy * dy) / 2.0
            if covers(cx, cy, r):
                best = min(best, r)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                c = make_circumcircle(float(points[i][0]), float(points[i][1]),
                                      float(points[j][0]), float(points[j][1]),
                                      float(points[k][0]), float(points[k][1]))
                if c is None:
                    continue
                if covers(c[0], c[1], c[2]):
                    best = min(best, c[2])
    return best


# ---- 輸出格式化（兩語言完全一致）----

def round2(x: float) -> float:
    """四捨五入到 2 位小數（遠離零），並把 -0.0 正規化成 0.0。"""
    if x >= 0.0:
        v = float(int(x * 100.0 + 0.5)) / 100.0
    else:
        v = float(int(x * 100.0 - 0.5)) / 100.0
    if v == 0.0:
        return 0.0
    return v


def fmt(x: float) -> str:
    return "%.2f" % round2(x)


def run_tests() -> None:
    # ---- 固定用例（先跑過拿到真實輸出再寫死）----
    assert smallest_enclosing_circle([(3, 4)]) == (3.0, 4.0, 0.0)
    assert smallest_enclosing_circle([(0, 0), (4, 0)]) == (2.0, 0.0, 2.0)
    # 直角三角形：外接圓以斜邊為直徑，圓心 (2, 1.5)、半徑 2.5
    assert smallest_enclosing_circle([(0, 0), (4, 0), (0, 3)]) == (2.0, 1.5, 2.5)
    # 正方形：圓心 (1,1)、半徑 sqrt(2)
    c = smallest_enclosing_circle([(0, 0), (2, 0), (2, 2), (0, 2)])
    assert abs(c[0] - 1.0) < 1e-9 and abs(c[1] - 1.0) < 1e-9
    assert abs(c[2] - 2 ** 0.5) < 1e-9
    # 共線三點：兩端點為直徑
    assert smallest_enclosing_circle([(-1, 0), (0, 0), (5, 0)]) == (2.0, 0.0, 3.0)
    # 重複點
    assert smallest_enclosing_circle([(2, 2), (2, 2), (2, 2)]) == (2.0, 2.0, 0.0)

    # ---- 與暴力對拍（自製 LCG，兩語言同一序列）----
    state = 20261103

    def nxt(mod: int) -> int:
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mod

    for _ in range(250):
        n = 1 + nxt(14)  # 1..14
        pts = [(nxt(41) - 20, nxt(41) - 20) for _ in range(n)]  # [-20, 20]
        c = smallest_enclosing_circle(pts)
        # 包容性：所有點都在圓內（含容差）
        for (x, y) in pts:
            dx = x - c[0]
            dy = y - c[1]
            assert dx * dx + dy * dy <= c[2] * c[2] * (1.0 + 1e-9), (pts, c, x, y)
        # 半徑與暴力一致
        want = brute_force(pts)
        assert abs(c[2] - want) <= max(want, 1.0) * 1e-6, (pts, c, want)

    # 座標範圍較大（測精度）
    for _ in range(80):
        n = 1 + nxt(10)
        pts = [(nxt(20001) - 10000, nxt(20001) - 10000) for _ in range(n)]
        c = smallest_enclosing_circle(pts)
        for (x, y) in pts:
            dx = x - c[0]
            dy = y - c[1]
            assert dx * dx + dy * dy <= c[2] * c[2] * (1.0 + 1e-9), (pts, c, x, y)
        want = brute_force(pts)
        assert abs(c[2] - want) <= max(want, 1.0) * 1e-6, (pts, c, want)


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


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    n_opt = _parse_int(nxt() or "")
    if n_opt is None or not (1 <= n_opt <= NMAX):
        bad()
        return
    n = n_opt

    pts: List[Point] = []
    for _ in range(n):
        x = _parse_int(nxt() or "")
        y = _parse_int(nxt() or "")
        if x is None or y is None or abs(x) > XYMAX or abs(y) > XYMAX:
            bad()
            return
        pts.append((x, y))

    c = smallest_enclosing_circle(pts)
    emit("%s %s %s\n" % (fmt(c[0]), fmt(c[1]), fmt(c[2])))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
