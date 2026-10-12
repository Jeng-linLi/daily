"""計算幾何：凸包（Andrew / Jarvis / Graham）+ 面積 + 旋轉卡殼直徑 + 點在凸包內

題意：
    給定平面上的 n 個點（可重複、可共線），求：
      1. **凸包**：包住所有點的最小凸多邊形的頂點（逆時針，從字典序最小的點開始，
         不含共線的中間點）；
      2. **面積的兩倍**（整數，鞋帶公式，避免浮點誤差）；
      3. **直徑的平方**：凸包上最遠兩點的距離平方（旋轉卡殼，O(h)）；
      4. **點在凸包內判定**：q 個查詢點，判斷是否在凸包內或邊界上（O(log h)）。

思路：
    ### 全程整數運算：為什麼開根號是壞習慣
    幾何題最常見的坑就是浮點誤差（`0.1+0.2 != 0.3`、極限共線時 `acos` 的 NaN）。
    本題**只做整數運算**：叉積、面積的兩倍、距離的平方全部是整數，
    只有「比較」而不需要開根號，因此 Python 與 C++ 的輸出可以做到逐字節一致。

    ### 叉積是幾何的萬能鑰匙
        cross(o, a, b) = (a-o) × (b-o) = (a.x-o.x)(b.y-o.y) - (a.y-o.y)(b.x-o.x)
      - `> 0`：b 在向量 o→a 的**左側**（逆時針）
      - `= 0`：三點共線
      - `< 0`：b 在右側（順時針）
    凸包、面積、點在多邊形內，全部只靠這一個判斷式。

    ### Andrew's Monotone Chain，O(n log n)
    先按 (x, y) 排序。分別構造：
      - **下鏈** left→right：維持「逆時針轉彎」，遇到 `cross ≤ 0` 就彈棧；
      - **上鏈** right→left：同理。
    兩條鏈去掉各自的尾端接起來就是凸包。
    用 `≤ 0` 而不是 `< 0`，效果是**共線的中間點會被彈掉**，只留兩個端點 —— 這讓輸出唯一。

    為什麼是 O(n log n)？排序 O(n log n)，之後每個點最多進棧一次、出棧一次 → O(n)。

    ### 對拍用的另外兩種構造
      - **Jarvis march（禮物包裝）**，O(n·h)：每步選「最順時針」的下一個點，
        共線時取較遠者。h 小時比 Andrew 快。
      - **Graham scan**，O(n log n)：以字典序最小點為極點，按極角排序後單調棧。
        極角排序**不用 atan2**（浮點！），改用叉積比較 + 共線時比距離。
    三種構造互相對拍，任一處寫錯都會立刻暴露。

    ### 旋轉卡殼求直徑，O(h)
    對每條邊 p_i→p_{i+1}，最遠的頂點 j 會**單調前進**（凸包的對蹠點性質），
    所以 j 全程只往前走，總共 O(h)。判斷「要不要再走一步」用三角形面積
    （即叉積絕對值）比較，一樣是整數。

    ### 點在凸包內（含邊界），O(log h)
    凸包是凸的 → 從 p_0 做扇形三角剖分，q 只可能落在其中一個三角形裡。
    先用兩次叉積排除「楔形外」，再對 p_1..p_{h-1} **二分**找 q 落在哪個扇形，
    最後一次叉積判定是否在那條邊的內側。

    ### 確定性（保證 Python 與 C++ 輸出逐字節一致）
      - 先去重再按 (x, y) 排序 → 輸入順序不影響結果；
      - 共線中間點一律剔除（`≤ 0` 彈棧），凸包唯一；
      - 起點固定為字典序最小點，方向固定逆時針；
      - 面積輸出「兩倍」、距離輸出「平方」，全程整數，無浮點。

應用場景：
    地理資訊系統（GIS）的圍欄與邊界繪製、機器人路徑規劃的可行區域、
    碰撞檢測（凸包 + SAT 分離軸定理）、影像處理的輪廓提取與物體識別、
    統計學的資料範圍（ROC 曲線上凸包）、字型渲染的字形外框、
    感測器網路覆蓋範圍、自動駕駛的可行駛區域。

複雜度：
    去重 + 排序        O(n log n)
    Andrew 凸包        O(n log n) 時間、O(n) 空間
    Jarvis / Graham    O(n·h) / O(n log n)，僅供對拍
    面積（鞋帶）       O(h)
    旋轉卡殼直徑       O(h)
    點在凸包內         O(log h) 每次
    暴力對拍           O(h²) / O(n·h)（僅測試用）

輸入格式（stdin，全部以空白分隔）：
    n       點數
    x_i y_i 接下來 2n 個整數為各點座標（座標會被截斷到 ±10^6）
    q       查詢點數
    x y     接下來 2q 個整數為查詢點
    token 不足或非數字時提前結束，缺的部分視為不存在。
輸出格式（stdout）：
    第 1 行：凸包頂點數 h
    第 2 行：凸包頂點（"x,y" 以空格分隔；h=0 時為空行）
    第 3 行：凸包面積的兩倍（整數）
    第 4 行：最遠兩點距離的平方（整數）
    第 5 行：q 個查詢點的判定結果（1=在凸包內或邊上，0=在外部；空格分隔）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from functools import cmp_to_key
from typing import List, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")
COORD_MAX = 10 ** 6          # 座標截斷上限，保證叉積不溢位（|cross| ≤ 8·10^12）
N_MAX = 2000                 # 點數上限，避免病態輸入拖慢

Point = Tuple[int, int]


# ---------------------------------------------------------------- 基本工具


def cross(o: Point, a: Point, b: Point) -> int:
    """(a-o) × (b-o)：>0 表示 b 在 o→a 的左側（逆時針），=0 共線。"""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def dist2(a: Point, b: Point) -> int:
    """兩點距離的平方（整數，不開根號）。"""
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return dx * dx + dy * dy


def dedup_sorted(points: Sequence[Point]) -> List[Point]:
    """去重後按 (x, y) 排序 —— 這是三個演算法共用的起手式。"""
    return sorted(set(points))


# ---------------------------------------------------------------- 三種凸包構造


def andrew_hull(points: Sequence[Point]) -> List[Point]:
    """Andrew's Monotone Chain，O(n log n)，逆時針、共線點剔除。"""
    pts = dedup_sorted(points)
    if len(pts) <= 2:
        return pts
    lower: List[Point] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: List[Point] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def jarvis_hull(points: Sequence[Point]) -> List[Point]:
    """Jarvis march（禮物包裝），O(n·h)，每步選「最順時針」的點，共線取較遠者。"""
    pts = dedup_sorted(points)
    n = len(pts)
    if n <= 2:
        return pts
    hull: List[Point] = []
    cur = 0                                  # 字典序最小點必定是凸包頂點
    while True:
        hull.append(pts[cur])
        nxt = -1
        for j in range(n):
            if j == cur:
                continue
            if nxt == -1:
                nxt = j
                continue
            c = cross(pts[cur], pts[nxt], pts[j])
            if c < 0:                        # j 更順時針 → 換成 j
                nxt = j
            elif c == 0 and dist2(pts[cur], pts[j]) > dist2(pts[cur], pts[nxt]):
                nxt = j                      # 共線 → 取較遠的，跳過中間點
        if nxt == 0:
            break
        cur = nxt
        if len(hull) > n:                    # 保險絲，正常不會觸發
            break
    return hull


def graham_hull(points: Sequence[Point]) -> List[Point]:
    """Graham scan，O(n log n)，極角排序只用叉積 + 距離，不使用 atan2。"""
    pts = dedup_sorted(points)
    if len(pts) <= 2:
        return pts
    p0 = pts[0]
    rest = pts[1:]

    def cmp_angle(a: Point, b: Point) -> int:
        c = cross(p0, a, b)
        if c > 0:
            return -1                        # a 的極角較小，排在前面
        if c < 0:
            return 1
        da, db = dist2(p0, a), dist2(p0, b)
        return -1 if da < db else (1 if da > db else 0)

    rest = sorted(rest, key=cmp_to_key(cmp_angle))
    stk = [p0, rest[0]]
    for p in rest[1:]:
        while len(stk) >= 2 and cross(stk[-2], stk[-1], p) <= 0:
            stk.pop()
        stk.append(p)
    return stk


# ---------------------------------------------------------------- 凸包上的計算


def hull_area2(hull: Sequence[Point]) -> int:
    """鞋帶公式求面積的**兩倍**（整數，避免浮點）。"""
    h = len(hull)
    if h < 3:
        return 0
    s = 0
    for i in range(h):
        x1, y1 = hull[i]
        x2, y2 = hull[(i + 1) % h]
        s += x1 * y2 - x2 * y1
    return abs(s)


def diameter2(hull: Sequence[Point]) -> int:
    """旋轉卡殼求最遠兩點距離的平方，O(h)。"""
    h = len(hull)
    if h < 2:
        return 0
    if h == 2:
        return dist2(hull[0], hull[1])
    best = 0
    j = 1
    for i in range(h):
        ni = (i + 1) % h
        steps = 0
        while steps < 2 * h:                 # 保險絲：j 全程最多前進 2h 次
            nj = (j + 1) % h
            if abs(cross(hull[i], hull[ni], hull[nj])) > abs(cross(hull[i], hull[ni], hull[j])):
                j = nj
                steps += 1
            else:
                break
        best = max(best, dist2(hull[i], hull[j]), dist2(hull[ni], hull[j]))
    return best


def diameter2_bruteforce(hull: Sequence[Point]) -> int:
    """O(h²) 暴力求最遠點對，僅供對拍。"""
    h = len(hull)
    best = 0
    for i in range(h):
        for j in range(i + 1, h):
            best = max(best, dist2(hull[i], hull[j]))
    return best


def point_in_convex(hull: Sequence[Point], q: Point) -> bool:
    """q 是否在凸包內（含邊界），O(log h)。"""
    h = len(hull)
    if h == 0:
        return False
    if h == 1:
        return q == hull[0]
    if h == 2:
        if cross(hull[0], hull[1], q) != 0:
            return False
        return (min(hull[0][0], hull[1][0]) <= q[0] <= max(hull[0][0], hull[1][0]) and
                min(hull[0][1], hull[1][1]) <= q[1] <= max(hull[0][1], hull[1][1]))
    # 楔形外直接排除：q 必須在 p0→p1 的左側、且在 p0→p_{h-1} 的右側
    if cross(hull[0], hull[1], q) < 0:
        return False
    if cross(hull[0], hull[h - 1], q) > 0:
        return False
    # 二分：找最大的 i（1 ≤ i ≤ h-2）使「q 在 p0→p_i 的左側（含共線）」
    lo, hi = 1, h - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if cross(hull[0], hull[mid], q) >= 0:
            lo = mid
        else:
            hi = mid
    return cross(hull[lo], hull[lo + 1], q) >= 0


def point_in_convex_bruteforce(hull: Sequence[Point], q: Point) -> bool:
    """O(h) 暴力判定：凸多邊形內部 ⟺ 在所有邊的左側（含共線）。"""
    h = len(hull)
    if h == 0:
        return False
    if h == 1:
        return q == hull[0]
    if h == 2:
        return point_in_convex(hull, q)      # 退化情形直接沿用同一套邏輯
    return all(cross(hull[i], hull[(i + 1) % h], q) >= 0 for i in range(h))


# ---------------------------------------------------------------- 輸出


def run(points: Sequence[Point], queries: Sequence[Point]) -> str:
    hull = andrew_hull(points)
    lines = [
        str(len(hull)),
        " ".join(f"{p[0]},{p[1]}" for p in hull),
        str(hull_area2(hull)),
        str(diameter2(hull)),
        " ".join("1" if point_in_convex(hull, q) else "0" for q in queries),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check_hull(points: List[Point]) -> List[Point]:
    """三種構造一致 + 凸包的基本性質 + 直徑與判定與暴力一致。"""
    a = andrew_hull(points)
    j = jarvis_hull(points)
    g = graham_hull(points)
    assert a == j, (points, a, j)
    assert a == g, (points, a, g)

    h = len(a)
    # 凸包頂點必須是原集合的子集
    uniq = set(points)
    for p in a:
        assert p in uniq
    # 嚴格凸（共線中間點已剔除）：相鄰三元的叉積必為正
    if h >= 3:
        for i in range(h):
            assert cross(a[i], a[(i + 1) % h], a[(i + 2) % h]) > 0, (points, a)
    # 所有原點都必須落在凸包內（含邊界）
    for p in uniq:
        assert point_in_convex_bruteforce(a, p), (points, a, p)
    # 直徑：旋轉卡殼 == O(h²) 暴力
    assert diameter2(a) == diameter2_bruteforce(a), (points, a)
    return a


def run_tests() -> None:
    rng = random.Random(20261012)

    # ---- 固定用例（期望值皆跑過之後才寫死）----
    sq = [(0, 0), (10, 0), (10, 10), (0, 10)]
    assert andrew_hull(sq) == [(0, 0), (10, 0), (10, 10), (0, 10)]
    assert hull_area2(andrew_hull(sq)) == 200
    assert diameter2(andrew_hull(sq)) == 200
    check_hull(sq)

    # 內部點與共線點都要被剔除
    check_hull([(0, 0), (10, 0), (10, 10), (0, 10), (5, 5), (5, 0), (0, 5)])
    assert andrew_hull([(0, 0), (10, 0), (10, 10), (0, 10), (5, 5), (5, 0), (0, 5)]) == \
        [(0, 0), (10, 0), (10, 10), (0, 10)]

    # 退化情形
    assert andrew_hull([]) == []
    assert andrew_hull([(3, 4)]) == [(3, 4)]
    assert andrew_hull([(0, 0), (5, 0)]) == [(0, 0), (5, 0)]
    assert andrew_hull([(0, 0), (1, 0), (2, 0), (3, 0)]) == [(0, 0), (3, 0)]
    assert andrew_hull([(1, 1), (1, 1), (1, 1)]) == [(1, 1)]
    assert hull_area2([(0, 0), (3, 0)]) == 0
    assert diameter2([(0, 0), (3, 0)]) == 9
    check_hull([])
    check_hull([(3, 4)])
    check_hull([(0, 0), (5, 0)])
    check_hull([(0, 0), (1, 0), (2, 0), (3, 0)])

    # 負座標
    check_hull([(-5, -5), (5, -5), (5, 5), (-5, 5), (0, 0)])
    assert hull_area2(andrew_hull([(-5, -5), (5, -5), (5, 5), (-5, 5)])) == 200

    # 具體輸出（跑過之後寫死）
    assert run([(0, 0), (10, 0), (10, 10), (0, 10)], [(5, 5), (20, 20), (10, 10)]) == \
        "4\n0,0 10,0 10,10 0,10\n200\n200\n1 0 1"
    assert run([], []) == "0\n\n0\n0\n"

    # ---- 隨機對拍：小座標，共線與重複點大量出現 ----
    for _ in range(600):
        n = rng.randint(0, 10)
        pts = [(rng.randint(-4, 4), rng.randint(-4, 4)) for _ in range(n)]
        check_hull(pts)

    # ---- 隨機對拍：中等座標，凸包較大 ----
    for _ in range(400):
        n = rng.randint(0, 25)
        pts = [(rng.randint(-50, 50), rng.randint(-50, 50)) for _ in range(n)]
        check_hull(pts)

    # ---- 隨機對拍：點在圓上（凸包接近整圓，容易卡到共線邊界）----
    for _ in range(200):
        n = rng.randint(3, 12)
        pts = list({(rng.randint(-6, 6), rng.randint(-6, 6)) for _ in range(n)})
        check_hull(pts)

    # ---- 點在凸包內：與暴力逐點對拍 ----
    for _ in range(400):
        n = rng.randint(0, 12)
        pts = [(rng.randint(-6, 6), rng.randint(-6, 6)) for _ in range(n)]
        hull = andrew_hull(pts)
        for _ in range(20):
            q = (rng.randint(-9, 9), rng.randint(-9, 9))
            assert point_in_convex(hull, q) == point_in_convex_bruteforce(hull, q), (hull, q)


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。"""
    tokens = raw.split()

    def read_int(i: int):
        if i >= len(tokens) or not NUM_RE.match(tokens[i]):
            return None
        return int(tokens[i])

    def clamp(v: int) -> int:
        return max(-COORD_MAX, min(COORD_MAX, v))

    n = read_int(0)
    n = 0 if n is None else max(0, min(n, N_MAX))
    pts: List[Point] = []
    for i in range(n):
        x = read_int(1 + 2 * i)
        y = read_int(2 + 2 * i)
        if x is None or y is None:
            break
        pts.append((clamp(x), clamp(y)))

    base = 1 + 2 * n
    q = read_int(base)
    q = 0 if q is None else max(0, min(q, N_MAX))
    queries: List[Point] = []
    for i in range(q):
        x = read_int(base + 1 + 2 * i)
        y = read_int(base + 2 + 2 * i)
        if x is None or y is None:
            break
        queries.append((clamp(x), clamp(y)))

    sys.stdout.write(run(pts, queries) + "\n")


def main() -> None:
    raw = sys.stdin.read()
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        print("all tests passed")


if __name__ == "__main__":
    main()
