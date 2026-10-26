"""半平面交（Half-Plane Intersection, HPI）

題意：
    給 n 個半平面，第 i 個寫成 a_i·x + b_i·y ≤ c_i。求它們交集的邊界——
    一個（可能空的、可能退化的）凸區域。為了讓輸出永遠是**有界**多邊形，
    題目約定把結果裁剪到包圍盒 [-10⁹, 10⁹]² 內輸出。

    IO 格式（有 stdin 時）：
        第 1 行    n
        接下來 n 行  a b c          半平面 a·x + b·y ≤ c
    輸出：
        第 1 行    k               交集中「頂點」的個數（k = 0 表示空集）
        接下來 k 行  xn/xd yn/yd   頂點座標，以最簡分數輸出（分母 > 0），
                                  按逆時針順序、從字典序最小的頂點開始

    例：三個半平面 x ≥ 0、y ≥ 0、x + y ≤ 10 的交集是三角形，
    輸出 3 與 (0,0)、(10,0)、(0,10)。

思路：
    ### 1. 交集長什麼樣
    半平面交 = 凸多邊形（或空集、或線段 / 點這類退化情形）。
    有界凸多邊形的每個頂點，都是「兩條約束邊界線的交點」：
        頂點 p 處至少有兩條約束線同時取等號，且兩條線不平行。
    反過來也成立：任取兩條不平行的約束線，求交，若交點滿足全部約束，
    它一定是極值點（頂點）——因為在交點附近，可行域被夾在兩條交叉線
    形成的楔形裡，交點這個「尖角」抹不平。

    ### 2. 本題的演算法：候選枚舉（O(n³)，精確有理數）
    把全部約束線（n 個輸入半平面 + 4 條包圍盒線）兩兩求交，
    用**精確有理數**檢查交點是否滿足所有約束，倖存者就是全部頂點。
    這個方法直覺、好驗證，且所有座標都是精確分數，兩個語言可以逐位元組一致。

    包围盒的作用：無界區域（例如只剩一個半平面）被盒子裁成有界多邊形，
    盒子的四條線也參與「兩兩求交」，於是盒角與「輸入線 ∩ 盒線」都會被枚舉到。
    輸入範圍保證了：任何有界區域的頂點座標 |value| ≤ 10⁹ < 盒子邊長，
    盒子不會誤裁有界區域；而無界區域必然被盒子裁到（凸集無界 ⟹ 某方向延伸到無窮）。

    ### 3. 兩個對拍參照（Sutherland–Hodgman 裁剪 + 經典雙端棧）
    候選枚舉好懂但慢，正式賽事另有兩種主流做法，本題都實作作為內置對拍：
    - **多邊形裁剪**：以包圍盒為初始凸多邊形，依次被每個半平面裁剪
      （Sutherland–Hodgman，精確有理數）。直接算出「交集裁到盒子」的結果，
      天然處理空集 / 點 / 線段 / 多邊形全部退化情形，O(m·k)。
    - **排序 + 增量雙端棧**：按邊界線極角排序後增量維護，O(m log m)，
      常數小但平行 / 重複 / 空集等邊界情況極易寫錯，因此加了有效性護欄，
      只在結果可信時參與比對。

    ### 4. 退化與一致性細節
    - (0, 0, c)：c ≥ 0 是全平面（跳過）；c < 0 是空集（直接輸出 0 個頂點）。
    - 同向（法向量同向）的半平面只保留最緊的一個；法向量反向的半平面都要保留。
    - 頂點排序：取字典序最小的頂點 v0，其餘頂點按極角繞 v0 排序即得逆時針顺序
      （凸多邊形從任一頂點看其餘頂點，張角 < π，純叉積比較就是全序）。
    - 座標一律以最簡分數輸出，兩個語言都用精確有理數，不經過浮點。

複雜度（m 條約束線，含包圍盒）：
    候選枚舉        O(m³) 時間、O(m) 空間（m ≤ 64，實測毫秒級）
    多邊形裁剪      O(m·k) 時間、O(k) 空間（k 為交集頂點數）
    雙端棧版        O(m log m) 時間
    輸出規模        O(m)

應用場景：
    - 電腦圖學：視錐裁剪（frustum culling）、任意凸多邊形裁剪、CSG 布林運算
    - 機器人與遊戲：雷達 / 視野範圍（visibility region）、導航可行區域
    - 機器學習：線性分類器的可行域、支援向量機的間隔區域
    - 作業研究：線性規劃可行域的外顯表示（半平面交就是 2D 的 H-多面體）
    - 地理資訊：緩衝區交集、行政區合併的幾何核
"""

import functools
import random
import re
import sys
from collections import deque
from fractions import Fraction
from typing import List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 40            # IO 模式的半平面數上限
AMAX = 100           # |a|、|b| 上限
CMAX = 10 ** 4       # |c| 上限
CLAMP_ABS = 10 ** 18
BOX = 10 ** 9        # 包圍盒邊長（有界性保證：|頂點| ≤ 2·10⁶ < BOX）

Plane = Tuple[int, int, int]
Point = Tuple[Fraction, Fraction]


# ---------------------------------------------------------------- 基礎幾何


def line_intersect(p: Plane, q: Plane) -> Optional[Point]:
    """兩條直線 a·x + b·y = c 的交點；平行（det = 0）回傳 None。"""
    a1, b1, c1 = p
    a2, b2, c2 = q
    det = a1 * b2 - a2 * b1
    if det == 0:
        return None
    x = Fraction(c1 * b2 - c2 * b1, det)
    y = Fraction(a1 * c2 - a2 * c1, det)
    return (x, y)


def inside(pl: Plane, pt: Point) -> bool:
    """點 pt 是否在半平面 a·x + b·y ≤ c 內（含邊界）。"""
    a, b, c = pl
    return a * pt[0] + b * pt[1] <= c


def gcd3(a: int, b: int, c: int) -> int:
    from math import gcd
    return gcd(gcd(abs(a), abs(b)), abs(c))


def same_dir(p: Plane, q: Plane) -> bool:
    """兩條約束線的法向量是否同向（cross = 0 且 dot > 0）。"""
    a1, b1, _ = p
    a2, b2, _ = q
    return a1 * b2 - a2 * b1 == 0 and a1 * a2 + b1 * b2 > 0


def tighter(p: Plane, q: Plane) -> bool:
    """同向平行線中，p 的半平面是否比 q 更緊（p 的區域 ⊆ q 的區域）。"""
    a1, b1, c1 = p
    a2, b2, c2 = q
    # (a2, b2) = λ·(a1, b1)，λ > 0；p 更緊 ⟺ c1/λ ≤ c2
    lam = Fraction(a2, a1) if a1 != 0 else Fraction(b2, b1)
    return c1 * lam <= c2


def normalize_planes(planes: List[Plane]) -> Tuple[List[Plane], bool]:
    """清整 + 去重 + 加包圍盒。回傳 (約束線列表, 是否空集)。"""
    empty = False
    cleaned: List[Plane] = []
    for (a, b, c) in planes:
        if a == 0 and b == 0:
            if c < 0:
                empty = True
            continue
        g = gcd3(a, b, c)
        if g > 1:
            a //= g
            b //= g
            c //= g
        cleaned.append((a, b, c))
    if empty:
        return [], True
    # 包圍盒：把無界區域裁成有界多邊形
    cleaned += [(1, 0, BOX), (-1, 0, BOX), (0, 1, BOX), (0, -1, BOX)]
    # 同向去重，保留最緊的
    deduped: List[Plane] = []
    for pl in cleaned:
        merged = False
        for i, ex in enumerate(deduped):
            if same_dir(pl, ex):
                if tighter(pl, ex):
                    deduped[i] = pl
                merged = True
                break
        if not merged:
            deduped.append(pl)
    return deduped, False


# ---------------------------------------------------------------- 候選枚舉


def hpi_vertices(planes: List[Plane]) -> List[Point]:
    """兩兩求交 + 全約束檢查，收集全部頂點（已按座標去重排序）。"""
    verts: List[Point] = []
    m = len(planes)
    for i in range(m):
        ai, bi, _ = planes[i]
        for j in range(i + 1, m):
            pt = line_intersect(planes[i], planes[j])
            if pt is None:
                continue
            ok = True
            for k in range(m):
                if not inside(planes[k], pt):
                    ok = False
                    break
            if ok:
                verts.append(pt)
    verts.sort()
    out: List[Point] = []
    for v in verts:
        if not out or v != out[-1]:
            out.append(v)
    return out


# ---------------------------------------------------------------- 雙端棧版


def angle_cmp(p: Plane, q: Plane) -> int:
    """按邊界線方向 d = (-b, a) 的極角排序（上半平面優先，同半平面內比叉積）。"""
    a1, b1, _ = p
    a2, b2, _ = q
    h1 = 0 if (a1 > 0 or (a1 == 0 and b1 < 0)) else 1
    h2 = 0 if (a2 > 0 or (a2 == 0 and b2 < 0)) else 1
    if h1 != h2:
        return -1 if h1 < h2 else 1
    cr = a1 * b2 - a2 * b1  # cross(d1, d2) > 0 ⟹ d1 極角更小
    if cr > 0:
        return -1
    if cr < 0:
        return 1
    return 0


def hpi_deque(planes: List[Plane]) -> Optional[List[Point]]:
    """經典 O(n log n) 半平面交（排序 + 增量雙端棧）；退化時回傳 None。"""
    ps = sorted(planes, key=functools.cmp_to_key(angle_cmp))
    dq: deque = deque()
    for q in ps:
        while len(dq) >= 2:
            pt = line_intersect(dq[-2], dq[-1])
            if pt is not None and inside(q, pt):
                break
            dq.pop()
        while len(dq) >= 2:
            pt = line_intersect(dq[0], dq[1])
            if pt is not None and inside(q, pt):
                break
            dq.popleft()
        dq.append(q)
    while len(dq) >= 3:
        pt = line_intersect(dq[-2], dq[-1])
        if pt is not None and inside(dq[0], pt):
            break
        dq.pop()
    while len(dq) >= 3:
        pt = line_intersect(dq[0], dq[1])
        if pt is not None and inside(dq[-1], pt):
            break
        dq.popleft()
    if len(dq) < 3:
        return None
    k = len(dq)
    vs: List[Point] = []
    for i in range(k):
        pt = line_intersect(dq[i], dq[(i + 1) % k])
        if pt is None:
            return None
        vs.append(pt)
    vs.sort()
    out: List[Point] = []
    for v in vs:
        if not out or v != out[-1]:
            out.append(v)
    # 有效性護欄：交集為空時經典雙端棧可能吐出幽靈頂點，逐個驗證
    for v in out:
        for pl in planes:
            if not inside(pl, v):
                return None
    return out


# ---------------------------------------------------------------- 多邊形裁剪


def clip_half(poly: List[Point], pl: Plane) -> List[Point]:
    """Sutherland–Hodgman：用半平面 pl 裁剪凸多邊形 poly（精確有理數）。

    凸區域裁剪半平面是精確的；退化輸入（0/1/2 個頂點）也按同一規則處理，
    連續重複點在函數末尾去除，保證線段 / 點區域不會腐化。
    """
    a, b, c = pl
    out: List[Point] = []
    k = len(poly)
    for i in range(k):
        p, q = poly[i], poly[(i + 1) % k]
        d1 = a * p[0] + b * p[1] - c
        d2 = a * q[0] + b * q[1] - c
        pin = d1 <= 0
        qin = d2 <= 0
        if pin:
            out.append(p)
        if pin != qin:
            t = Fraction(d1, d1 - d2)  # p + t·(q - p)，t ∈ (0, 1)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    dedup: List[Point] = []
    for v in out:
        if not dedup or v != dedup[-1]:
            dedup.append(v)
    return dedup


def hpi_clip(planes: List[Plane]) -> List[Point]:
    """以包圍盒為初始多邊形，依次被每個半平面精確裁剪。

    直接算出「交集裁剪到包圍盒」的結果，天然處理空集 / 點 / 線段 / 多邊形，
    回傳去重排序後的頂點列表，作為候選枚舉的對拍參照。
    """
    poly: List[Point] = [
        (Fraction(-BOX), Fraction(-BOX)), (Fraction(BOX), Fraction(-BOX)),
        (Fraction(BOX), Fraction(BOX)), (Fraction(-BOX), Fraction(BOX)),
    ]
    for (a, b, c) in planes:
        if a == 0 and b == 0:
            if c < 0:
                return []
            continue
        poly = clip_half(poly, (a, b, c))
        if not poly:
            return []
    vs = sorted(poly)
    out: List[Point] = []
    for v in vs:
        if not out or v != out[-1]:
            out.append(v)
    return out


# ---------------------------------------------------------------- 輸出排序


def order_vertices(verts: List[Point]) -> List[Point]:
    """凸多邊形頂點：逆時針、從字典序最小點開始。"""
    if len(verts) <= 2:
        return sorted(verts)
    v0 = min(verts)
    rest = [v for v in verts if v != v0]
    vx, vy = v0

    def cmp(p: Point, q: Point) -> int:
        # cross(p - v0, q - v0) > 0 ⟹ q 的極角更大 ⟹ p 在前
        cr = (p[0] - vx) * (q[1] - vy) - (p[1] - vy) * (q[0] - vx)
        if cr > 0:
            return -1
        if cr < 0:
            return 1
        return 0

    rest.sort(key=functools.cmp_to_key(cmp))
    return [v0] + rest


def fmt_frac(f: Fraction) -> str:
    return f"{f.numerator}/{f.denominator}"


def run(planes: List[Plane]) -> str:
    L, empty = normalize_planes(planes)
    if empty:
        return "0\n"
    verts = hpi_vertices(L)
    verts = order_vertices(verts)
    out = [str(len(verts))]
    for (x, y) in verts:
        out.append(f"{fmt_frac(x)} {fmt_frac(y)}")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------- 結構檢查


def on_segment(p: Point, q: Point, r: Point) -> bool:
    """r 是否在線段 pq 上（精確判斷）。"""
    cr = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    if cr != 0:
        return False
    return (min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and
            min(p[1], q[1]) <= r[1] <= max(p[1], q[1]))


def point_in_poly(verts: List[Point], p: Point) -> bool:
    """點是否在（ccw 凸）多邊形內，含邊界。"""
    k = len(verts)
    if k == 0:
        return False
    if k == 1:
        return verts[0] == p
    if k == 2:
        return on_segment(verts[0], verts[1], p)
    for i in range(k):
        x1, y1 = verts[i]
        x2, y2 = verts[(i + 1) % k]
        cr = (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1)
        if cr < 0:
            return False
    return True


def check_structure(L: List[Plane], verts: List[Point]) -> None:
    """輸出多邊形的結構檢查：嚴格凸 + ccw；每條邊貼著某條約束線。"""
    k = len(verts)
    if k >= 3:
        for i in range(k):
            ax = verts[(i + 1) % k][0] - verts[i][0]
            ay = verts[(i + 1) % k][1] - verts[i][1]
            bx = verts[(i + 2) % k][0] - verts[(i + 1) % k][0]
            by = verts[(i + 2) % k][1] - verts[(i + 1) % k][1]
            cr = ax * by - ay * bx
            assert cr > 0, ("not strictly convex ccw", verts)
        for i in range(k):
            p, q = verts[i], verts[(i + 1) % k]
            found = False
            for (a, b, c) in L:
                if (a * p[0] + b * p[1] == c and
                        a * q[0] + b * q[1] == c):
                    found = True
                    break
            assert found, ("edge not on a constraint line", p, q)


def check_membership(orig: List[Plane], verts: List[Point],
                     rng: random.Random) -> None:
    """在區域附近抽點：落在輸出多邊形內 ⟺ 滿足全部原始半平面。"""
    if not verts:
        return
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    lo_x = min(xs) - 3
    hi_x = max(xs) + 3
    lo_y = min(ys) - 3
    hi_y = max(ys) + 3
    ix = (int(lo_x) if lo_x == int(lo_x) else int(lo_x) - 1)
    ax = (int(hi_x) if hi_x == int(hi_x) else int(hi_x) + 1)
    iy = (int(lo_y) if lo_y == int(lo_y) else int(lo_y) - 1)
    ay = (int(hi_y) if hi_y == int(hi_y) else int(hi_y) + 1)
    nondeg = [pl for pl in orig if not (pl[0] == 0 and pl[1] == 0)]
    for _ in range(120):
        x = Fraction(rng.randint(ix, ax))
        y = Fraction(rng.randint(iy, ay))
        inpoly = point_in_poly(verts, (x, y))
        feasible = all(inside(pl, (x, y)) for pl in nondeg)
        assert inpoly == feasible, ("membership", (x, y), verts)


def parse_verts(text: str) -> List[Point]:
    lines = text.rstrip("\n").split("\n")
    k = int(lines[0])
    assert len(lines) == k + 1, ("line count", text)
    out: List[Point] = []
    for i in range(k):
        a, b = lines[1 + i].split()
        xn, xd = a.split("/")
        yn, yd = b.split("/")
        out.append((Fraction(int(xn), int(xd)), Fraction(int(yn), int(yd))))
    return out


# ---------------------------------------------------------------- 內置測試


def gen_planes(rng: random.Random, style: int) -> List[Plane]:
    n = rng.randint(0, 24)
    planes: List[Plane] = []
    if style == 0:
        # 保證非空：先挑內點，再環繞它產生半平面
        x0 = rng.randint(-50, 50)
        y0 = rng.randint(-50, 50)
        dirs = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (1, -1),
                (-1, 1), (-1, -1), (2, 1), (1, 2)]
        for _ in range(n):
            a, b = rng.choice(dirs)
            slack = rng.randint(0, 60)
            planes.append((a, b, a * x0 + b * y0 + slack))
    elif style == 1:
        # 純隨機（多數為空集）
        for _ in range(n):
            while True:
                a = rng.randint(-5, 5)
                b = rng.randint(-5, 5)
                if a or b:
                    break
            planes.append((a, b, rng.randint(-30, 30)))
    elif style == 2:
        # 小係數、強制平行 / 同向重複
        dirs = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (2, 2),
                (-1, -1), (1, 2), (2, 4)]
        for _ in range(n):
            a, b = rng.choice(dirs)
            planes.append((a, b, rng.randint(-8, 8)))
    else:
        # 夾帶退化平面 (0, 0, c)
        for _ in range(n):
            if rng.random() < 0.15:
                planes.append((0, 0, rng.randint(-5, 5)))
            else:
                a = rng.randint(-4, 4)
                b = rng.randint(-4, 4)
                if a == 0 and b == 0:
                    b = 1
                planes.append((a, b, rng.randint(-20, 20)))
    return planes


def run_tests() -> None:
    B = BOX

    # ---- 固定用例：三角形 x ≥ 0、y ≥ 0、x + y ≤ 10 ----
    planes = [(-1, 0, 0), (0, -1, 0), (1, 1, 10)]
    out = run(planes)
    assert out == "3\n0/1 0/1\n10/1 0/1\n0/1 10/1\n", out

    # ---- 固定用例：正方形 [0,10]² ----
    planes = [(1, 0, 10), (-1, 0, 0), (0, 1, 10), (0, -1, 0)]
    out = run(planes)
    assert out == ("4\n0/1 0/1\n10/1 0/1\n10/1 10/1\n0/1 10/1\n"), out

    # ---- 固定用例：退化成一個點 ----
    planes = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    assert run(planes) == "1\n0/1 0/1\n"

    # ---- 固定用例：空集（x ≤ 0 且 x ≥ 1）----
    assert run([(1, 0, 0), (-1, 0, -1)]) == "0\n"

    # ---- 固定用例：退化平面 ----
    assert run([(0, 0, -1)]) == "0\n"          # 0 ≤ -1，空集
    assert run([(0, 0, 5)]) == (f"4\n-{B}/1 -{B}/1\n{B}/1 -{B}/1\n"
                                  f"{B}/1 {B}/1\n-{B}/1 {B}/1\n")

    # ---- 固定用例：無界半平面 x ≤ 5 → 被盒子裁成四邊形 ----
    out = run([(1, 0, 5)])
    verts = parse_verts(out)
    assert len(verts) == 4, out
    assert sorted(verts) == sorted([
        (Fraction(5), Fraction(-B)), (Fraction(5), Fraction(B)),
        (Fraction(-B), Fraction(-B)), (Fraction(-B), Fraction(B))]), out

    # ---- 固定用例：長條 3 ≤ x ≤ 5 ----
    out = run([(1, 0, 5), (-1, 0, -3)])
    verts = parse_verts(out)
    assert sorted(verts) == sorted([
        (Fraction(3), Fraction(-B)), (Fraction(3), Fraction(B)),
        (Fraction(5), Fraction(-B)), (Fraction(5), Fraction(B))]), out

    # ---- 隨機對拍：候選枚舉 vs 多邊形裁剪 + 雙端棧 + 結構 / 成員檢查 ----
    rng = random.Random(20261023)
    for t in range(220):
        style = t % 4
        planes = gen_planes(rng, style)
        L, empty = normalize_planes(planes)
        if empty:
            assert run(planes) == "0\n"
            assert hpi_clip(planes) == []
            continue
        verts = parse_verts(run(planes))
        # 參照 1：Sutherland–Hodgman 精確裁剪（空 / 點 / 線段 / 多邊形全覆蓋）
        assert sorted(verts) == hpi_clip(planes), (t, style, planes, verts)
        check_structure(L, verts)
        check_membership(planes, verts, rng)
        # 參照 2：經典雙端棧版（退化 / 空集時回傳 None，直接跳過）
        dq = hpi_deque(L)
        if dq is not None:
            assert sorted(verts) == dq, (t, style, planes, verts, dq)
    # ---- 隨機對拍（較大規模，候選枚舉 vs 兩個參照）----
    for t in range(30):
        style = t % 4
        rng2 = random.Random(90000 + t)
        planes = gen_planes(rng2, style)
        # 放大到 40 個平面
        while len(planes) < 40:
            planes += gen_planes(rng2, style)
        planes = planes[:40]
        L, empty = normalize_planes(planes)
        if empty:
            assert hpi_clip(planes) == []
            continue
        verts = parse_verts(run(planes))
        assert sorted(verts) == hpi_clip(planes), (t, style)
        dq = hpi_deque(L)
        if dq is not None:
            assert sorted(verts) == dq, (t, style, len(verts), len(dq))

    # ---- 空輸入（沒有任何半平面）→ 整個包圍盒 ----
    assert run([]) == (f"4\n-{B}/1 -{B}/1\n{B}/1 -{B}/1\n"
                       f"{B}/1 {B}/1\n-{B}/1 {B}/1\n")


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出，避免 Windows 把 \n 換成 \r\n。"""
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
        return max(-CLAMP_ABS, min(CLAMP_ABS, v))

    def cl(v: int, lo: int, hi: int) -> int:
        return max(lo, min(hi, v))

    n = cl(next_int(0), 0, NMAX)
    planes: List[Plane] = []
    for _ in range(n):
        a = cl(next_int(0), -AMAX, AMAX)
        b = cl(next_int(0), -AMAX, AMAX)
        c = cl(next_int(0), -CMAX, CMAX)
        planes.append((a, b, c))
    emit(run(planes))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()

