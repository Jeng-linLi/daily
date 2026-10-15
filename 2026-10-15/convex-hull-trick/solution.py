"""凸包優化（Convex Hull Trick）與李超線段樹（Li Chao Segment Tree）

題意：
    給定 k 條直線 y = m·x + b、q 個查詢點 x，以及一個長度 n 的陣列 a。
    需要回答三件事：
      1. 對每個 x，求所有直線在 x 上的最小值 —— 用**李超線段樹** O(log C)，
         並與 O(k) 的暴力逐點驗證。
      2. 把直線依斜率排序後建**單調凸包**（deque），查詢點排序後 O(1) 均攤回答。
      3. 斜率優化 DP：把 a 切成若干段，每段成本 = (段內元素和)² + C，求最小總成本。
         dp[i] = min_{j<i} { dp[j] + (pre[i] - pre[j])² + C }，
         暴力 O(n²) 與李超線段樹 O(n log C) 互相驗證。

思路：
    ### 為什麼需要「凸包優化」
    一大類 DP 可以寫成
        dp[i] = min_{j<i} { dp[j] + f(j)·g(i) + h(i) + w(j) }
    把 j 視成一條直線：**斜率 m_j = f(j)、截距 b_j = dp[j] + w(j)**，
    而 i 的轉移就是在 x = g(i) 上取所有直線的最小值。
    於是把「O(n) 找最小值」換成「一堆直線求 min」，就有了優化的空間。

    ### 李超線段樹：把「直線」當成區間上的懶標記
    定義域 x ∈ [L, R] 建立在線段樹上，每個節點存一條「在該區間中點最優」的直線。
    插入新直線 nw 時，跟節點當前的直線 cur 比中點：
      - 若 nw 在中點更優 → 交換，nw 變成「被擠掉的那一條」往下傳；
      - 若 nw 在左端點更優 → 它只可能在左半區間勝出，往左遞迴；
      - 否則若 nw 在右端點更優 → 往右遞迴；
      - 兩端都不優 → nw 在整個區間都被 cur 壓住，直接丟掉。
    正確性來自「兩條直線最多交一次」：在一個區間裡，若 nw 在兩端都輸給 cur，
    則它在區間內部也不可能贏。每個節點只存一條線，插入 O(log C)、查詢 O(log C)。
    優點是**完全不要求斜率或查詢點單調**，缺點是需要離散（或座標壓縮）的定義域。

    ### 單調凸包（dq-CHT）：斜率與查詢都單調時可以做到 O(1)
    把直線依**斜率遞減**排列後建下凸包：斜率最大者在 x → -∞ 時最優，
    斜率最小者在 x → +∞ 時最優，所以查詢點 x 遞增時，最佳直線在 hull 上的
    位置只會單調往前走，用一根指標就夠，全部查詢 O(q + k)。
      - 加入：若最後兩條與新線的交點不滿足「交點遞增」，彈掉最後一條；
      - 查詢：若下一條比當前好，指標就前進（指標不回頭）。
    判斷「中間那條是不是多餘的」用整數交叉相乘，完全不碰浮點數
    （斜率遞減：m₁ > m₂ > m₃）：
        交點順序亂掉 ⟺ (b₂-b₁)(m₂-m₃) ≥ (b₃-b₂)(m₁-m₂)

    ### 兩者的取捨
      李超線段樹：任意順序、任意斜率，O(log C)，程式稍長，定義域要離散。
      單調凸包　：要求斜率單調 + 查詢單調，O(1) 均攤，常數極小。
      實務上：DP 轉移順序天然單調時優先單調凸包；否則用李超。

應用場景：
    - 斜率優化 DP：序列分段代價（本題）、土地購買、任務排程、IOI 的「印刷廠」
    - 計算幾何：一堆直線求上下包絡、半平面交的前置處理
    - 圖論最短路的「勢能」轉換、費用流的凸費用邊建模
    - 機器學習：分段線性函數（ReLU 組合）的最小值查詢

複雜度（k 條線、q 次查詢、定義域寬度 C、DP 長度 n）：
    李超線段樹插入 / 查詢   O(log C)，空間 O(k log C)
    單調凸包建構 O(k log k)（排序）/ O(k)（已排序）；全部查詢 O(q + k)
    斜率優化 DP             O(n log C)（李超）或 O(n)（單調凸包）
"""

import re
import sys
from typing import List, Optional, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

INF = 4 * 10 ** 18          # 整數哨兵（不要用 float，會污染整數輸出）

XLO, XHI = -10 ** 7, 10 ** 7   # 李超線段樹的定義域（離散整數座標）

MAXK = 200                  # 直線數上限
MAXQ = 200                  # 查詢數上限
MAXN = 150                  # DP 陣列長度上限
MLIM = 10 ** 6              # 斜率截斷
BLIM = 10 ** 9              # 截距截斷
XLIM = 10 ** 6              # 查詢點截斷
ALIM = 10 ** 4              # DP 陣列元素截斷
COST = 10                   # 每段固定成本 C

Line = Tuple[int, int]


def ev(ln: Line, x: int) -> int:
    return ln[0] * x + ln[1]


# ---------------------------------------------------------------- 李超線段樹


class LiChao:
    """動態開點的李超線段樹，維護 x ∈ [xlo, xhi] 上所有直線的最小值。"""

    __slots__ = ("xlo", "xhi", "line", "lc", "rc")

    def __init__(self, xlo: int, xhi: int) -> None:
        self.xlo = xlo
        self.xhi = xhi
        self.line: List[Optional[Line]] = [None]
        self.lc: List[int] = [-1]
        self.rc: List[int] = [-1]

    def _node(self) -> int:
        self.line.append(None)
        self.lc.append(-1)
        self.rc.append(-1)
        return len(self.line) - 1

    def add(self, nw: Line) -> None:
        self._add(nw, 0, self.xlo, self.xhi)

    def _add(self, nw: Line, v: int, l: int, r: int) -> None:
        cur = self.line[v]
        if cur is None:
            self.line[v] = nw
            return
        # 用 l + (r - l) // 2 保證兩個語言在負座標下算出同一個中點
        mid = l + (r - l) // 2
        if ev(nw, mid) < ev(cur, mid):
            nw, cur = cur, nw
            self.line[v] = cur
        if l == r:
            return
        if ev(nw, l) < ev(cur, l):
            if self.lc[v] == -1:
                self.lc[v] = self._node()
            self._add(nw, self.lc[v], l, mid)
        elif ev(nw, r) < ev(cur, r):
            if self.rc[v] == -1:
                self.rc[v] = self._node()
            self._add(nw, self.rc[v], mid + 1, r)
        # 兩端都不優 → nw 在整個區間都被 cur 壓住（兩線只交一次），丟掉

    def query(self, x: int) -> int:
        return self._query(x, 0, self.xlo, self.xhi)

    def _query(self, x: int, v: int, l: int, r: int) -> int:
        cur = self.line[v]
        res = INF if cur is None else ev(cur, x)
        if l == r:
            return res
        mid = l + (r - l) // 2
        if x <= mid:
            c = self.lc[v]
            if c != -1:
                t = self._query(x, c, l, mid)
                if t < res:
                    res = t
        else:
            c = self.rc[v]
            if c != -1:
                t = self._query(x, c, mid + 1, r)
                if t < res:
                    res = t
        return res


# ---------------------------------------------------------------- 單調凸包


def bad(l1: Line, l2: Line, l3: Line) -> bool:
    """斜率嚴格**遞減**時，判斷 l2 是否被 l1 與 l3 夾成多餘（全整數運算）。"""
    # 交點(l1,l2) 應小於 交點(l2,l3)；用交叉相乘避免浮點
    return (l2[1] - l1[1]) * (l2[0] - l3[0]) >= (l3[1] - l2[1]) * (l1[0] - l2[0])


def cht_build(lines: List[Line]) -> List[Line]:
    """依斜率**遞減**建下凸包；同斜率只保留截距最小的那條。"""
    # 排序鍵 (-m, b)：斜率由大到小，同斜率時截距小的排前面（直接捨棄後面的）
    ls = sorted(lines, key=lambda t: (-t[0], t[1]))
    hull: List[Line] = []
    for ln in ls:
        if hull and hull[-1][0] == ln[0]:
            continue          # 同斜率且截距不是最小的，直接丟掉
        while len(hull) >= 2 and bad(hull[-2], hull[-1], ln):
            hull.pop()
        hull.append(ln)
    return hull


def cht_query(hull: List[Line], xs: List[int]) -> List[int]:
    """xs 必須**升序**；指標只前進，全部 O(q + |hull|)。"""
    res: List[int] = []
    if not hull:
        return [INF] * len(xs)
    j = 0
    for x in xs:
        while j + 1 < len(hull) and ev(hull[j], x) >= ev(hull[j + 1], x):
            j += 1
        res.append(ev(hull[j], x))
    return res


# ---------------------------------------------------------------- 斜率優化 DP


def dp_li_chao(a: List[int]) -> List[int]:
    """dp[i] = min_{j<i} { dp[j] + (pre[i]-pre[j])² + C }，回傳 dp[1..n]。"""
    n = len(a)
    pre = [0] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] + a[i]
    lc = LiChao(XLO, XHI)
    dp = [0] * (n + 1)
    # j = 0 的線：斜率 -2·pre[0] = 0，截距 dp[0] + pre[0]² + C = C
    lc.add((0, COST))
    for i in range(1, n + 1):
        dp[i] = pre[i] * pre[i] + lc.query(pre[i])
        lc.add((-2 * pre[i], dp[i] + pre[i] * pre[i] + COST))
    return dp[1:]


def dp_brute(a: List[int]) -> List[int]:
    n = len(a)
    pre = [0] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] + a[i]
    dp = [INF] * (n + 1)
    dp[0] = 0
    for i in range(1, n + 1):
        best = INF
        for j in range(i):
            v = dp[j] + (pre[i] - pre[j]) ** 2 + COST
            if v < best:
                best = v
        dp[i] = best
    return dp[1:]


# ---------------------------------------------------------------- 暴力對拍


def brute_min(lines: List[Line], xs: List[int]) -> List[int]:
    out: List[int] = []
    for x in xs:
        best = INF
        for ln in lines:
            v = ev(ln, x)
            if v < best:
                best = v
        out.append(best)
    return out


# ---------------------------------------------------------------- 驅動


def join(nums: List[int]) -> str:
    return " ".join(str(v) for v in nums)


def run(lines: List[Line], xs: List[int], a: List[int]) -> str:
    lc = LiChao(XLO, XHI)
    for ln in lines:
        lc.add(ln)
    a1 = [lc.query(x) for x in xs]
    a2 = brute_min(lines, xs)

    hull = cht_build(lines)
    sx = sorted(xs)
    a3 = cht_query(hull, sx)
    a4 = brute_min(lines, sx)

    d1 = dp_li_chao(a)
    d2 = dp_brute(a)

    ok = 1 if (a1 == a2 and a3 == a4 and d1 == d2) else 0
    return "\n".join([join(a1), join(a2), join(a3), join(d1), join(d2), str(ok)])


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    import random

    # ---- 固定用例 ----
    lines = [(1, 0), (-1, 5), (0, 3)]
    xs = [-5, 0, 5, 10]
    assert run(lines, xs, [1, 2, 3]).endswith("\n1")

    # 只有一條線：最小值永遠是它自己
    assert run([(2, 7)], [0, 1, -1], []).split("\n")[0] == "7 9 5"

    # 空輸入
    assert run([], [], []) == "\n\n\n\n\n1"

    # 兩條過原點的線 y = x 與 y = -x：最小值恆為 -|x|
    assert cht_query(cht_build([(1, 0), (-1, 0)]), [-3, -1, 0, 1, 3]) == [-3, -1, 0, -1, -3]

    # DP：李超版必須與暴力版一致
    assert dp_brute([1, 2, 3]) == dp_li_chao([1, 2, 3])

    # ---- 隨機對拍 ----
    rng = random.Random(20261015)
    for _ in range(300):
        k = rng.randint(0, 30)
        lines = [(rng.randint(-20, 20), rng.randint(-50, 50)) for _ in range(k)]
        q = rng.randint(0, 30)
        xs = [rng.randint(-30, 30) for _ in range(q)]
        n = rng.randint(0, 25)
        a = [rng.randint(-30, 30) for _ in range(n)]
        assert run(lines, xs, a).endswith("\n1")

    # ---- 大座標只用李超（不走 O(k·q) 暴力）----
    big = [(rng.randint(-MLIM, MLIM), rng.randint(-BLIM, BLIM)) for _ in range(MAXK)]
    lc = LiChao(XLO, XHI)
    for ln in big:
        lc.add(ln)
    for x in (XLO, -1, 0, 1, XHI):
        assert lc.query(x) == brute_min(big, [x])[0]
    # 單調凸包在大樣本上與李超一致（查詢點升序）
    sx = list(range(-XLIM, XLIM + 1, 20000))
    hull = cht_build(big)
    assert cht_query(hull, sx) == [lc.query(x) for x in sx]


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出（latin-1），避免 Windows 把 \\n 換成 \\r\\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。"""
    tokens = raw.split()
    pos = 0

    def read_int() -> Optional[int]:
        nonlocal pos
        if pos < len(tokens) and NUM_RE.match(tokens[pos]):
            v = int(tokens[pos])
            pos += 1
            return v
        return None

    def clamp(v: Optional[int], lim: int) -> int:
        if v is None:
            return 0
        return max(-lim, min(lim, v))

    def clamp_to(v: Optional[int], hi: int) -> int:
        if v is None:
            return 0
        return max(0, min(hi, v))

    k = clamp_to(read_int(), MAXK)
    lines: List[Line] = []
    for _ in range(k):
        m = read_int()
        b = read_int()
        if m is None or b is None:
            break
        lines.append((clamp(m, MLIM), clamp(b, BLIM)))

    q = clamp_to(read_int(), MAXQ)
    xs: List[int] = []
    for _ in range(q):
        v = read_int()
        if v is None:
            break
        xs.append(clamp(v, XLIM))

    n = clamp_to(read_int(), MAXN)
    a: List[int] = []
    for _ in range(n):
        v = read_int()
        if v is None:
            break
        a.append(clamp(v, ALIM))

    emit(run(lines, xs, a) + "\n")


def main() -> None:
    sys.setrecursionlimit(1000000)
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
