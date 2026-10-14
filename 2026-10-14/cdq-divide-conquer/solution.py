"""CDQ 分治：逆序對計數 / 三維偏序（陌上花開）/ 分治優化 DP

題意：
    給定三個長度為 n 的整數陣列 a、b、c（第 i 個元素組成的三元組為 (a_i, b_i, c_i)）
    以及一個權重陣列 w。輸出四項結果：
      1. 陣列 b 的**逆序對個數**；
      2. 對每個 i，三元組偏序計數 f(i) = #{ j ≠ i : a_j ≤ a_i, b_j ≤ b_i, c_j ≤ c_i }；
      3. 分治優化 DP：dp[i] = w_i + max{ dp[j] : j < i 且 a_j ≤ a_i }
         （沒有合法的 j 時 dp[i] = w_i）；
      4. dp 的最大值（n = 0 時記 0）。
    以上全部用 CDQ 分治完成，並與 O(n²) 暴力對拍。

思路：
    ### CDQ 分治在做什麼
    普通的「分治」（如合併排序）是：把問題切成兩半、各自遞迴、**合併**。
    CDQ 分治多了一步：**統計「跨左右兩半」的貢獻**。
    也就是在遞迴解決左半與右半之後，額外回答一個問題：
        「有多少對 (u, v) 滿足 u 在左半、v 在右半，且 (u, v) 滿足某個條件？」
    只要每一對 (u, v) 在遞迴樹上恰好在某一個節點被統計一次
    （就是在它們「分開的那一層」），全部跨半貢獻加起來就是全體配對的答案。
    這就是把一個 O(n²) 的配對統計降到 O(n log n)（或 O(n log² n)）的關鍵。

    ### ① 逆序對：一維 CDQ
    求 #{(i, j) : i < j 且 b_i > b_j}。
    分治時左右兩半都已經「依 b 排序」，用兩個指標把左右歸併：
    右半的元素 b_j 在左半指標停下來時，左半剩下的每一個都比它大 → 一次加一整段。
    這其實就是「合併排序計逆序對」，只是從 CDQ 的視角看它。

    ### ② 三維偏序：CDQ + 樹狀陣列
    先依 (a, b, c) 字典序排序。因為已依 a 升序，對任意「左半位置 u < 右半位置 v」
    必有 a_u ≤ a_v —— a 這一維**自動滿足**，只剩下 b、c 兩維要處理。
    CDQ 的不變量是「遞迴返回後，區段內依 b 排序」，於是可以：
      - 兩個指標掃描（左半依 b、右半依 b）；
      - 左半走過的元素把它的 c 加到樹狀陣列（c 需離散化）；
      - 右半每個 v 查詢「c ≤ c_v 的個數」，就是左半對它的貢獻。
    這樣跨半統計是 O(n log n)，遞迴 log n 層 → 總共 O(n log² n)。

    **相同三元組要特別處理**：CDQ 只數到「排序後位置 u < v」的配對，
    但題意要求 j ≠ i 的全部配對。若 (a,b,c) 完全相同的一組在排序後佔了
    位置 p .. p+s-1，那麼對位置 v 的元素，組內還有 s-1-(v-p) 個排在它後面、
    沒被 CDQ 數到的同伴。補上這一項就與暴力完全一致。

    ### ③ 分治優化 DP
    遞推 dp[i] = w_i + max{ dp[j] : j < i 且 a_j ≤ a_i }。
    暴力是 O(n²)：對每個 i 掃所有 j。CDQ 的做法是把它看成
    「二維偏序上的最大值」——第一維是下標、第二維是 a。
    CDQ(l, r) 的順序（這是本題最容易寫錯的地方）：
      1. 先遞迴左半 cdq(l, m)：左半的 dp 全部算完，且左半**依 a 排序**；
      2. 跨半更新：左半依 a 排序後可以預處理「前綴最大值」，
         右半每個 i 用二分找到「a ≤ a_i 的最後一個位置」，
         一次 O(log n) 拿到 max{dp[j]}，更新 dp[i]；
         （右半此時仍是原始下標順序，正好符合「j < i」的要求）
      3. 再遞迴右半 cdq(m, r)：右半用左半更新後的 dp 繼續算；
      4. 歸併，讓整個區段依 a 排序，維持不變量。
    第 2 步必須在第 3 步之前——右半的 dp 需要左半的貢獻。

應用場景：
    - 逆序對：排序代價、推薦系統的「偏好一致性」度量（Kendall tau 距離）
    - 三維偏序：多維度支配關係統計（「各方面都不比我差的人有多少」）、
      資料庫的 skyline / 支配計數查詢、遊戲戰力排行榜
    - 分治優化 DP：一維 DP 加上一個偏序約束時（LIS 帶權、嵌套矩形、
      多條件的最長鏈），把 O(n²) 降到 O(n log n)
    - 進階：CDQ 分治 + FFT（多項式卷積）、CDQ 套 CDQ（四維偏序）

複雜度（n 為元素個數）：
    逆序對（一維 CDQ）   O(n log n)   空間 O(n)
    三維偏序（CDQ + BIT） O(n log² n)  空間 O(n)
    分治優化 DP          O(n log n)   空間 O(n)
    暴力對拍             O(n²)

輸入格式（stdin，全部以空白分隔的整數）：
    n
    a_1 ... a_n
    b_1 ... b_n
    c_1 ... c_n
    w_1 ... w_n
    非數字 token 一律視為「輸入到此為止」；讀不滿的那一列直接丟棄。
    n 截斷到 [0, 3000]；數值截斷到 [-10^12, 10^12]。
輸出格式（stdout，共 4 行）：
    第 1 行  逆序對個數（陣列 b）
    第 2 行  f(1) ... f(n)（三元組偏序計數，原輸入順序）
    第 3 行  dp_1 ... dp_n
    第 4 行  max(dp)，n = 0 時為 0
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import List, Optional, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 3000
VMAX = 10 ** 12
NEG_INF = -(10 ** 18)


# ---------------------------------------------------------------- 樹狀陣列


class Fenwick:
    """1-indexed 樹狀陣列（前綴和）。"""

    __slots__ = ("n", "t")

    def __init__(self, n: int) -> None:
        self.n = n
        self.t = [0] * (n + 1)

    def add(self, i: int, delta: int) -> None:
        while i <= self.n:
            self.t[i] += delta
            i += i & (-i)

    def sum(self, i: int) -> int:
        """前綴和 [1, i]；i ≤ 0 時為 0，i > n 時為全部。"""
        if i > self.n:
            i = self.n
        res = 0
        while i > 0:
            res += self.t[i]
            i -= i & (-i)
        return res

    def reset(self, i: int) -> None:
        """把位置 i 歸零（用過的位置逐一同 add(-值) 太慢，這裡直接加 -1 計數）。"""
        self.add(i, -1)


# ---------------------------------------------------------------- ① 逆序對


def inversion_count(arr: Sequence[int]) -> int:
    """CDQ 分治數逆序對：#{(i,j) : i<j 且 arr[i] > arr[j]}。"""
    a = list(arr)
    tmp = [0] * len(a)

    def cdq(l: int, r: int) -> int:
        if r - l <= 1:
            return 0
        m = (l + r) // 2
        cnt = cdq(l, m) + cdq(m, r)
        i, j, k = l, m, l
        while i < m and j < r:
            if a[i] <= a[j]:
                tmp[k] = a[i]
                i += 1
            else:
                cnt += m - i          # 左半剩下的全部 > a[j]
                tmp[k] = a[j]
                j += 1
            k += 1
        while i < m:
            tmp[k] = a[i]
            i += 1
            k += 1
        while j < r:
            tmp[k] = a[j]
            j += 1
            k += 1
        a[l:r] = tmp[l:r]
        return cnt

    return cdq(0, len(a))


# ---------------------------------------------------------------- ② 三維偏序


def compress(vals: Sequence[int]) -> Tuple[List[int], int]:
    """離散化成 1-indexed 的名次（相同值同名次）。"""
    uniq = sorted(set(vals))
    return [bisect_left(uniq, v) + 1 for v in vals], len(uniq)


def bisect_left(arr: Sequence[int], x: int) -> int:
    lo, hi = 0, len(arr)
    while lo < hi:
        mid = (lo + hi) // 2
        if arr[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def partial_order_3d(a: Sequence[int], b: Sequence[int],
                     c: Sequence[int]) -> List[int]:
    """對每個原始下標 i，回傳 #{ j ≠ i : a_j≤a_i, b_j≤b_i, c_j≤c_i }。"""
    n = len(a)
    if n == 0:
        return []
    cr, m = compress(c)

    # 依 (a, b, c) 排序，記住原始下標
    order = sorted(range(n), key=lambda i: (a[i], b[i], c[i]))
    sa = [a[i] for i in order]
    sb = [b[i] for i in order]
    sc = [cr[i] for i in order]
    g = [0] * n          # g[v] = 排序後位置 < v 且支配 v 的個數

    # 相同三元組的分組必須在 CDQ 打亂順序「之前」算好，
    # 而且要以**原始下標**為鍵存起來（CDQ 返回時位置 v 已經換成別的元素了）
    pos_orig = [0] * n        # 原始下標 i → 排序後的位置
    grp_start = [0] * n       # 原始下標 i → 同組的起始位置
    grp_size = [0] * n        # 原始下標 i → 同組的大小
    for v in range(n):
        pos_orig[order[v]] = v
    p = 0
    while p < n:
        q = p + 1
        while q < n and (sa[q], sb[q], sc[q]) == (sa[p], sb[p], sc[p]):
            q += 1
        for v in range(p, q):
            grp_start[order[v]] = p
            grp_size[order[v]] = q - p
        p = q

    bit = Fenwick(m)
    tmp_a = [0] * n
    tmp_b = [0] * n
    tmp_c = [0] * n
    tmp_g = [0] * n
    tmp_o = [0] * n

    def cdq(l: int, r: int) -> None:
        if r - l <= 1:
            return
        m = (l + r) // 2
        cdq(l, m)
        cdq(m, r)
        # 跨半統計：左右都已依 b 排序
        i = l
        for v in range(m, r):
            while i < m and sb[i] <= sb[v]:
                bit.add(sc[i], 1)
                i += 1
            g[v] += bit.sum(sc[v])
        for t in range(l, i):
            bit.add(sc[t], -1)
        # 依 b 歸併（相同 b 時依 c，保證確定性）
        p, q, k = l, m, l
        while p < m and q < r:
            if (sb[p], sc[p]) <= (sb[q], sc[q]):
                tmp_a[k], tmp_b[k], tmp_c[k], tmp_o[k] = sa[p], sb[p], sc[p], order[p]
                tmp_g[k] = g[p]
                p += 1
            else:
                tmp_a[k], tmp_b[k], tmp_c[k], tmp_o[k] = sa[q], sb[q], sc[q], order[q]
                tmp_g[k] = g[q]
                q += 1
            k += 1
        while p < m:
            tmp_a[k], tmp_b[k], tmp_c[k], tmp_o[k] = sa[p], sb[p], sc[p], order[p]
            tmp_g[k] = g[p]
            p += 1
            k += 1
        while q < r:
            tmp_a[k], tmp_b[k], tmp_c[k], tmp_o[k] = sa[q], sb[q], sc[q], order[q]
            tmp_g[k] = g[q]
            q += 1
            k += 1
        for k in range(l, r):
            sa[k], sb[k], sc[k], order[k], g[k] = (
                tmp_a[k], tmp_b[k], tmp_c[k], tmp_o[k], tmp_g[k])

    cdq(0, n)

    # 相同三元組：補上「排在它後面、沒被 CDQ 數到的組內同伴」
    res = [0] * n
    for v in range(n):
        i = order[v]
        res[i] = g[v] + (grp_size[i] - 1 - (pos_orig[i] - grp_start[i]))
    return res


# ---------------------------------------------------------------- ③ 分治優化 DP


def cdq_dp(a: Sequence[int], w: Sequence[int]) -> List[int]:
    """dp[i] = w[i] + max{ dp[j] : j < i 且 a[j] ≤ a[i] }，無合法 j 時 dp[i] = w[i]。"""
    n = len(a)
    if n == 0:
        return []
    idx = list(range(n))          # 區段內的下標（會在歸併時重排）
    aval = list(a)
    # acc[v] = 目前為止「合法 j 的 dp[j]」的最大值；完全沒有合法 j 時為 NEG_INF
    acc = [NEG_INF] * n
    # 「沒有合法 j」與「有一個 dp[j] = 0 的合法 j」語義不同，
    # 所以不能讓 w[i] 直接參與 max —— 必須用 acc 是否為 NEG_INF 來區分
    dval = [0] * n
    tmp_i = [0] * n
    tmp_a = [0] * n
    tmp_c = [NEG_INF] * n

    def cdq(l: int, r: int) -> None:
        if r - l <= 1:
            return
        m = (l + r) // 2
        # 1) 左半先算完（返回時左半依 a 排序），左半的 dp 已是最終值
        cdq(l, m)
        for t in range(l, m):
            gain = acc[t]
            dval[t] = w[idx[t]] if gain == NEG_INF else w[idx[t]] + gain
        # 2) 跨半更新：左半依 a 排序 → 前綴最大值 + 二分
        pref = [NEG_INF] * (m - l)
        best = NEG_INF
        for t in range(l, m):
            if dval[t] > best:
                best = dval[t]
            pref[t - l] = best
        keys = [aval[t] for t in range(l, m)]
        for v in range(m, r):
            # 最後一個 aval ≤ aval[v] 的位置
            lo, hi = 0, m - l
            while lo < hi:
                mid = (lo + hi) // 2
                if keys[mid] <= aval[v]:
                    lo = mid + 1
                else:
                    hi = mid
            if lo > 0 and pref[lo - 1] > acc[v]:
                acc[v] = pref[lo - 1]
        # 3) 右半繼續算（它已經拿到左半的貢獻）
        cdq(m, r)
        # 4) 依 a 歸併（相同 a 時依下標，保證確定性）
        p, q, k = l, m, l
        while p < m and q < r:
            if (aval[p], idx[p]) <= (aval[q], idx[q]):
                tmp_i[k], tmp_a[k], tmp_c[k] = idx[p], aval[p], acc[p]
                p += 1
            else:
                tmp_i[k], tmp_a[k], tmp_c[k] = idx[q], aval[q], acc[q]
                q += 1
            k += 1
        while p < m:
            tmp_i[k], tmp_a[k], tmp_c[k] = idx[p], aval[p], acc[p]
            p += 1
            k += 1
        while q < r:
            tmp_i[k], tmp_a[k], tmp_c[k] = idx[q], aval[q], acc[q]
            q += 1
            k += 1
        for k in range(l, r):
            idx[k], aval[k], acc[k] = tmp_i[k], tmp_a[k], tmp_c[k]

    cdq(0, n)

    out = [0] * n
    for k in range(n):
        gain = acc[k]
        out[idx[k]] = w[idx[k]] if gain == NEG_INF else w[idx[k]] + gain
    return out


# ---------------------------------------------------------------- 主流程


def run(a: Sequence[int], b: Sequence[int], c: Sequence[int],
        w: Sequence[int]) -> str:
    inv = inversion_count(b)
    f = partial_order_3d(a, b, c)
    dp = cdq_dp(a, w)
    best = max(dp) if dp else 0
    lines = [
        str(inv),
        " ".join(str(v) for v in f),
        " ".join(str(v) for v in dp),
        str(best),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 對拍用暴力


def bf_inversion(b: Sequence[int]) -> int:
    n = len(b)
    return sum(1 for i in range(n) for j in range(i + 1, n) if b[i] > b[j])


def bf_partial_order(a: Sequence[int], b: Sequence[int], c: Sequence[int]) -> List[int]:
    n = len(a)
    res = [0] * n
    for i in range(n):
        cnt = 0
        for j in range(n):
            if j != i and a[j] <= a[i] and b[j] <= b[i] and c[j] <= c[i]:
                cnt += 1
        res[i] = cnt
    return res


def bf_dp(a: Sequence[int], w: Sequence[int]) -> List[int]:
    n = len(a)
    dp = [0] * n
    for i in range(n):
        best = None
        for j in range(i):
            if a[j] <= a[i] and (best is None or dp[j] > best):
                best = dp[j]
        dp[i] = w[i] if best is None else w[i] + best
    return dp


def bf_run(a: Sequence[int], b: Sequence[int], c: Sequence[int],
           w: Sequence[int]) -> str:
    dp = bf_dp(a, w)
    best = max(dp) if dp else 0
    lines = [
        str(bf_inversion(b)),
        " ".join(str(v) for v in bf_partial_order(a, b, c)),
        " ".join(str(v) for v in dp),
        str(best),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check(a: Sequence[int], b: Sequence[int], c: Sequence[int],
          w: Sequence[int]) -> None:
    got = run(a, b, c, w)
    want = bf_run(a, b, c, w)
    assert got == want, (list(a), list(b), list(c), list(w), got, want)


def run_tests() -> None:
    # ---- 固定用例 ----
    assert run([], [], [], []) == "0\n\n\n0"
    assert run([1], [5], [7], [3]) == "0\n0\n3\n3"
    assert run([1, 2], [2, 1], [1, 2], [1, 2]) == "1\n0 0\n1 3\n3"
    assert run([1, 1], [1, 1], [1, 1], [5, 5]) == "0\n1 1\n5 10\n10"

    check([], [], [], [])
    check([1], [1], [1], [1])
    check([1, 2], [2, 1], [1, 2], [1, 2])
    check([1, 1], [1, 1], [1, 1], [5, 5])
    check([3, 1, 2], [1, 2, 3], [2, 2, 2], [4, -1, 7])
    check([1, 2, 3], [3, 2, 1], [1, 1, 1], [1, 1, 1])

    rng = random.Random(20261014)

    # ---- 隨機對拍：小規模全值域 ----
    for _ in range(320):
        n = rng.randint(0, 12)
        a = [rng.randint(0, 3) for _ in range(n)]
        b = [rng.randint(0, 3) for _ in range(n)]
        c = [rng.randint(0, 3) for _ in range(n)]
        w = [rng.randint(-5, 5) for _ in range(n)]
        check(a, b, c, w)

    # ---- 隨機對拍：大值域（離散化路徑）----
    for _ in range(120):
        n = rng.randint(0, 40)
        a = [rng.randint(-50, 50) for _ in range(n)]
        b = [rng.randint(-50, 50) for _ in range(n)]
        c = [rng.randint(-50, 50) for _ in range(n)]
        w = [rng.randint(-30, 30) for _ in range(n)]
        check(a, b, c, w)

    # ---- 隨機對拍：大量重複值（相同三元組分組邏輯）----
    for _ in range(160):
        n = rng.randint(1, 30)
        a = [rng.randint(0, 1) for _ in range(n)]
        b = [rng.randint(0, 1) for _ in range(n)]
        c = [rng.randint(0, 1) for _ in range(n)]
        w = [rng.randint(-4, 4) for _ in range(n)]
        check(a, b, c, w)

    # ---- 遞增 / 遞減極端 ----
    for n in (1, 2, 3, 17, 64):
        check(list(range(n)), list(range(n)), list(range(n)), [1] * n)
        check(list(range(n)), list(range(n - 1, -1, -1)), [0] * n,
              list(range(1, n + 1)))

    # ---- 中量級效能（不跑 O(n²) 暴力，只確認不崩且自洽）----
    n = 3000
    a = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(n)]
    b = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(n)]
    c = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(n)]
    w = [rng.randint(-100, 100) for _ in range(n)]
    out = run(a, b, c, w).split("\n")
    assert len(out) == 4
    assert len(out[1].split()) == n and len(out[2].split()) == n
    assert int(out[0]) == bf_inversion(b)
    assert [int(v) for v in out[2].split()] == bf_dp(a, w)


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    tokens = raw.split()
    pos = 0

    def read_int() -> Optional[int]:
        nonlocal pos
        if pos < len(tokens) and NUM_RE.match(tokens[pos]):
            v = int(tokens[pos])
            pos += 1
            return v
        return None

    def clamp_val(v: Optional[int]) -> int:
        if v is None:
            return 0
        return max(-VMAX, min(VMAX, v))

    n = 0
    v = read_int()
    if v is not None:
        n = max(0, min(MAXN, v))

    cols: List[List[int]] = []
    for _ in range(4):
        col: List[int] = []
        for _ in range(n):
            x = read_int()
            if x is None:
                break
            col.append(clamp_val(x))
        cols.append(col)

    # 任何一列讀不滿 → 統統截斷到最短長度
    m = min(len(x) for x in cols) if cols else 0
    a, b, c, w = (cols[0][:m], cols[1][:m], cols[2][:m], cols[3][:m])

    emit(run(a, b, c, w) + "\n")


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
