"""莫隊算法：區間查詢（普通莫隊 + 奇偶優化）與帶修莫隊

題意：
    維護一個整數陣列 a（長度 n）與兩組操作：
      A. q 個**靜態區間查詢** [l, r]：輸出該區間內
           ① 不同數的個數（distinct）
           ② 區間和（sum）
           ③ 眾數的出現次數（mode，出現最多次的那個數字出現了幾次）
      B. m 個**帶修改的操作**：`1 l r` 是查詢（在「執行到這裡為止的陣列」上回答），
           `2 p x` 是把 a[p] 改成 x。對每個查詢輸出同樣的三個量。
    A 用普通莫隊（塊排序 + 奇偶優化），B 用帶修莫隊（三維排序），
    並與「每次 O(n) 從頭掃一遍」的暴力對拍。

思路：
    ### 莫隊的核心：把「區間查詢」變成「指標在陣列上走路」
    如果我知道區間 [L, R] 的答案，那麼只要能 O(1)（或 O(log n)）地
      - 在右邊加一個元素 / 在右邊刪一個元素
      - 在左邊加一個元素 / 在左邊刪一個元素
    就可以 O(|Δ|) 地移動到任意相鄰區間。本題三個統計量都滿足這個條件：
      - `sum`      ：add 就 += v，remove 就 -= v；
      - `distinct` ：用 cnt[v] 記出現次數，cnt 從 0 變 1 時 +1、從 1 變 0 時 -1；
      - `mode`     ：再用一層「出現次數的出現次數」freqCnt[c] = 有多少個不同的值
                     恰好出現 c 次。加入時 cnt[v] 增加，若超過當前最大值就更新；
                     刪除時若最大值那一層被清空就把最大值往下退
                     （`while maxc > 0 and freqCnt[maxc] == 0: maxc -= 1`）。
                     這是「眾數出現次數」的經典 O(1) 維護法。

    ### 排序策略決定總移動量
    把所有查詢排序。最樸素是按 l 排序，r 可能每次都從頭跑到尾 → O(n·q)。
    莫隊的做法是**分塊**：塊長 B ≈ √n，排序鍵為
        (l 所在塊, r)     且奇數塊的 r 用**降序**（奇偶優化）
    這樣：
      - 同一塊內 r 單調移動 → 每塊 r 總共走 O(n)，共 O(√n · n) = O(n^1.5)；
      - 換塊時 r 最多跳回 n，共 O(√n) 次 → O(n^1.5)；
      - l 在塊內每次最多走 B → 每次查詢 O(√n)，共 O(q·√n)。
    奇偶優化讓 r 在塊與塊之間「來回掃」而不是每次都跳回起點，常數減半。

    ### 帶修莫隊：多一維「時間」
    現在陣列還會被修改。把每個查詢看成三維點 (l 塊, r 塊, t)，
    其中 t 是「這個查詢之前發生過幾次修改」。塊長取 B ≈ n^(2/3)，
    排序鍵 (l 塊, r 塊, t)。移動時除了動 L、R，還要把時間從 curT 推向 t：
      - 時間往前：套用修改（若被改的位置在當前區間內，要先 remove 舊值、add 新值）；
      - 時間往後：撤銷修改（同理）。
    複雜度 O(n^(5/3))。

    ### 為什麼兩個語言能逐位元組一致
    塊長用**整數立方根**（二分找最大的 b 使 b³ ≤ n²）算，不用浮點；
    排序鍵是純整數且加上查詢編號做最終 tie-break，因此排序結果唯一。
    而答案與處理順序無關，所以兩邊必然一致。

應用場景：
    - 靜態陣列的大量區間查詢：顏色計數、區間眾數、區間不同數、
      「區間內有多少對 (i,j) 滿足 a_i = a_j」、區間 mex
    - 帶修改版本：日誌回放分析、資料庫在「某個時間點」的區間統計
    - 變體：樹上莫隊（歐拉序把路徑轉成區間）、回滾莫隊（只加不刪，
      處理「區間最大值 / 需要 O(1) 撤銷」的問題）、二維莫隊

複雜度（n 為陣列長度，q 為查詢數）：
    普通莫隊（塊長 √n）    O((n + q)·√n)
    帶修莫隊（塊長 n^(2/3)） O(n^(5/3))
    暴力對拍                O(n·q)
    空間                    O(n + q)

輸入格式（stdin，全部以空白分隔的整數）：
    n
    a_1 ... a_n            初始陣列
    q                      靜態查詢個數
    l_1 r_1                q 組，1-indexed 閉區間
    ...
    m                      帶修改操作數
    接下來 m 組，每組 3 個整數：
        1 l r              查詢 [l, r]
        2 p x              把 a[p] 改成 x（1-indexed）
    非數字 token 一律視為「輸入到此為止」；讀不滿的那一組直接丟棄。
    n、q、m 截斷到 [0, 3000]；數值截斷到 [-10^9, 10^9]。
輸出格式（stdout）：
    第 1 行      q（靜態查詢個數）
    接下來 q 行  每行 `distinct sum mode`
    第 q+2 行    q2（帶修查詢個數）
    接下來 q2 行 每行 `distinct sum mode`
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import List, Optional, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 3000
VMAX = 10 ** 9


# ---------------------------------------------------------------- 整數方根


def isqrt_floor(v: int) -> int:
    """整數平方根下取整（二分，兩個語言一致）。"""
    lo, hi = 0, 1
    while hi * hi <= v:
        hi *= 2
    lo = hi // 2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid * mid <= v:
            lo = mid
        else:
            hi = mid - 1
    return lo


def icbrt_floor(v: int) -> int:
    """整數立方根下取整（二分）。"""
    lo, hi = 0, 1
    while hi * hi * hi <= v:
        hi *= 2
    lo = hi // 2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid * mid * mid <= v:
            lo = mid
        else:
            hi = mid - 1
    return lo


def block_size_plain(n: int) -> int:
    """普通莫隊塊長：√n（至少 1）。"""
    return max(1, isqrt_floor(n))


def block_size_modified(n: int) -> int:
    """帶修莫隊塊長：n^(2/3)（至少 1）。"""
    if n <= 1:
        return 1
    return max(1, icbrt_floor(n * n))


# ---------------------------------------------------------------- 區間統計量


class Counter:
    """維護 distinct / sum / mode（眾數出現次數）。"""

    __slots__ = ("cnt", "freq", "distinct", "total", "maxc")

    def __init__(self, m: int, n: int) -> None:
        # m = 離散化後的值的個數；n = 陣列長度（出現次數最多到 n）
        self.cnt = [0] * (m + 2)      # 以離散化後的編號為索引
        self.freq = [0] * (n + 4)     # freq[c] = 出現次數恰為 c 的值有幾個
        self.distinct = 0
        self.total = 0
        self.maxc = 0

    def add(self, cid: int, val: int) -> None:
        c = self.cnt[cid]
        self.freq[c] -= 1
        c += 1
        self.cnt[cid] = c
        self.freq[c] += 1
        if c == 1:
            self.distinct += 1
        if c > self.maxc:
            self.maxc = c
        self.total += val

    def remove(self, cid: int, val: int) -> None:
        c = self.cnt[cid]
        self.freq[c] -= 1
        c -= 1
        self.cnt[cid] = c
        self.freq[c] += 1
        if c == 0:
            self.distinct -= 1
        while self.maxc > 0 and self.freq[self.maxc] == 0:
            self.maxc -= 1
        self.total -= val

    def triple(self) -> Tuple[int, int, int]:
        return (self.distinct, self.total, self.maxc)


# ---------------------------------------------------------------- ① 普通莫隊


def mo_queries(a: Sequence[int], queries: Sequence[Tuple[int, int]],
               cid: Sequence[int], m: int) -> List[Tuple[int, int, int]]:
    """queries 為 0-indexed [l, r] 閉區間；回傳與輸入同順序的 (distinct, sum, mode)。"""
    n = len(a)
    q = len(queries)
    ans: List[Optional[Tuple[int, int, int]]] = [None] * q
    if n == 0 or q == 0:
        return [(0, 0, 0)] * q

    B = block_size_plain(n)
    order = list(range(q))
    order.sort(key=lambda i: (queries[i][0] // B,
                              -queries[i][1] if (queries[i][0] // B) % 2 == 1
                              else queries[i][1],
                              i))

    ctr = Counter(m, n)
    L, R = 0, -1
    for idx in order:
        l, r = queries[idx]
        while L > l:
            L -= 1
            ctr.add(cid[L], a[L])
        while R < r:
            R += 1
            ctr.add(cid[R], a[R])
        while L < l:
            ctr.remove(cid[L], a[L])
            L += 1
        while R > r:
            ctr.remove(cid[R], a[R])
            R -= 1
        ans[idx] = ctr.triple()
    return [x if x is not None else (0, 0, 0) for x in ans]


# ---------------------------------------------------------------- ② 帶修莫隊


def mo_with_updates(a: Sequence[int], ops: Sequence[Tuple[int, int, int]],
                    idof: "dict", m: int) -> List[Tuple[int, int, int]]:
    """ops：`(1, l, r)` 查詢（1-indexed 閉區間）或 `(2, p, x)` 修改（1-indexed 位置）。

    回傳依序的查詢結果 (distinct, sum, mode)。
    """
    n = len(a)
    # 先正規化操作串：查詢記 (l0, r0, t)；修改記 (p0, nv, nid, ov, oid)
    cur_val = list(a)
    cur_id = [idof[v] for v in a]
    qs: List[Tuple[int, int, int]] = []          # (l0, r0, t)
    mods: List[Tuple[int, int, int, int, int]] = []
    t = 0
    for op, x, y in ops:
        if op == 1:
            if n == 0:
                qs.append((0, -1, t))            # 空區間
            else:
                lo = max(0, min(x - 1, n - 1))
                hi = max(0, min(y - 1, n - 1))
                if lo > hi:
                    lo, hi = hi, lo
                qs.append((lo, hi, t))
        elif op == 2:
            if 1 <= x <= n:
                p = x - 1
                mods.append((p, y, idof[y], cur_val[p], cur_id[p]))
                cur_val[p] = y
                cur_id[p] = idof[y]
                t += 1
    if n == 0:
        return [(0, 0, 0)] * len(qs)
    if not qs:
        return []

    # 復原到初始狀態（正規化時模擬過一次修改，這裡重來）
    cur_val = list(a)
    cur_id = [idof[v] for v in a]

    B = block_size_modified(n)
    order = list(range(len(qs)))
    order.sort(key=lambda i: (qs[i][0] // B, qs[i][1] // B, qs[i][2], i))

    ctr = Counter(m, n)
    L, R, T = 0, -1, 0
    res: List[Optional[Tuple[int, int, int]]] = [None] * len(qs)

    def change_time(target: int) -> None:
        nonlocal T
        while T < target:
            p, nv, nid, ov, oid = mods[T]
            if L <= p <= R:
                ctr.remove(oid, ov)
                ctr.add(nid, nv)
            cur_val[p] = nv
            cur_id[p] = nid
            T += 1
        while T > target:
            T -= 1
            p, nv, nid, ov, oid = mods[T]
            if L <= p <= R:
                ctr.remove(nid, nv)
                ctr.add(oid, ov)
            cur_val[p] = ov
            cur_id[p] = oid

    for idx in order:
        l, r, tt = qs[idx]
        change_time(tt)
        while L > l:
            L -= 1
            ctr.add(cur_id[L], cur_val[L])
        while R < r:
            R += 1
            ctr.add(cur_id[R], cur_val[R])
        while L < l:
            ctr.remove(cur_id[L], cur_val[L])
            L += 1
        while R > r:
            ctr.remove(cur_id[R], cur_val[R])
            R -= 1
        res[idx] = ctr.triple()

    return [x if x is not None else (0, 0, 0) for x in res]


# ---------------------------------------------------------------- 主流程


def run(a: Sequence[int], queries: Sequence[Tuple[int, int]],
        ops: Sequence[Tuple[int, int, int]]) -> str:
    # 離散化：初始陣列 + 所有修改的新值（兩個語言的集合完全相同）
    pool = list(a) + [y for op, x, y in ops if op == 2]
    uniq = sorted(set(pool))
    idof = {v: i for i, v in enumerate(uniq)}
    cid = [idof[v] for v in a]
    m = len(uniq)

    ra = mo_queries(a, queries, cid, m)
    rb = mo_with_updates(a, ops, idof, m)

    lines = [str(len(ra))]
    lines += [f"{d} {s} {m}" for d, s, m in ra]
    lines.append(str(len(rb)))
    lines += [f"{d} {s} {m}" for d, s, m in rb]
    return "\n".join(lines)


# ---------------------------------------------------------------- 對拍用暴力


def bf_triple(vals: Sequence[int]) -> Tuple[int, int, int]:
    if not vals:
        return (0, 0, 0)
    cnt: dict = {}
    for v in vals:
        cnt[v] = cnt.get(v, 0) + 1
    return (len(cnt), sum(vals), max(cnt.values()))


def bf_run(a: Sequence[int], queries: Sequence[Tuple[int, int]],
           ops: Sequence[Tuple[int, int, int]]) -> str:
    ra = [bf_triple(a[l:r + 1]) for l, r in queries]

    cur = list(a)
    rb: List[Tuple[int, int, int]] = []
    for op, x, y in ops:
        if op == 1:
            if len(cur) == 0:
                rb.append((0, 0, 0))
            else:
                l = max(0, min(x - 1, len(cur) - 1))
                r = max(0, min(y - 1, len(cur) - 1))
                if l > r:
                    l, r = r, l
                rb.append(bf_triple(cur[l:r + 1]))
        elif op == 2:
            if 1 <= x <= len(cur):
                cur[x - 1] = y

    lines = [str(len(ra))]
    lines += [f"{d} {s} {m}" for d, s, m in ra]
    lines.append(str(len(rb)))
    lines += [f"{d} {s} {m}" for d, s, m in rb]
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check(a: Sequence[int], queries: Sequence[Tuple[int, int]],
          ops: Sequence[Tuple[int, int, int]]) -> None:
    got = run(a, queries, ops)
    want = bf_run(a, queries, ops)
    assert got == want, (list(a), list(queries), list(ops), got, want)


def run_tests() -> None:
    # ---- 固定用例 ----
    assert run([], [], []) == "0\n0"
    assert run([1, 2, 3], [], []) == "0\n0"
    assert run([1, 1, 2], [(0, 2)], []) == "1\n2 4 2\n0"
    assert run([5], [(0, 0)], [(1, 1, 1)]) == "1\n1 5 1\n1\n1 5 1"
    assert run([1, 2], [(0, 1)], [(2, 1, 9), (1, 1, 2)]) == (
        "1\n2 3 1\n1\n2 11 1")
    assert run([1, 2, 3, 4], [(0, 3)],
               [(1, 1, 4), (2, 2, 9), (1, 1, 4), (2, 4, 9), (1, 1, 4)]) == (
        "1\n4 10 1\n3\n4 10 1\n4 17 1\n3 22 2")

    check([], [], [])
    check([1], [], [])
    check([1], [(0, 0)], [])
    check([1, 2, 3], [(0, 2), (1, 1), (0, 0)], [])
    check([1, 1, 1, 2, 2], [(0, 4), (1, 3), (2, 2), (0, 1)], [])
    check([3, 1, 4, 1, 5, 9, 2, 6],
          [(0, 7), (2, 5), (0, 0), (3, 3), (1, 6)], [])
    check([1, 2, 3], [], [(1, 1, 3)])
    check([1, 2, 3], [], [(2, 1, 5), (1, 1, 3)])
    check([1, 2, 3, 4], [(0, 3)],
          [(1, 1, 4), (2, 2, 9), (1, 1, 4), (2, 4, 9), (1, 1, 4)])

    rng = random.Random(20261014)

    # ---- 隨機對拍：只有靜態查詢 ----
    for _ in range(260):
        n = rng.randint(0, 12)
        a = [rng.randint(0, 4) for _ in range(n)]
        q = rng.randint(0, 10)
        queries = []
        for _ in range(q):
            l = rng.randint(-2, n + 2)
            r = rng.randint(-2, n + 2)
            if n > 0:
                l = max(0, min(l, n - 1))
                r = max(0, min(r, n - 1))
            else:
                l = r = 0
            if l > r:
                l, r = r, l
            queries.append((l, r))
        check(a, queries, [])

    # ---- 隨機對拍：帶修改 ----
    for _ in range(260):
        n = rng.randint(0, 12)
        a = [rng.randint(0, 4) for _ in range(n)]
        m = rng.randint(0, 12)
        ops = []
        for _ in range(m):
            if rng.random() < 0.5:
                ops.append((1, rng.randint(-1, n + 1), rng.randint(-1, n + 1)))
            else:
                ops.append((2, rng.randint(0, n + 1), rng.randint(0, 4)))
        check(a, [], ops)

    # ---- 隨機對拍：兩者混合 ----
    for _ in range(160):
        n = rng.randint(1, 14)
        a = [rng.randint(0, 6) for _ in range(n)]
        queries = []
        for _ in range(rng.randint(0, 8)):
            l = rng.randint(0, n - 1)
            r = rng.randint(0, n - 1)
            if l > r:
                l, r = r, l
            queries.append((l, r))
        ops = []
        for _ in range(rng.randint(0, 12)):
            if rng.random() < 0.5:
                ops.append((1, rng.randint(1, n), rng.randint(1, n)))
            else:
                ops.append((2, rng.randint(1, n), rng.randint(0, 6)))
        check(a, queries, ops)

    # ---- 區塊大小的邊界 ----
    for n in (1, 2, 3, 4, 8, 9, 16, 27, 64, 100):
        a = [rng.randint(0, 5) for _ in range(n)]
        queries = [(0, n - 1), (n // 2, n - 1), (0, 0)] if n else []
        ops = [(1, 1, n), (2, 1, 7), (1, 1, n)] if n else []
        check(a, queries, ops)

    # ---- 中量級效能 ----
    n = 1500
    a = [rng.randint(0, 300) for _ in range(n)]
    queries = []
    for _ in range(400):
        l = rng.randint(0, n - 1)
        r = rng.randint(0, n - 1)
        if l > r:
            l, r = r, l
        queries.append((l, r))
    ops = []
    for _ in range(400):
        if rng.random() < 0.5:
            ops.append((1, rng.randint(1, n), rng.randint(1, n)))
        else:
            ops.append((2, rng.randint(1, n), rng.randint(0, 300)))
    out = run(a, queries, ops)
    assert out.split("\n")[0] == "400"
    # 抽幾個靜態查詢與暴力對照
    for _ in range(5):
        i = rng.randrange(len(queries))
        l, r = queries[i]
        tri = bf_triple(a[l:r + 1])
        assert f"{tri[0]} {tri[1]} {tri[2]}" in out


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

    def clamp_n(v: Optional[int]) -> int:
        if v is None:
            return 0
        return max(0, min(MAXN, v))

    n = clamp_n(read_int())
    a: List[int] = []
    for _ in range(n):
        v = read_int()
        if v is None:
            break
        a.append(clamp_val(v))

    q = clamp_n(read_int())
    queries: List[Tuple[int, int]] = []
    for _ in range(q):
        l = read_int()
        r = read_int()
        if l is None or r is None:
            break
        if n == 0:
            queries.append((0, 0))
            continue
        lo = max(0, min(clamp_val(l) - 1, n - 1))
        hi = max(0, min(clamp_val(r) - 1, n - 1))
        if lo > hi:
            lo, hi = hi, lo
        queries.append((lo, hi))

    m = clamp_n(read_int())
    ops: List[Tuple[int, int, int]] = []
    for _ in range(m):
        o = read_int()
        x = read_int()
        y = read_int()
        if o is None or x is None or y is None:
            break
        oc = max(0, min(2, clamp_val(o)))
        ops.append((oc, clamp_val(x), clamp_val(y)))

    emit(run(a, queries, ops) + "\n")


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
