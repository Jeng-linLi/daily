"""整體二分（Overall Binary Search / 離線整體二分）：單點修改 + 區間第 k 小

題意：
    給定一個長度為 n 的數列 a，初始值為 a[1..n]。接著依序執行 m 個操作：
        1 p x   —— 把 a[p] 改成 x
        2 l r k —— 查詢「當下」a[l..r] 中第 k 小的值（1 <= k <= r-l+1）
    要求把所有查詢的答案依序輸出。
    限制：1 <= n, m <= 100000，|x| <= 10^9。

思路：
    ### 1. 為什麼不能用線段樹套值域樹
    「動態第 k 小」最直覺的做法是主席樹（可持久化線段樹），但那是**線上**解法：
    每個版本要 O(log V) 空間，n + m 個版本就是 O((n+m) log V) 個節點。
    本題是純離線的，所以可以用更省空間的 **整體二分**。

    ### 2. 離線二分的基本套路
    若我們能回答「查詢範圍 [l,r] 中，值 <= mid 的元素有幾個」，
    就能判斷第 k 小的值落在 [L,mid] 還是 [mid+1,R]（落在右邊時把 k 減掉左邊的個數）。
    問題是：每一層二分都要重新回答一次「有幾個」，
    若每次都重掃全部操作，複雜度會變成 O(m log V * m)，不可行。

    ### 3. 整體二分的核心：所有查詢共用同一棵二分樹
    把所有操作（含修改）依時間順序排成一個事件序列，遞迴處理值域 [vl, vr]：
        mid = (vl + vr) / 2
        從頭到尾掃描事件序列，同時用**樹狀數組**維護「目前值落在 [vl,mid] 的位置數」：
            * 修改事件：若新值的秩 <= mid，就 BIT[pos] += delta（並把事件留在左半邊）
                         否則留在右半邊（不動 BIT）
            * 查詢事件：cnt = BIT.sum(r) - BIT.sum(l-1)
                         若 k <= cnt  → 事件留在左半邊
                         否則          → k -= cnt，事件留在右半邊
        掃完後把這一層加進 BIT 的量**全部回滾**，再各自遞迴。
    因為事件序列在遞迴中保持時間順序，BIT 在任何時刻都精確反映
    「該時刻值落在 [vl,mid] 的位置數」，所以每一層只需掃一次序列，
    總複雜度降為 O((n+m) log V log n)。

    ### 4. 「修改」怎麼變成可以被 BIT 累加的事件
    BIT 只能做「加一點 / 減一點」，所以把賦值拆成「刪舊值 + 加新值」：
        a[p]: old → new   等價於   DELETE(p, old)  然後  INSERT(p, new)
    並且在最前面補上所有初始值的 INSERT(p, a[p])。
    這樣每個位置的「當前值是否落在 [vl,mid]」永遠是 0 或 1，
    INSERT/DELETE 把它維護得對稱，兩邊相消就是 0。
    注意 DELETE 事件必須帶「舊值的秩」，所以要先掃一遍操作序列求出每個位置的前值。

    ### 5. 為什麼答案不會有平手問題
    第 k 小指的是**值**（不是位置），值相同也只佔一格，
    所以答案在兩個語言裡都唯一，不需要任何正規化處理。

複雜度：
    值域離散化（排序）      O((n+m) log(n+m))
    整體二分（每層掃一次序列）
        共 O(log V) 層，每層 O((n+m) log n)   →  O((n+m) log V log n)
    空間                    O(n + m + V)

應用場景（Applications）：
    - **即時排行榜 / 賽事名次**：分數被更新後查詢「第 k 名」，
      是整體二分最典型的線上查詢改寫（每輪結算後重跑一次離線批次）
    - **資料庫的 ORDER BY ... LIMIT k**：批次匯出前 N 名
    - **版本控制 / 分散式版本號排序**：找出第 k 舊的版本
    - **監控系統**：取一段時間視窗內的第 k 大延遲（p99 latency）
    - **資本市場**：組合持股第 k 大的權重（風險集中度監控）
"""

import sys
from typing import List, Optional, Sequence, Tuple

NMAX = 100000
MMAX = 100000
VAL_MAX = 10 ** 9


class Fenwick:
    """樹狀數組（Fenwick Tree）：單點加、前綴和。"""

    __slots__ = ("n", "tree")

    def __init__(self, n: int) -> None:
        self.n = n
        self.tree = [0] * (n + 1)

    def add(self, i: int, delta: int) -> None:
        tree = self.tree
        n = self.n
        while i <= n:
            tree[i] += delta
            i += i & -i

    def sum(self, i: int) -> int:
        """前綴和 sum(1..i)。"""
        tree = self.tree
        s = 0
        while i > 0:
            s += tree[i]
            i -= i & -i
        return s


# 事件格式：
#   [0, pos, delta, rank, 0]  —— 修改事件（INSERT: delta=+1，DELETE: delta=-1）
#   [1, l, r, k, qid]          —— 查詢事件（k 會在遞迴中被就地修改）
E_UPDATE = 0
E_QUERY = 1


def overall_binary_search(n: int, initial: Sequence[int],
                          ops: Sequence[Tuple]) -> List[int]:
    """整體二分求解。

    參數：
        n       —— 數列長度
        initial —— 初始值（長度 n，索引 0 起算）
        ops     —— 操作序列，元素為 ("U", pos, val) 或 ("Q", l, r, k)，
                   pos / l / r 皆為 1 起算
    回傳：每個查詢的答案（實際數值），順序與查詢出現順序相同。
    """
    # ---- 1. 值域離散化：初始值 + 所有修改後的值 ----
    all_vals = list(initial)
    for op in ops:
        if op[0] == "U":
            all_vals.append(op[2])
    vals = sorted(set(all_vals))
    vcount = len(vals)
    rank_of = {v: i for i, v in enumerate(vals)}

    # ---- 2. 建立事件序列（保持時間順序）----
    events: List[list] = []
    cur_rank = [0] * (n + 1)
    for i in range(1, n + 1):
        r = rank_of[initial[i - 1]]
        cur_rank[i] = r
        events.append([E_UPDATE, i, 1, r, 0])          # INSERT(p, a[p])

    qid = 0
    for op in ops:
        if op[0] == "U":
            _, pos, val = op
            new_r = rank_of[val]
            old_r = cur_rank[pos]
            if old_r != new_r:
                events.append([E_UPDATE, pos, -1, old_r, 0])   # DELETE(p, old)
                events.append([E_UPDATE, pos, 1, new_r, 0])    # INSERT(p, new)
                cur_rank[pos] = new_r
        else:
            _, l, r, k = op
            events.append([E_QUERY, l, r, k, qid])
            qid += 1

    answers = [0] * qid
    if qid == 0:
        return answers

    bit = Fenwick(n)

    # ---- 3. 遞迴 ----
    def solve(vl: int, vr: int, evs: List[list]) -> None:
        if not evs:
            return
        if vl == vr:
            value = vals[vl]
            for e in evs:
                if e[0] == E_QUERY:
                    answers[e[4]] = value
            return

        mid = (vl + vr) >> 1
        left: List[list] = []
        right: List[list] = []
        for e in evs:
            if e[0] == E_UPDATE:
                if e[3] <= mid:
                    left.append(e)
                    bit.add(e[1], e[2])
                else:
                    right.append(e)
            else:
                cnt = bit.sum(e[2]) - bit.sum(e[1] - 1)
                if e[3] <= cnt:
                    left.append(e)
                else:
                    e[3] -= cnt
                    right.append(e)

        # 回滾這一層對 BIT 的所有影響
        for e in left:
            if e[0] == E_UPDATE:
                bit.add(e[1], -e[2])
        evs.clear()  # 釋放父層的事件列表，控制記憶體

        solve(vl, mid, left)
        solve(mid + 1, vr, right)

    solve(0, vcount - 1, events)
    return answers


def brute_force(n: int, initial: Sequence[int],
                ops: Sequence[Tuple]) -> List[int]:
    """參考實作：直接維護數列 + 排序取第 k 小（O(n log n) 每次查詢，僅供對拍）。"""
    a = list(initial)
    out: List[int] = []
    for op in ops:
        if op[0] == "U":
            _, pos, val = op
            a[pos - 1] = val
        else:
            _, l, r, k = op
            out.append(sorted(a[l - 1:r])[k - 1])
    return out


def run_tests() -> None:
    # ---- 固定用例 ----
    # 1) 沒有修改，查詢 [1,3] 第 2 小 → 2
    assert overall_binary_search(3, [3, 1, 2], [("Q", 1, 3, 2)]) == [2]
    # 2) 全部相同
    assert overall_binary_search(4, [7, 7, 7, 7], [("Q", 1, 4, 3)]) == [7]
    # 3) 修改後查詢：a = [1,2,3]，a[1] := 10，問 [1,3] 第 1 小 → 2
    assert overall_binary_search(3, [1, 2, 3], [("U", 1, 10), ("Q", 1, 3, 1)]) == [2]
    # 4) 修改成最小值
    assert overall_binary_search(3, [5, 6, 7], [("U", 3, -1), ("Q", 1, 3, 1)]) == [-1]
    # 5) 修改成最大值：a = [5,6,7] -> [100,6,7]，問 [1,3] 第 3 小 → 100
    assert overall_binary_search(3, [5, 6, 7], [("U", 1, 100), ("Q", 1, 3, 3)]) == [100]
    # 6) 負數與 0
    assert overall_binary_search(5, [-3, 0, -7, 4, 0], [("Q", 1, 5, 1)]) == [-7]
    assert overall_binary_search(5, [-3, 0, -7, 4, 0], [("Q", 1, 5, 5)]) == [4]
    # 7) 區間只有一個元素
    assert overall_binary_search(3, [9, 1, 4], [("Q", 2, 2, 1)]) == [1]
    # 8) 連續多次修改同一位置：a[1] 依序 8 -> 2 -> 5，問 [1,3] 第 2 小 → 3
    ops = [("U", 1, 8), ("U", 1, 2), ("U", 1, 5), ("Q", 1, 3, 2)]
    assert overall_binary_search(3, [1, 2, 3], ops) == [3]
    # 9) 修改成原本相同的值（no-op）
    assert overall_binary_search(3, [1, 2, 3], [("U", 2, 2), ("Q", 1, 3, 1)]) == [1]
    # 10) 查詢在所有修改之後
    assert overall_binary_search(4, [10, 20, 30, 40],
                                 [("U", 1, 5), ("U", 2, 6), ("Q", 1, 4, 2)]) == [6]
    # 11) 連續兩次查詢中間夾修改：第一次問 [1,3] 第 3 小 → 3，
    #     修改 a[1]:=100 後再問 → 100
    ops = [("Q", 1, 3, 3), ("U", 1, 100), ("Q", 1, 3, 3)]
    assert overall_binary_search(3, [1, 2, 3], ops) == [3, 100]
    # 12) 只有查詢沒有修改，n = 1
    assert overall_binary_search(1, [42], [("Q", 1, 1, 1)]) == [42]

    # ---- 與參考實作對拍 ----
    import random

    rng = random.Random(20261102)
    for _ in range(300):
        n = rng.randint(1, 9)
        initial = [rng.randint(-12, 12) for _ in range(n)]
        ops = []
        for _ in range(rng.randint(0, 12)):
            if rng.random() < 0.45:
                ops.append(("U", rng.randint(1, n), rng.randint(-12, 12)))
            else:
                l = rng.randint(1, n)
                r = rng.randint(l, n)
                k = rng.randint(1, r - l + 1)
                ops.append(("Q", l, r, k))
        got = overall_binary_search(n, initial, ops)
        want = brute_force(n, initial, ops)
        assert got == want, (initial, ops, got, want)

    # 值域較大且修改密集的情況
    for _ in range(120):
        n = rng.randint(1, 7)
        initial = [rng.randint(-1000, 1000) for _ in range(n)]
        ops = []
        for _ in range(rng.randint(1, 16)):
            if rng.random() < 0.6:
                ops.append(("U", rng.randint(1, n), rng.randint(-1000, 1000)))
            else:
                l = rng.randint(1, n)
                r = rng.randint(l, n)
                ops.append(("Q", l, r, rng.randint(1, r - l + 1)))
        got = overall_binary_search(n, initial, ops)
        want = brute_force(n, initial, ops)
        assert got == want, (initial, ops, got, want)

    # 大量查詢、零修改（退化成靜態區間第 k 小）
    for _ in range(60):
        n = rng.randint(1, 20)
        initial = [rng.randint(-50, 50) for _ in range(n)]
        ops = []
        for _ in range(rng.randint(1, 20)):
            l = rng.randint(1, n)
            r = rng.randint(l, n)
            ops.append(("Q", l, r, rng.randint(1, r - l + 1)))
        got = overall_binary_search(n, initial, ops)
        want = brute_force(n, initial, ops)
        assert got == want, (initial, ops, got, want)


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

    n_opt = _parse_int(nxt() or "")
    if n_opt is None or not (1 <= n_opt <= NMAX):
        bad()
        return
    m_opt = _parse_int(nxt() or "")
    if m_opt is None or not (0 <= m_opt <= MMAX):
        bad()
        return
    n = n_opt
    m = m_opt

    initial: List[int] = []
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None or abs(v) > VAL_MAX:
            bad()
            return
        initial.append(v)

    ops: List[Tuple] = []
    for _ in range(m):
        kind = nxt()
        if kind is None:
            bad()
            return
        if kind == "1":
            p = _parse_int(nxt() or "")
            x = _parse_int(nxt() or "")
            if p is None or x is None or not (1 <= p <= n) or abs(x) > VAL_MAX:
                bad()
                return
            ops.append(("U", p, x))
        elif kind == "2":
            l = _parse_int(nxt() or "")
            r = _parse_int(nxt() or "")
            k = _parse_int(nxt() or "")
            if l is None or r is None or k is None:
                bad()
                return
            if not (1 <= l <= r <= n) or not (1 <= k <= r - l + 1):
                bad()
                return
            ops.append(("Q", l, r, k))
        else:
            bad()
            return

    for ans in overall_binary_search(n, initial, ops):
        out.append(str(ans))
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
