"""重鏈剖分 Heavy-Light Decomposition（HLD）：樹上路徑修改與查詢

題意：
    給定 n 個點的樹（根為 1 號點），每個點有點權。支援 q 筆操作：
        1 u v x —— 把 u 到 v 路徑上所有點的點權加上 x
        2 u v   —— 查詢 u 到 v 路徑上點權的最大值
        3 u v   —— 查詢 u 到 v 路徑上點權之和
    輸入：第一行 n q，第二行 n 個整數（點權，1 號點到 n 號點），
          接下來 n-1 行每行兩個整數表示一條樹邊，最後 q 行每行一筆操作。
    輸出：對每筆 2 / 3 操作各輸出一行。

思路：
    ### 1. 重輕劃分
    對每個點，其子節點中子樹最大的那個叫「重兒子」，其餘是「輕兒子」。
    重兒子之間的邊是「重邊」，由重邊連成的極長鏈叫「重鏈」。
    關鍵性質：從任意點走到根，經過的**輕邊**不超過 O(log n) 條——
    每走一條輕邊，所在子樹大小至少翻倍。

    ### 2. DFS 序上開線段樹
    對整棵樹 DFS（**重兒子優先**），給每個點標一個連續位置 pos[u]。
    每條重鏈內部位置連續，於是「鏈上區間」=「線段樹上連續區間」，
    樹上路徑問題就落到序列問題：區間加 / 區間最大值 / 區間和。

    ### 3. 路徑拆解
    求 u 到 v 的路徑：只要兩點不在同一條重鏈上，就把「鏈頭較深」的
    那個點跳到其鏈頭的父親（這一段鏈對應線段樹上一段連續區間），
    重複直到同鏈；最後處理同鏈上兩點之間的區間。
    每次跳一條輕邊，總跳數 O(log n)，每段在線段樹上 O(log n)，
    單次操作總複雜度 O(log^2 n)。

    ### 4. 對拍
    內建測試用暴力：路徑用「較深點逐步上爬到 LCA」逐點收集，
    直接對點權數組做加 / max / sum，與 HLD + 線段樹結果比對。

複雜度：
    預處理   O(n)
    單次操作 O(log^2 n)
    空間     O(n)

應用場景（Applications）：
    - **組織架構 / 匯報鏈**：任意兩人之間的晉升鏈上批量調整、統計
    - **網路路由**：兩節點間路徑上的鏈路批量加權 / 找瓶頸（最大值）
    - **版本依賴樹**：兩版本間的依賴鏈統計
    - HLD 也是「樹上問題轉序列問題」的通用框架，可套任何支援
      區間操作的線段樹（區間賦值、矩陣乘、Segment Tree Beats 等）
"""

import sys
from typing import List, Optional, Sequence, Tuple

NMAX = 20000
QMAX = 20000
NEG_INF = -(1 << 62)
ABS_LIM = 10 ** 18


class SegTree:
    """線段樹：區間加、區間和、區間最大值（lazy 下傳）。"""

    def __init__(self, n: int, vals: Sequence[int]) -> None:
        self.n = n
        self.sm = [0] * (4 * n)
        self.mx = [NEG_INF] * (4 * n)
        self.lz = [0] * (4 * n)
        self._vals = vals
        if n > 0:
            self._build(1, 0, n - 1)

    def _build(self, o: int, l: int, r: int) -> None:
        if l == r:
            self.sm[o] = self._vals[l]
            self.mx[o] = self._vals[l]
            return
        m = (l + r) >> 1
        self._build(o * 2, l, m)
        self._build(o * 2 + 1, m + 1, r)
        self._pull(o)

    def _pull(self, o: int) -> None:
        self.sm[o] = self.sm[o * 2] + self.sm[o * 2 + 1]
        self.mx[o] = max(self.mx[o * 2], self.mx[o * 2 + 1])

    def _apply(self, o: int, l: int, r: int, x: int) -> None:
        self.sm[o] += x * (r - l + 1)
        self.mx[o] += x
        self.lz[o] += x

    def _push(self, o: int, l: int, r: int) -> None:
        if self.lz[o]:
            m = (l + r) >> 1
            self._apply(o * 2, l, m, self.lz[o])
            self._apply(o * 2 + 1, m + 1, r, self.lz[o])
            self.lz[o] = 0

    def add(self, o: int, l: int, r: int, ql: int, qr: int, x: int) -> None:
        if qr < l or r < ql:
            return
        if ql <= l and r <= qr:
            self._apply(o, l, r, x)
            return
        self._push(o, l, r)
        m = (l + r) >> 1
        self.add(o * 2, l, m, ql, qr, x)
        self.add(o * 2 + 1, m + 1, r, ql, qr, x)
        self._pull(o)

    def query_sum(self, o: int, l: int, r: int, ql: int, qr: int) -> int:
        if qr < l or r < ql:
            return 0
        if ql <= l and r <= qr:
            return self.sm[o]
        self._push(o, l, r)
        m = (l + r) >> 1
        return (self.query_sum(o * 2, l, m, ql, qr) +
                self.query_sum(o * 2 + 1, m + 1, r, ql, qr))

    def query_max(self, o: int, l: int, r: int, ql: int, qr: int) -> int:
        if qr < l or r < ql:
            return NEG_INF
        if ql <= l and r <= qr:
            return self.mx[o]
        self._push(o, l, r)
        m = (l + r) >> 1
        return max(self.query_max(o * 2, l, m, ql, qr),
                   self.query_max(o * 2 + 1, m + 1, r, ql, qr))


class HLD:
    """重鏈剖分 + 線段樹（0 起算，根為 0 號點）。"""

    def __init__(self, n: int, parents: Sequence[int],
                 vals: Sequence[int]) -> None:
        self.n = n
        self.par = list(parents)
        self.dep = [0] * n
        self.sz = [1] * n
        self.heavy = [-1] * n
        self.head = [0] * n
        self.pos = [0] * n
        # 子節點表
        ch: List[List[int]] = [[] for _ in range(n)]
        for v in range(1, n):
            ch[self.par[v]].append(v)
        # BFS 定深
        order = [0]
        for u in order:            # order 邊遍歷邊擴展（BFS）
            for w in ch[u]:
                self.dep[w] = self.dep[u] + 1
                order.append(w)
        # 逆序求子樹大小與重兒子
        for u in reversed(order):
            for w in ch[u]:
                self.sz[u] += self.sz[w]
            best, bw = -1, -1
            for w in ch[u]:
                if self.sz[w] > best:
                    best, bw = self.sz[w], w
            self.heavy[u] = bw
        # 重兒子優先的迭代 DFS，分配 pos 與 head
        pos = 0
        stack = [0]
        while stack:
            u = stack.pop()
            self.pos[u] = pos
            pos += 1
            self.head[u] = u if self.par[u] == -1 or \
                self.heavy[self.par[u]] != u else self.head[self.par[u]]
            # 先壓輕兒子（後處理），再壓重兒子（先處理）
            for w in ch[u]:
                if w != self.heavy[u]:
                    stack.append(w)
            if self.heavy[u] != -1:
                stack.append(self.heavy[u])
        # 按 pos 重排點權，建線段樹
        ordered = [0] * n
        for u in range(n):
            ordered[self.pos[u]] = vals[u]
        self.seg = SegTree(n, ordered)

    # ---- 三種路徑操作（u、v 為 0 起算）----
    def _jump(self, u: int, v: int) -> List[Tuple[int, int]]:
        """把 u-v 路徑拆成若干 [pos_l, pos_r] 區間（pos_l <= pos_r）。"""
        segs: List[Tuple[int, int]] = []
        while self.head[u] != self.head[v]:
            if self.dep[self.head[u]] < self.dep[self.head[v]]:
                u, v = v, u
            segs.append((self.pos[self.head[u]], self.pos[u]))
            u = self.par[self.head[u]]
        if self.dep[u] > self.dep[v]:
            u, v = v, u
        segs.append((self.pos[u], self.pos[v]))
        return segs

    def path_add(self, u: int, v: int, x: int) -> None:
        for l, r in self._jump(u, v):
            self.seg.add(1, 0, self.n - 1, l, r, x)

    def path_sum(self, u: int, v: int) -> int:
        return sum(self.seg.query_sum(1, 0, self.n - 1, l, r)
                   for l, r in self._jump(u, v))

    def path_max(self, u: int, v: int) -> int:
        return max(self.seg.query_max(1, 0, self.n - 1, l, r)
                   for l, r in self._jump(u, v))


def _path_nodes(par: List[int], dep: List[int], u: int, v: int) -> List[int]:
    """暴力：收集 u-v 路徑上所有點。"""
    out: List[int] = []
    while dep[u] > dep[v]:
        out.append(u)
        u = par[u]
    while dep[v] > dep[u]:
        out.append(v)
        v = par[v]
    while u != v:
        out.append(u)
        out.append(v)
        u, v = par[u], par[v]
    out.append(u)
    return out


def run_tests() -> None:
    # ---- 固定用例：一條鏈 0-1-2-3-4 ----
    h = HLD(5, [-1, 0, 1, 2, 3], [1, 2, 3, 4, 5])
    assert h.path_sum(0, 4) == 15
    assert h.path_max(0, 4) == 5
    h.path_add(1, 3, 10)
    assert h.path_sum(0, 4) == 45
    assert h.path_max(2, 4) == 14
    assert h.path_max(0, 1) == 12

    # ---- 星形：0 連 1..4 ----
    s = HLD(5, [-1, 0, 0, 0, 0], [5, 1, 4, 2, 3])
    assert s.path_sum(1, 2) == 10   # 1 + 5 + 4（路徑過根）
    assert s.path_max(1, 3) == 5
    assert s.path_sum(0, 3) == 7
    s.path_add(0, 4, 3)
    assert s.path_sum(0, 4) == 14   # 8 + 6
    assert s.path_max(0, 4) == 8

    # ---- 與暴力對拍（自製 LCG，兩語言同一序列）----
    state = 20261104

    def nxt(mod: int) -> int:
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mod

    for _ in range(120):
        n = 2 + nxt(60)                       # 2..61
        par = [0] * n
        vals = [nxt(201) - 100 for _ in range(n)]  # [-100, 100]
        par[0] = -1
        for v in range(1, n):
            par[v] = nxt(v)                   # 隨機父親（< v）
        h = HLD(n, par, vals)
        dep = h.dep
        for _ in range(30):
            u, v = nxt(n), nxt(n)
            op = nxt(3)
            if op == 0:
                x = nxt(21) - 10
                h.path_add(u, v, x)
                for w in _path_nodes(par, dep, u, v):
                    vals[w] += x
            else:
                got = h.path_max(u, v) if op == 1 else h.path_sum(u, v)
                nodes = _path_nodes(par, dep, u, v)
                want = max(vals[w] for w in nodes) if op == 1 \
                    else sum(vals[w] for w in nodes)
                assert got == want, (n, par, u, v, got, want)


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
    if val > ABS_LIM:
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

    def pint() -> Optional[int]:
        return _parse_int(nxt() or "")

    n = pint()
    if n is None or not (1 <= n <= NMAX):
        bad()
        return
    q = pint()
    if q is None or not (0 <= q <= QMAX):
        bad()
        return

    vals: List[int] = []
    for _ in range(n):
        x = pint()
        if x is None or abs(x) > 10 ** 9:
            bad()
            return
        vals.append(x)

    edges: List[Tuple[int, int]] = []
    for _ in range(n - 1):
        a = pint()
        b = pint()
        if a is None or b is None or not (1 <= a <= n and 1 <= b <= n):
            bad()
            return
        edges.append((a - 1, b - 1))

    # 由無向邊構建以 0 為根的有根樹（BFS）
    adj: List[List[int]] = [[] for _ in range(n)]
    for a, b in edges:
        if a == b:
            bad()
            return
        adj[a].append(b)
        adj[b].append(a)
    parents = [-1] * n
    seen = [False] * n
    seen[0] = True
    bfs = [0]
    for u in bfs:
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parents[w] = u
                bfs.append(w)
    if not all(seen):              # 不連通
        bad()
        return

    h = HLD(n, parents, vals)
    out: List[str] = []
    for _ in range(q):
        op = pint()
        if op == 1:
            u = pint()
            v = pint()
            x = pint()
            if None in (u, v, x) or abs(x) > 10 ** 9:
                bad()
                return
            if not (1 <= u <= n and 1 <= v <= n):
                bad()
                return
            h.path_add(u - 1, v - 1, x)
        elif op in (2, 3):
            u = pint()
            v = pint()
            if u is None or v is None or not (1 <= u <= n and 1 <= v <= n):
                bad()
                return
            r = h.path_max(u - 1, v - 1) if op == 2 \
                else h.path_sum(u - 1, v - 1)
            out.append(str(r))
        else:
            bad()
            return
    for line in out:
        emit(line + "\n")


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
