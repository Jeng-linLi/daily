"""線段樹分治 + 可撤銷並查集：離線動態連通性

題意：
    n 個點，T 個時刻（1..T）。m 條邊，第 i 條邊 (u, v) 在時刻區間
    [l_i, r_i] 內存在。q 次詢問，每次給 (t, a, b)：問時刻 t 時 a 與 b
    是否連通。輸出 yes / no。

    輸入：第一行 n T m q；接下來 m 行每行 u v l r（1 <= l <= r <= T）；
          接下來 q 行每行 t a b（1-indexed）。
    輸出：每個詢問一行 yes 或 no。

思路：
    ### 1. 為什麼不能一邊加一邊查
    「支持撤銷」的並查集做不到路徑壓縮（撤銷會破壞壓縮結果），
    只能按秩/按大小合併，單次操作 O(log n)。若每個時刻都重算一次
    全圖連通性，總代價 O(T * m * log n)，動態加減邊更是要重來。

    ### 2. 線段樹分治（Segment Tree Divide & Conquer / Offline DSU）
    每條邊只在 [l, r] 內有效——把「時間軸」建成一顆線段樹，
    邊的存活區間 [l, r] 拆成 O(log T) 個線段樹節點，掛到節點上。
    然後 DFS 整顆線段樹：
    - 進入節點時：把掛在該節點的所有邊並入並查集（記錄棧深）；
    - 葉子節點：此時並查集恰好 = 該時刻的邊集，回答掛在此時刻的詢問；
    - 離開節點時：按棧回滾到進入前的狀態（只撤銷 merge，不真刪）。

    ### 3. 可撤銷並查集
    按大小合併（size 大的當根），每次真正改動的是 fa[小根] 和 size，
    用棧記下被改的 (節點, 舊父親)，撤銷時還原即可。**不路徑壓縮**，
    find 走純路徑遞歸/迭代，單次 O(log n)。

    ### 4. 對拍
    內建測試對每個 (邊集, 詢問) 用「逐時刻重建活躍邊集 + BFS」的
    暴力比對；另含固定用例（邊出現又消失、橋接邊後到等）。

複雜度：
    每條邊掛到 O(log T) 個節點，DFS 總並查集操作 O((m log T + q) log n)，
    時間 O(m log T log n + q log n)，空間 O(m log T + q)。

應用場景（Applications）：
    - **動態連通性**：社交網路中好友關係時斷時續，查詢某時刻兩人是否
      在同一圈子；網路鏈路按計畫時段啟停的故障影響分析
    - 同框架可換資料結構解「區間刪邊最短路」「線段樹分治 + 線性基」等
    - 若需要強制在線，得用 LCT（Link-Cut Tree），代價高得多
"""

import sys
from typing import Dict, List, Optional, Sequence, Tuple

NMAX = 300
TMAX = 300
VAL_LIM = 10 ** 18


class RollbackDSU:
    """可撤銷並查集：按大小合併，不路徑壓縮。"""

    def __init__(self, n: int) -> None:
        self.fa = list(range(n))
        self.sz = [1] * n
        self.stack: List[Tuple[int, int]] = []  # (被改的子節點, 其原父親)

    def find(self, x: int) -> int:
        while self.fa[x] != x:
            x = self.fa[x]
        return x

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.sz[ra] < self.sz[rb]:
            ra, rb = rb, ra
        self.stack.append((rb, self.fa[rb]))
        self.fa[rb] = ra
        self.sz[ra] += self.sz[rb]
        return True

    def snapshot(self) -> int:
        return len(self.stack)

    def rollback(self, mark: int) -> None:
        while len(self.stack) > mark:
            node, old_fa = self.stack.pop()
            self.sz[self.fa[node]] -= self.sz[node]
            self.fa[node] = old_fa


def dynamic_connectivity(
    n: int,
    T: int,
    edges: Sequence[Tuple[int, int, int, int]],
    queries: Sequence[Tuple[int, int, int]],
) -> List[str]:
    """edges: (u, v, l, r) 0-indexed 端點、時刻閉區間 [l, r]（0..T-1）。
    queries: (t, a, b)。回傳每個詢問 "yes"/"no"。"""
    # 線段樹：節點 1 為根，覆蓋 [0, T-1]
    size = 1
    while size < T:
        size *= 2
    tree: List[List[Tuple[int, int]]] = [[] for _ in range(2 * size)]

    def add_interval(node: int, nl: int, nr: int, l: int, r: int, e: Tuple[int, int]) -> None:
        if r < nl or nr < l:
            return
        if l <= nl and nr <= r:
            tree[node].append(e)
            return
        mid = (nl + nr) // 2
        add_interval(2 * node, nl, mid, l, r, e)
        add_interval(2 * node + 1, mid + 1, nr, l, r, e)

    for u, v, l, r in edges:
        if l > r:
            continue
        add_interval(1, 0, size - 1, l, r, (u, v))

    # 詢問掛到葉子（帶原下標，DFS 訪問順序 != 詢問順序）
    q_at: Dict[int, List[Tuple[int, int, int]]] = {}
    for idx, (t, a, b) in enumerate(queries):
        q_at.setdefault(t, []).append((idx, a, b))

    dsu = RollbackDSU(n)
    out: Dict[int, str] = {}

    def dfs(node: int, nl: int, nr: int) -> None:
        mark = dsu.snapshot()
        for u, v in tree[node]:
            dsu.union(u, v)
        if nl == nr:
            for idx, a, b in q_at.get(nl, []):
                out[idx] = "yes" if dsu.find(a) == dsu.find(b) else "no"
        else:
            mid = (nl + nr) // 2
            dfs(2 * node, nl, mid)
            dfs(2 * node + 1, mid + 1, nr)
        dsu.rollback(mark)

    if T > 0:
        dfs(1, 0, size - 1)
    return [out[i] for i in range(len(queries))]


# --------------------------------------------------------------- 暴力對拍

def brute_dynamic(
    n: int,
    T: int,
    edges: Sequence[Tuple[int, int, int, int]],
    queries: Sequence[Tuple[int, int, int]],
) -> List[str]:
    """逐時刻重建活躍邊集 + BFS 求連通。"""
    res: List[str] = []
    for t, a, b in queries:
        adj: List[List[int]] = [[] for _ in range(n)]
        for u, v, l, r in edges:
            if l <= t <= r:
                adj[u].append(v)
                adj[v].append(u)
        seen = [False] * n
        stack = [a]
        seen[a] = True
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if not seen[y]:
                    seen[y] = True
                    stack.append(y)
        res.append("yes" if seen[b] else "no")
    return res


# --------------------------------------------------------------- 內建測試

def run_tests() -> None:
    # 固定用例 1：邊 (0,1) 全程在，邊 (1,2) 只在時刻 2-3
    # 時刻 0: {0,1} 連通，2 不連；時刻 2: 三點全連通
    edges1 = [(0, 1, 0, 3), (1, 2, 2, 3)]
    qs1 = [(0, 0, 2), (2, 0, 2), (3, 1, 2), (1, 1, 2)]
    got = dynamic_connectivity(3, 4, edges1, qs1)
    assert got == ["no", "yes", "yes", "no"], f"case1 got {got}"

    # 固定用例 2：橋接邊後到 —— a-b 直連邊 [0,1]，繞路鏈 [2,4]
    edges2 = [(0, 1, 0, 1), (0, 2, 2, 4), (2, 1, 2, 4)]
    qs2 = [(0, 0, 1), (1, 0, 1), (2, 0, 1), (4, 0, 1), (3, 1, 2)]
    got = dynamic_connectivity(3, 5, edges2, qs2)
    assert got == ["yes", "yes", "yes", "yes", "yes"], f"case2 got {got}"

    # 固定用例 3：零詢問、零邊
    got = dynamic_connectivity(2, 1, [], [])
    assert got == [], "case3"

    # 隨機對拍：n <= 8, T <= 8
    rng_state = 24681357

    def rnd(mod: int) -> int:
        nonlocal rng_state
        rng_state = (rng_state * 1103515245 + 12345) % (1 << 31)
        return rng_state % mod

    for it in range(200):
        n = rnd(8) + 1
        T = rnd(8) + 1
        m = rnd(13)
        edges = []
        for _ in range(m):
            u = rnd(n)
            v = rnd(n)
            if u == v:
                continue
            l = rnd(T)
            r = l + rnd(T - l)
            edges.append((u, v, l, r))
        q = rnd(11)
        queries = [(rnd(T), rnd(n), rnd(n)) for _ in range(q)]
        got = dynamic_connectivity(n, T, edges, queries)
        want = brute_dynamic(n, T, edges, queries)
        assert got == want, f"random#{it} n={n} T={T} edges={edges} queries={queries} got={got} want={want}"


# --------------------------------------------------------------- IO 模式

def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


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
    if val > VAL_LIM:
        return None
    return val


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

    def get_int(lo: int, hi: int) -> Optional[int]:
        v = _parse_int(nxt() or "")
        if v is None or not (lo <= v <= hi):
            return None
        return v

    n = get_int(1, NMAX)
    if n is None:
        bad()
        return
    T = get_int(1, TMAX)
    if T is None:
        bad()
        return
    m = get_int(0, NMAX * NMAX)
    if m is None:
        bad()
        return
    q = get_int(0, 100000)
    if q is None:
        bad()
        return
    edges: List[Tuple[int, int, int, int]] = []
    for _ in range(m):
        u = get_int(1, n)
        v = get_int(1, n)
        l = get_int(1, T)
        r = get_int(1, T)
        if u is None or v is None or l is None or r is None or u == v:
            bad()
            return
        edges.append((u - 1, v - 1, l - 1, r - 1))
    queries: List[Tuple[int, int, int]] = []
    for _ in range(q):
        t = get_int(1, T)
        a = get_int(1, n)
        b = get_int(1, n)
        if t is None or a is None or b is None:
            bad()
            return
        queries.append((t - 1, a - 1, b - 1))
    for line in dynamic_connectivity(n, T, edges, queries):
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
