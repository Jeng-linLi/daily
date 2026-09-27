"""最小生成树（Kruskal 算法 + 并查集）

题意：
    给定一个 n 个点、m 条边的**无向带权图**（顶点编号 0..n-1，边权为整数，可负），
    求一棵最小生成树：选出若干条边把全部点连通，且边权之和最小。
    若图本身不连通（不存在生成树），输出 disconnected。

思路：
    Kruskal 是**贪心 + 并查集**：把所有边按权值从小到大排序，依次尝试加入；
    如果这条边的两个端点当前还不连通，就把它选进生成树并合并两端点所在集合；
    否则这条边会和已选的边成环，直接丢弃。选够 n-1 条边就结束。

    为什么这样是对的（切分性质）：扫描到某条边 (u, v, w) 且 u、v 尚未连通时，
    把「u 所在的连通块」看作一个切分的一侧，那么所有跨这个切分的边里，
    (u, v, w) 是当前未处理的最小者 —— 因为更小的边全都已经被考虑过了，
    它们要么没能跨这个切分，要么会把点合并进来（那样 u、v 就已经连通了）。
    跨切分的最小边必然属于某棵最小生成树，所以选它不会错。

    并查集负责「u、v 是否已连通」的判定与合并，用路径压缩 + 按大小合并，
    均摊接近 O(1)。注意排序要**稳定**（同权边保持输入顺序），
    这样 Python 与 C++ 两版输出才能逐字节一致。

输入格式（stdin，所有数字按空白分隔）：
    n m
    u1 v1 w1
    u2 v2 w2
    ...（共 m 行，顶点 0-based）
输出格式（stdout）：
    连通时：第一行为最小生成树的总权重，随后每行一条被选中的边 `u v w`
            （按 Kruskal 的选中顺序输出，即权值升序、同权值按输入顺序）
    不连通时：只输出一行 disconnected
无 stdin 输入时运行内置断言测试并输出 `all tests passed`。
"""

import sys
from typing import Iterator, List, Optional, Tuple

Edge = Tuple[int, int, int]


class DSU:
    """并查集：路径压缩 + 按大小合并，同时维护连通块个数。"""

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n
        self.components = n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]   # 路径压缩（折半）
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> bool:
        """合并 a、b 所在集合；返回是否真的合并了（原本不连通才返回 True）。"""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.size[ra] < self.size[rb]:     # 按大小合并，小树挂到大树下
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        self.components -= 1
        return True


def kruskal(n: int, edges: List[Edge]) -> Tuple[int, List[Edge], bool]:
    """返回 (总权重, 被选中的边, 是否连通)。

    时间 O(m log m)（瓶颈在排序），空间 O(n + m)。
    """
    dsu = DSU(n)
    chosen: List[Edge] = []
    total = 0
    # 稳定排序：只按权值排序，同权值保持输入先后顺序
    for idx in sorted(range(len(edges)), key=lambda i: edges[i][2]):
        u, v, w = edges[idx]
        if dsu.union(u, v):
            chosen.append((u, v, w))
            total += w
            if len(chosen) == n - 1:
                break                          # 已经是一棵生成树，可以提前结束
    connected = (n <= 1) or (len(chosen) == n - 1)
    return total, chosen, connected


def prim(n: int, edges: List[Edge]) -> Tuple[int, bool]:
    """对照用的 Prim（邻接表 + O(n^2) 选最小），用于与 Kruskal 交叉验证总权重。"""
    if n == 0:
        return 0, True
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        if u == v:
            continue                           # 自环对 MST 无意义
        adj[u].append((v, w))
        adj[v].append((u, w))

    INF = float("inf")
    dist = [INF] * n
    used = [False] * n
    dist[0] = 0
    total = 0
    picked = 0
    for _ in range(n):
        best = -1
        for i in range(n):
            if not used[i] and (best == -1 or dist[i] < dist[best]):
                best = i
        if best == -1 or dist[best] == INF:
            return total, False                # 剩下的点都够不着，图不连通
        used[best] = True
        total += dist[best]
        picked += 1
        for v, w in adj[best]:
            if not used[v] and w < dist[v]:
                dist[v] = w
    return total, picked == n


def mst_brute(n: int, edges: List[Edge]) -> Tuple[Optional[int], bool]:
    """对照用的指数级枚举：枚举所有边子集，挑出权和最小的生成树。仅用于极小规模测试。"""
    if n <= 1:
        return 0, True
    m = len(edges)
    best: Optional[int] = None
    for mask in range(1 << m):
        if bin(mask).count("1") != n - 1:      # 生成树恰好 n-1 条边
            continue
        dsu = DSU(n)
        ok = True
        total = 0
        for i in range(m):
            if (mask >> i) & 1:
                u, v, w = edges[i]
                if dsu.union(u, v):
                    total += w
                else:
                    ok = False                 # 成环，不是树
                    break
        if ok and dsu.components == 1:
            if best is None or total < best:
                best = total
    return best, best is not None


def total_weight(edges: List[Edge]) -> int:
    return sum(w for _, _, w in edges)


# ---------------- IO ----------------

def _next_int(it: Iterator[str], default: int = 0) -> int:
    """取下一个整数；输入被截断时用默认值兜底，避免直接抛异常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    it = iter(data.split())
    n = _next_int(it)
    m = _next_int(it)
    edges: List[Edge] = []
    for _ in range(m):
        u = _next_int(it)
        v = _next_int(it)
        w = _next_int(it)
        edges.append((u, v, w))

    total, chosen, connected = kruskal(n, edges)
    if not connected:
        print("disconnected")
        return
    print(total)
    for u, v, w in chosen:
        print(f"{u} {v} {w}")


def run_tests() -> None:
    # README 示例：4 点 5 边，MST = (0,1,1) + (1,2,2) + (2,3,3) = 6
    edges = [(0, 1, 1), (0, 2, 4), (1, 2, 2), (1, 3, 5), (2, 3, 3)]
    total, chosen, ok = kruskal(4, edges)
    assert ok is True
    assert total == 6
    assert chosen == [(0, 1, 1), (1, 2, 2), (2, 3, 3)]
    assert prim(4, edges) == (6, True)
    assert mst_brute(4, edges) == (6, True)

    # 不连通：两个点之间只有一条边，第三个点孤立
    total2, chosen2, ok2 = kruskal(3, [(0, 1, 5)])
    assert ok2 is False
    assert prim(3, [(0, 1, 5)])[1] is False
    assert mst_brute(3, [(0, 1, 5)])[1] is False

    # 退化情形
    assert kruskal(0, []) == (0, [], True)            # 空图视为连通，权重 0
    assert kruskal(1, []) == (0, [], True)            # 单点，不需要边
    assert kruskal(2, [(0, 1, 7)]) == (7, [(0, 1, 7)], True)
    assert kruskal(2, []) == (0, [], False)           # 两点无边，不连通

    # 自环与重边：自环必被丢弃，重边只留一条
    assert kruskal(2, [(0, 0, 1), (0, 1, 3), (0, 1, 3)]) == (3, [(0, 1, 3)], True)
    assert kruskal(3, [(0, 1, 2), (1, 2, 2), (0, 2, 2)])[0] == 4

    # 负权边同样成立
    assert kruskal(3, [(0, 1, -5), (1, 2, -1), (0, 2, 10)]) == (-6, [(0, 1, -5), (1, 2, -1)], True)

    # 同权边按输入顺序稳定选中（保证两语言输出一致）
    same = [(2, 3, 1), (0, 1, 1), (1, 2, 1)]
    assert kruskal(4, same)[1] == [(2, 3, 1), (0, 1, 1), (1, 2, 1)]

    import random

    random.seed(20260927)

    # 随机对拍一：小规模图上 Kruskal、Prim、指数级枚举三者结果一致
    for _ in range(150):
        n = random.randint(1, 6)
        m = random.randint(0, 8)
        cand = [(u, v) for u in range(n) for v in range(u + 1, n)]
        random.shuffle(cand)
        edges = [(u, v, random.randint(1, 20)) for u, v in cand[:m]]
        total, chosen, ok = kruskal(n, edges)
        pt, pok = prim(n, edges)
        bt, bok = mst_brute(n, edges)
        assert ok == pok == bok
        if ok:
            assert total == pt == bt, (n, edges, total, pt, bt)
            assert len(chosen) == n - 1
            assert total_weight(chosen) == total
            # 选出来的边确实构成生成树：n-1 条且把所有点连成一块
            dsu = DSU(n)
            for u, v, _w in chosen:
                assert dsu.union(u, v) is True
            assert dsu.components == 1
            # 权值和不超过随便一棵生成树：与 Prim 再比一次即可
        else:
            assert total == 0 or len(chosen) < n - 1

    # 随机对拍二：保证连通的随机图（先造一条链再补随机边），Kruskal 与 Prim 对拍
    for _ in range(150):
        n = random.randint(2, 9)
        edges: List[Edge] = []
        for i in range(1, n):                  # 先连成链，确保一定连通
            j = random.randint(0, i - 1)
            edges.append((j, i, random.randint(1, 30)))
        for _ in range(random.randint(0, 6)):  # 再补一些随机边
            u, v = random.sample(range(n), 2)
            edges.append((u, v, random.randint(1, 30)))
        total, chosen, ok = kruskal(n, edges)
        pt, pok = prim(n, edges)
        assert ok is True and pok is True
        assert total == pt
        assert len(chosen) == n - 1
        dsu = DSU(n)
        for u, v, _w in chosen:
            assert dsu.union(u, v) is True
        assert dsu.components == 1

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
