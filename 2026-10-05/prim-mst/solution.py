"""Prim 最小生成樹（樸素 O(V^2) + 堆優化 O(E log V)，並與 Kruskal 交叉驗證）

題意：
    給定一張 **無向帶權連通（或不連通）圖**，求最小生成樹：
      1. MST 的總權重；
      2. 圖是否連通（不連通則不存在生成樹）；
      3. MST 具體選了哪些邊（輸出時按 `(u, v)` 字典序，`u < v`）；
      4. 若圖不連通，額外報告「連通分量個數」與「每個分量的最小生成樹權重」（最小生成森林）。
    同時提供三種互相獨立的實現彼此對拍：
      - `prim_naive`：鄰接矩陣 + 每輪線性掃描，O(V^2)，稠密圖最優；
      - `prim_heap`：鄰接表 + 二元堆（懶刪除），O(E log V)，稀疏圖最優；
      - `kruskal`：並查集 + 邊排序，O(E log E)，作為第三方基準。

思路：
    ### Prim 的貪心原理（cut property）
    維護已選點集 S。對任意跨越 S 與 V−S 的邊（稱為「割邊」），
    **權重最小的那條必定屬於某棵 MST**。Prim 每輪就是取這條最小割邊，
    把對應的新點併入 S，共做 V−1 輪。
    與 Dijkstra 的形式非常像：都是「每次取當前距離最小的未確定點」，
    差別在於 Dijkstra 的 dist 是「到起點的路徑總長」，而 Prim 的 dist 是「到當前樹的單邊權」。

    ### 樸素版 vs 堆優化版怎麼選
      - 樸素版每輪 O(V) 掃描找最小、O(V) 更新鄰點 → 總 O(V^2)，**與邊數無關**，
        所以稠密圖（E ≈ V^2）用它最划算，也不需要額外的堆結構。
      - 堆版每條邊最多進堆一次（懶刪除，可能有過期條目）→ O(E log V)，
        稀疏圖（E ≈ V）明顯更快。
      - 這裡兩個版本同時實現並斷言結果一致。

    ### 為什麼兩個版本的結果能逐字節一致
    兩者都採用同一套確定性平手規則：
      - 選點：先比 dist，相同則比 **點編號**（小者優先）；
      - 鬆弛：只有嚴格 `w < dist[v]` 才更新 parent（先來先佔）。
    堆版用 `(w, v)` 作為堆鍵，彈出時跳過已訪問的點，
    等價於「每次在未訪問且 dist 有限的點中取 (dist, 編號) 最小者」，與樸素版完全相同。

    ### 不連通圖
    Prim 從 0 號點出發只能覆蓋一個連通分量。偵測方式是「最終選中的點數 < V」。
    最小生成森林則是對每個未訪問的點再跑一次 Prim，把各分量的 MST 權重加總。

輸入格式（stdin，全部以空白分隔）：
    n m
    m 行：u v w      （無向邊，0-indexed，w 可為負）
輸出格式（stdout）：
    第 1 行：是否連通（1 / 0）
    第 2 行：prim_naive 的 MST 總權重（不連通為 -1）
    第 3 行：prim_heap  的 MST 總權重（不連通為 -1）
    第 4 行：kruskal    的 MST 總權重（不連通為 -1）
    第 5 行：最小生成森林的總權重（各分量 MST 權重之和）
    第 6 行：連通分量個數
    第 7 行：MST 邊數
    第 8 行起：每條 MST 邊一行 `u v w`（u < v，按 (u, v) 排序）
輸入被截斷時，有多少邊讀多少邊。無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import heapq
import random
import re
import sys
from typing import List, Optional, Tuple

INF = 10 ** 18

# 只接受 [+-]?digits，與 C++ 版本的 tryLL 完全一致；
# 其它 token 一律視為「輸入到此為止」，避免兩語言一個崩一個不崩。
_INT_RE = re.compile(r"^[+-]?[0-9]+$")


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；不是十進制整數則返回 None。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


# ---------------------------------------------------------------- 工具
class UnionFind:
    """並查集：路徑壓縮 + 按秩合併。"""

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.rank = [0] * n
        self.count = n          # 連通分量個數

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        self.count -= 1
        return True


def normalize_edges(n: int, edges: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
    """過濾掉非法邊（自環與越界），保留重邊（Prim / Kruskal 都能正確處理）。"""
    out = []
    for u, v, w in edges:
        if u == v:
            continue                       # 自環不可能是任何生成樹的邊
        if 0 <= u < n and 0 <= v < n:
            out.append((u, v, w))
    return out


def sort_mst_edges(mst: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
    """統一成 u < v 並按 (u, v) 排序，讓輸出確定。"""
    norm = [(min(u, v), max(u, v), w) for u, v, w in mst]
    norm.sort(key=lambda e: (e[0], e[1], e[2]))
    return norm


# ---------------------------------------------------------------- Prim：樸素 O(V^2)
def prim_naive(n: int, edges: List[Tuple[int, int, int]],
               root: int = 0) -> Tuple[bool, int, List[Tuple[int, int, int]]]:
    """鄰接矩陣 + 每輪線性掃描。適合稠密圖。返回 (是否連通, 總權重, MST 邊)。"""
    if n <= 0:
        return True, 0, []
    adj = [[INF] * n for _ in range(n)]
    for u, v, w in edges:
        if w < adj[u][v]:
            adj[u][v] = adj[v][u] = w

    used = [False] * n
    dist = [INF] * n
    parent = [-1] * n
    dist[root] = 0
    total = 0
    mst: List[Tuple[int, int, int]] = []
    picked = 0

    for _ in range(n):
        best = -1
        for v in range(n):
            if used[v] or dist[v] == INF:
                continue
            # 平手取編號小者，讓結果唯一
            if best < 0 or dist[v] < dist[best] or (dist[v] == dist[best] and v < best):
                best = v
        if best < 0:
            return False, -1, []            # 剩下的點都不可達
        used[best] = True
        picked += 1
        total += dist[best]
        if parent[best] >= 0:
            mst.append((parent[best], best, dist[best]))
        for v in range(n):
            if not used[v] and adj[best][v] < dist[v]:
                dist[v] = adj[best][v]
                parent[v] = best

    return (picked == n), total, sort_mst_edges(mst)


# ---------------------------------------------------------------- Prim：堆優化 O(E log V)
def prim_heap(n: int, edges: List[Tuple[int, int, int]],
              root: int = 0) -> Tuple[bool, int, List[Tuple[int, int, int]]]:
    """鄰接表 + 二元堆（懶刪除）。適合稀疏圖。返回 (是否連通, 總權重, MST 邊)。"""
    if n <= 0:
        return True, 0, []
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))

    used = [False] * n
    dist = [INF] * n
    parent = [-1] * n
    dist[root] = 0
    heap = [(0, root)]
    total = 0
    mst: List[Tuple[int, int, int]] = []
    picked = 0

    while heap:
        d, u = heapq.heappop(heap)
        if used[u]:
            continue                        # 過期條目（懶刪除）
        used[u] = True
        picked += 1
        total += d
        if parent[u] >= 0:
            mst.append((parent[u], u, d))
        for v, w in adj[u]:
            if not used[v] and w < dist[v]:
                dist[v] = w
                parent[v] = u
                heapq.heappush(heap, (w, v))

    return (picked == n), (total if picked == n else -1), (sort_mst_edges(mst) if picked == n else [])


# ---------------------------------------------------------------- Kruskal（第三方基準）
def kruskal(n: int, edges: List[Tuple[int, int, int]]) -> Tuple[bool, int, List[Tuple[int, int, int]]]:
    """並查集 + 邊排序。返回 (是否連通, 總權重, MST 邊)。"""
    if n <= 0:
        return True, 0, []
    ordered = sorted((w, min(u, v), max(u, v)) for u, v, w in edges)
    uf = UnionFind(n)
    total = 0
    mst: List[Tuple[int, int, int]] = []
    for w, u, v in ordered:
        if uf.union(u, v):
            total += w
            mst.append((u, v, w))
            if len(mst) == n - 1:
                break
    return (len(mst) == n - 1), (total if len(mst) == n - 1 else -1), mst


# ---------------------------------------------------------------- 連通分量 + 最小生成森林
def connected_components(n: int, edges: List[Tuple[int, int, int]]) -> int:
    uf = UnionFind(n)
    for u, v, _w in edges:
        uf.union(u, v)
    return uf.count


def minimum_spanning_forest(n: int, edges: List[Tuple[int, int, int]]) -> Tuple[int, int]:
    """返回 (森林總權重, 連通分量個數)。對每個分量各自跑一次堆版 Prim。"""
    comps = connected_components(n, edges)
    if n <= 0:
        return 0, 0
    visited = [False] * n
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))

    total = 0
    for s in range(n):
        if visited[s]:
            continue
        visited[s] = True
        heap = [(w, v) for v, w in adj[s]]
        heapq.heapify(heap)
        while heap:
            d, u = heapq.heappop(heap)
            if visited[u]:
                continue
            visited[u] = True
            total += d
            for v, w in adj[u]:
                if not visited[v]:
                    heapq.heappush(heap, (w, v))
    return total, comps


# ---------------------------------------------------------------- IO 模式
def run_io(data: str) -> None:
    toks = data.split()
    pos = 0

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos < len(toks):
            v = toks[pos]
            pos += 1
            return v
        return None

    def next_int(default: int) -> int:
        v = parse_int(nxt())
        return default if v is None else v

    n = max(next_int(0), 0)
    m = max(next_int(0), 0)
    raw: List[Tuple[int, int, int]] = []
    for _ in range(m):
        a, b, c = parse_int(nxt()), parse_int(nxt()), parse_int(nxt())
        if a is None or b is None or c is None:
            break                            # 輸入截斷或出現非數字 token
        raw.append((a, b, c))

    edges = normalize_edges(n, raw)
    ok_n, tot_n, mst_n = prim_naive(n, edges)
    ok_h, tot_h, mst_h = prim_heap(n, edges)
    ok_k, tot_k, _mst_k = kruskal(n, edges)
    forest, comps = minimum_spanning_forest(n, edges)

    print(1 if (ok_n and ok_h and ok_k) else 0)
    print(tot_n)
    print(tot_h)
    print(tot_k)
    print(forest)
    print(comps)
    print(len(mst_h))
    for u, v, w in mst_h:
        print(u, v, w)


# ---------------------------------------------------------------- 測試
def run_tests() -> None:
    # 空圖 / 單點
    assert prim_naive(0, []) == (True, 0, [])
    assert prim_heap(0, []) == (True, 0, [])
    assert kruskal(0, []) == (True, 0, [])
    assert prim_naive(1, []) == (True, 0, [])
    assert prim_heap(1, []) == (True, 0, [])
    assert kruskal(1, []) == (True, 0, [])

    # README 示例圖（0-1:2, 0-3:6, 1-2:3, 1-3:8, 1-4:5, 2-4:7, 3-4:9）
    edges = [(0, 1, 2), (0, 3, 6), (1, 2, 3), (1, 3, 8), (1, 4, 5), (2, 4, 7), (3, 4, 9)]
    ok_n, tot_n, mst_n = prim_naive(5, edges)
    ok_h, tot_h, mst_h = prim_heap(5, edges)
    ok_k, tot_k, _ = kruskal(5, edges)
    assert ok_n and ok_h and ok_k
    assert tot_n == tot_h == tot_k == 16
    assert mst_n == mst_h
    assert mst_h == [(0, 1, 2), (0, 3, 6), (1, 2, 3), (1, 4, 5)]
    assert connected_components(5, edges) == 1

    # 不連通：兩個三角形
    edges2 = [(0, 1, 1), (1, 2, 1), (0, 2, 5), (3, 4, 2), (4, 5, 2), (3, 5, 9)]
    assert prim_naive(6, edges2)[0] is False
    assert prim_naive(6, edges2)[1] == -1
    assert prim_heap(6, edges2)[:2] == (False, -1)
    assert kruskal(6, edges2)[:2] == (False, -1)
    forest, comps = minimum_spanning_forest(6, edges2)
    assert comps == 2
    assert forest == 2 + 4                 # 分量 A: 1+1=2；分量 B: 2+2=4

    # 自環與重邊：自環被忽略，重邊取小的
    edges3 = [(0, 0, -100), (0, 1, 4), (0, 1, 2), (1, 2, 3), (2, 2, 7)]
    assert prim_naive(3, normalize_edges(3, edges3)) == (True, 5, [(0, 1, 2), (1, 2, 3)])
    assert prim_heap(3, normalize_edges(3, edges3)) == (True, 5, [(0, 1, 2), (1, 2, 3)])
    assert kruskal(3, normalize_edges(3, edges3))[1] == 5

    # 負權邊：Prim / Kruskal 都照樣成立（最短路才怕負環）
    edges4 = [(0, 1, -5), (1, 2, -7), (0, 2, 100), (2, 3, -1)]
    assert prim_naive(4, edges4)[1] == -13
    assert prim_heap(4, edges4)[1] == -13
    assert kruskal(4, edges4)[1] == -13
    assert prim_naive(4, edges4)[2] == prim_heap(4, edges4)[2]

    # 孤立點（n 個點但邊為空）
    assert prim_naive(3, [])[0] is False
    assert prim_heap(3, [])[:2] == (False, -1)
    assert kruskal(3, [])[:2] == (False, -1)
    assert minimum_spanning_forest(3, []) == (0, 3)

    # 鏈狀圖：唯一生成樹
    chain = [(i, i + 1, i + 1) for i in range(9)]
    assert prim_naive(10, chain)[1] == sum(range(1, 10))
    assert prim_heap(10, chain)[1] == sum(range(1, 10))
    assert prim_naive(10, chain)[2] == prim_heap(10, chain)[2]

    # 完全圖 K6，邊權 = u+v（有大量平手，檢查確定性）
    full = [(u, v, u + v) for u in range(6) for v in range(u + 1, 6)]
    a_n, a_t, a_e = prim_naive(6, full)
    b_n, b_t, b_e = prim_heap(6, full)
    c_n, c_t, _ = kruskal(6, full)
    assert a_n and b_n and c_n
    assert a_t == b_t == c_t
    assert a_e == b_e                       # 平手規則一致 → 邊集也一致

    # 隨機對拍：三種實現的總權重必須一致（MST 權重唯一）
    random.seed(20261005)
    for _ in range(300):
        n = random.randint(1, 9)
        max_e = n * (n - 1) // 2
        m = random.randint(0, max_e)
        pool = [(u, v) for u in range(n) for v in range(u + 1, n)]
        random.shuffle(pool)
        es = [(u, v, random.randint(-9, 9)) for u, v in pool[:m]]
        ok_n, tot_n, mst_n = prim_naive(n, es)
        ok_h, tot_h, mst_h = prim_heap(n, es)
        ok_k, tot_k, _ = kruskal(n, es)
        assert ok_n == ok_h == ok_k
        assert tot_n == tot_h == tot_k
        assert mst_n == mst_h                # 平手規則一致，邊集也應相同
        assert len(mst_h) == (n - 1 if ok_n else 0)
        # 生成樹權重 = 最小生成森林權重（連通時）
        forest, comps = minimum_spanning_forest(n, es)
        if ok_n:
            assert comps == 1
            assert forest == tot_n
        else:
            assert comps > 1

    # 隨機對拍：最小生成森林 vs 對每個分量單獨加總
    for _ in range(120):
        n = random.randint(1, 10)
        m = random.randint(0, 12)
        es = []
        for _ in range(m):
            u = random.randint(0, n - 1)
            v = random.randint(0, n - 1)
            es.append((u, v, random.randint(-9, 9)))
        es = normalize_edges(n, es)
        forest, comps = minimum_spanning_forest(n, es)
        assert comps == connected_components(n, es)
        if comps == 1:
            assert forest == kruskal(n, es)[1]

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
