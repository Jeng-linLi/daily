"""最小費用最大流（Min-Cost Max-Flow, MCMF）

題意：
    給一張有向圖，每條邊有容量 cap 與單位費用 cost。
    在所有「從 s 到 t 流量最大」的流中，找出**總費用最小**的那一條，並輸出每條邊的流量。

思路：
    ### 連續最短路增廣（Successive Shortest Path, SSP）
    最短路增廣是最大流 Ford-Fulkerson 的推廣：每次不在「任意殘量路徑」上增廣，
    而是在「**費用最小**的殘量路徑」上增廣，能增廣多少就增廣多少（直到該路徑上某條邊飽和）。

    > **定理**：若每次都沿當前殘量網路中的最短路（以費用為邊權）增廣，
    > 則增廣過程中得到的每個流都是「同流量下費用最小」的流。
    > 因此算法結束時得到的就是最小費用最大流。

    直觀上，若某個流量 f 的流不是最小費用，它與真正的最小費用流之差可以分解成若干個
    「負費用環 + 路」，沿最短路增廣恰好不會引入這種可改進的結構。

    ### 殘量網路的反向邊
    正向邊 (u→v, cap, cost) 會配一條反向邊 (v→u, 0, -cost)。
    增廣時正向減、反向加；費用計算用「原始邊」的 cost。
    反向邊的負費用是它能「退流」的關鍵，也讓殘量網路出現負權邊。

    ### 兩個版本的差別
    | 版本 | 每次求最短路 | 適用 |
    |---|---|---|
    | SPFA / Bellman-Ford | O(VE) | 有負權邊（反向邊必然負權），寫起來最直接 |
    | Dijkstra + Johnson 勢能 | O(E log V) | 先用一次 SPFA 算出初始勢能 h，把邊權改為非負 |

    勢能的技巧：令初始 `h[v]` = 從 s 出發的最短路（SPFA 算一次，可處理負邊）。
    之後每次增廣後更新 `h[v] += dist[v]`，並用**既約費用**
        `c'(u,v) = cost(u,v) + h[u] - h[v] ≥ 0`
    當邊權跑 Dijkstra。這就是 Johnson 全點對最短路的同一套手法。
    兩個版本在本題中互相對拍，答案必須一致。

    ### 應用：指派問題（Assignment Problem）
    n 個人做 n 件工作，每人一件、每件一人，`c[i][j]` 是 i 做 j 的費用。
    建模：s → 每個人（cap 1, cost 0）、人 → 工作（cap 1, cost c[i][j]）、每件工作 → t（cap 1, cost 0）。
    最大流必為 n，其最小費用即最小總費用指派。匈牙利算法是它的特例（O(n³)），
    費用流版本是 O(流量 · E log V)，但能處理「一人可做多件」「有容量分組」等更一般的情況。

應用場景：
    - 物流與供應鏈：多倉到多店的運輸排程，在滿足需求的同時最小化運費
    - 雲端排程：把任務指派到機器，最小化成本或延遲懲罰
    - 電力調度、網路流量工程（TE）的最小成本路由
    - 二分圖帶權匹配：人員派工、廣告投放預算分配、網約車派單
    - 影像分割（graph cuts 的費用流變體）、庫存與生產計畫的滾動排程

複雜度：
    SSP + SPFA      O(流量 · V · E)，流量 ≤ Σ cap
    SSP + Dijkstra  O(流量 · E log V)
    空間            O(V + E)
    （若流量以二進位位數計，嚴格上界是 O(F · E log V)，F 為最大流值）

輸入格式（stdin，全部以空白分隔）：
    n m s t
    u v cap cost     接下來 m 條邊（u、v 為 0-indexed 點號）

    非數字 token 一律視為「輸入到此為止」。
    數值截斷：n ≤ 100、m ≤ 400、cap ∈ [0, 10^6]、cost ∈ [-10^6, 10^6]；
    點號越界的邊視為無效（流量記 0）。
輸出格式（stdout，共 3 行）：
    第 1 行  flow cost       SPFA 版的最大流值與最小費用
    第 2 行  m 個整數         每條輸入邊上的流量（按輸入順序）
    第 3 行  flow2 cost2     Dijkstra + 勢能版（應與第 1 行相同）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import heapq
import itertools
import random
import re
import sys
from collections import deque
from typing import List, Optional, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 100
MAXM = 400
CAP_MAX = 10 ** 6
COST_ABS = 10 ** 6
INF = 1 << 60            # 整數哨兵：避免浮點污染費用輸出（兩個語言一致）


class Edge:
    """殘量網路中的一條邊；`index ^ 1` 是它的反向邊。"""

    __slots__ = ("to", "cap", "cost")

    def __init__(self, to: int, cap: int, cost: int) -> None:
        self.to = to
        self.cap = cap
        self.cost = cost


class MinCostFlow:
    """最小費用最大流（兩種最短路增廣實作）。"""

    __slots__ = ("n", "edges", "adj")

    def __init__(self, n: int) -> None:
        self.n = n
        self.edges: List[Edge] = []
        self.adj: List[List[int]] = [[] for _ in range(n)]

    def add_edge(self, u: int, v: int, cap: int, cost: int) -> int:
        """加入一條邊，回傳正向邊在 `edges` 中的索引。"""
        idx = len(self.edges)
        self.edges.append(Edge(v, cap, cost))
        self.edges.append(Edge(u, 0, -cost))
        self.adj[u].append(idx)
        self.adj[v].append(idx ^ 1)
        return idx

    # ------------------------------------------------------------ SPFA 版
    def _spfa(self, s: int) -> Tuple[List[int], List[int], List[int]]:
        """從 s 出發的最短路（允許負權邊），回傳 (dist, prevv, preve)。"""
        n = self.n
        dist = [INF] * n
        prevv = [-1] * n
        preve = [-1] * n
        inq = [False] * n
        dist[s] = 0
        dq = deque([s])
        inq[s] = True
        while dq:
            u = dq.popleft()
            inq[u] = False
            du = dist[u]
            for ei in self.adj[u]:
                e = self.edges[ei]
                if e.cap > 0 and du + e.cost < dist[e.to]:
                    dist[e.to] = du + e.cost
                    prevv[e.to] = u
                    preve[e.to] = ei
                    if not inq[e.to]:
                        inq[e.to] = True
                        dq.append(e.to)
        return dist, prevv, preve

    def has_negative_cycle(self, s: int) -> bool:
        """Bellman-Ford 負費用環檢測（只看從 s 可達的子圖）。

        初始圖若含負環，增廣路徑的費用可以任意小（沿負環繞幾圈），
        最小費用沒有下界，SSP 也不會終止。本題偵測到此情況時一律回傳 (0, 0)，
        兩個語言用同一套判定，輸出必然一致。
        """
        n = self.n
        dist = [INF] * n
        dist[s] = 0
        for _ in range(max(0, n - 1)):
            changed = False
            for u in range(n):
                if dist[u] == INF:
                    continue
                for ei in self.adj[u]:
                    e = self.edges[ei]
                    if e.cap > 0 and dist[u] + e.cost < dist[e.to]:
                        dist[e.to] = dist[u] + e.cost
                        changed = True
            if not changed:
                break
        for u in range(n):
            if dist[u] == INF:
                continue
            for ei in self.adj[u]:
                e = self.edges[ei]
                if e.cap > 0 and dist[u] + e.cost < dist[e.to]:
                    return True
        return False

    def min_cost_flow_spfa(self, s: int, t: int, maxf: int = 10 ** 18) -> Tuple[int, int]:
        """SPFA（Bellman-Ford 佇列優化）版 SSP；回傳 (流量, 最小費用)。

        注意：s == t 時流量無上界，本題統一定義為 0（兩個語言一致）。
        """
        if s == t:
            return 0, 0
        if self.has_negative_cycle(s):
            return 0, 0
        flow = 0
        cost = 0
        while flow < maxf:
            dist, prevv, preve = self._spfa(s)
            if dist[t] == INF:
                break
            # 一次推滿整條路徑（阻塞流），直到某條邊飽和
            d = maxf - flow
            v = t
            while v != s:
                d = min(d, self.edges[preve[v]].cap)
                v = prevv[v]
            flow += d
            cost += d * dist[t]
            v = t
            while v != s:
                ei = preve[v]
                self.edges[ei].cap -= d
                self.edges[ei ^ 1].cap += d
                v = prevv[v]
        return flow, cost

    # ------------------------------------------------------------ Dijkstra + 勢能版
    def min_cost_flow_dijkstra(self, s: int, t: int, maxf: int = 10 ** 18) -> Tuple[int, int]:
        """Dijkstra + Johnson 勢能版；要求初始圖無負費用環。s == t 時回傳 (0, 0)。"""
        if s == t:
            return 0, 0
        if self.has_negative_cycle(s):
            return 0, 0
        n = self.n
        h, _, _ = self._spfa(s)                    # 初始勢能（可處理負邊）
        flow = 0
        cost = 0
        prevv = [-1] * n
        preve = [-1] * n
        while flow < maxf:
            dist = [INF] * n
            dist[s] = 0
            pq: List[Tuple[int, int]] = [(0, s)]
            while pq:
                d, u = heapq.heappop(pq)
                if d > dist[u]:
                    continue
                hu = h[u]
                for ei in self.adj[u]:
                    e = self.edges[ei]
                    if e.cap <= 0 or hu == INF or h[e.to] == INF:
                        continue
                    nd = d + e.cost + hu - h[e.to]  # 既約費用，必 ≥ 0
                    if nd < dist[e.to]:
                        dist[e.to] = nd
                        prevv[e.to] = u
                        preve[e.to] = ei
                        heapq.heappush(pq, (nd, e.to))
            if dist[t] == INF:
                break
            for v in range(n):                     # 更新勢能，維持既約費用的非負性
                if dist[v] < INF:
                    h[v] += dist[v]
            d = maxf - flow
            v = t
            while v != s:
                d = min(d, self.edges[preve[v]].cap)
                v = prevv[v]
            flow += d
            cost += d * (h[t] - h[s])              # 真實最短路長 = h[t] - h[s]
            v = t
            while v != s:
                ei = preve[v]
                self.edges[ei].cap -= d
                self.edges[ei ^ 1].cap += d
                v = prevv[v]
        return flow, cost


# ---------------------------------------------------------------- 最大流（對拍用）


class Dinic:
    """最大流（不考慮費用），用來對拍流量值是否相同。"""

    def __init__(self, n: int) -> None:
        self.n = n
        self.graph: List[List[int]] = [[] for _ in range(n)]
        self.cap: List[int] = []
        self.to: List[int] = []

    def add_edge(self, u: int, v: int, c: int) -> None:
        self.to.append(v)
        self.cap.append(c)
        self.graph[u].append(len(self.to) - 1)
        self.to.append(u)
        self.cap.append(0)
        self.graph[v].append(len(self.to) - 1)

    def max_flow(self, s: int, t: int) -> int:
        if s == t:
            return 0
        flow = 0
        while True:
            level = [-1] * self.n
            level[s] = 0
            q = deque([s])
            while q:
                u = q.popleft()
                for ei in self.graph[u]:
                    if self.cap[ei] > 0 and level[self.to[ei]] < 0:
                        level[self.to[ei]] = level[u] + 1
                        q.append(self.to[ei])
            if level[t] < 0:
                return flow
            it = [0] * self.n
            while True:
                f = self._dfs(s, t, 10 ** 18, level, it)
                if f == 0:
                    break
                flow += f

    def _dfs(self, u: int, t: int, f: int, level: List[int], it: List[int]) -> int:
        if u == t:
            return f
        while it[u] < len(self.graph[u]):
            ei = self.graph[u][it[u]]
            v = self.to[ei]
            if self.cap[ei] > 0 and level[v] == level[u] + 1:
                d = self._dfs(v, t, min(f, self.cap[ei]), level, it)
                if d > 0:
                    self.cap[ei] -= d
                    self.cap[ei ^ 1] += d
                    return d
            it[u] += 1
        return 0


def edmonds_karp_max_flow(n: int, edges: Sequence[Tuple[int, int, int]], s: int, t: int) -> int:
    """BFS 增廣的最大流（O(VE²)），僅供對拍。s == t 時回傳 0。"""
    if s == t:
        return 0
    cap = [[0] * n for _ in range(n)]
    for u, v, c in edges:
        cap[u][v] += c
    flow = 0
    while True:
        prev = [-1] * n
        prev[s] = s
        q = deque([s])
        while q and prev[t] == -1:
            u = q.popleft()
            for v in range(n):
                if prev[v] == -1 and cap[u][v] > 0:
                    prev[v] = u
                    q.append(v)
        if prev[t] == -1:
            return flow
        d = 10 ** 18
        v = t
        while v != s:
            d = min(d, cap[prev[v]][v])
            v = prev[v]
        v = t
        while v != s:
            cap[prev[v]][v] -= d
            cap[v][prev[v]] += d
            v = prev[v]
        flow += d


# ---------------------------------------------------------------- 暴力對拍


def reachable_from(n: int, edges: Sequence[Tuple[int, int, int, int]], s: int) -> List[bool]:
    """從 s 沿「容量 > 0」的邊能到達的點（費用流只會在這些點上推流）。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v, c, _ in edges:
        if c > 0:
            adj[u].append(v)
    seen = [False] * n
    seen[s] = True
    dq = deque([s])
    while dq:
        u = dq.popleft()
        for v in adj[u]:
            if not seen[v]:
                seen[v] = True
                dq.append(v)
    return seen


def bf_min_cost_max_flow(n: int, edges: Sequence[Tuple[int, int, int, int]],
                         s: int, t: int) -> Tuple[int, int]:
    """枚舉所有邊流量組合，找「最大流量下的最小費用」，僅供對拍（規模極小）。

    只允許在「從 s 可達」的點上推流：費用流不會跑去 s 到不了的負費用環上繞。
    """
    m = len(edges)
    if m == 0:
        return 0, 0
    reach = reachable_from(n, edges, s)
    maxflow = edmonds_karp_max_flow(n, [(u, v, c) for u, v, c, _ in edges], s, t)
    best = None
    ranges = [range(c + 1) if reach[u] else range(1) for u, _, c, _ in edges]
    for combo in itertools.product(*ranges):
        net = [0] * n
        total = 0
        for i, (u, v, _, w) in enumerate(edges):
            f = combo[i]
            net[u] -= f
            net[v] += f
            total += f * w
        if any(net[v] != 0 for v in range(n) if v != s and v != t):
            continue
        if -net[s] != maxflow:
            continue
        if best is None or total < best:
            best = total
    return maxflow, (0 if best is None else best)


# ---------------------------------------------------------------- 應用：指派問題


def min_cost_assignment(cost: Sequence[Sequence[int]]) -> Tuple[int, List[int]]:
    """n 人 n 工作的最小費用完美指派；回傳 (最小總費用, 每個人被分到的工作)。"""
    n = len(cost)
    if n == 0:
        return 0, []
    s, t = 2 * n, 2 * n + 1
    mcf = MinCostFlow(2 * n + 2)
    for i in range(n):
        mcf.add_edge(s, i, 1, 0)
    for i in range(n):
        for j in range(n):
            mcf.add_edge(i, n + j, 1, cost[i][j])
    for j in range(n):
        mcf.add_edge(n + j, t, 1, 0)
    flow, total = mcf.min_cost_flow_spfa(s, t)
    assert flow == n
    who = [-1] * n
    for i in range(n):
        for ei in mcf.adj[i]:
            e = mcf.edges[ei]
            v = e.to
            if n <= v < 2 * n and mcf.edges[ei ^ 1].cap > 0:   # 反向邊有容量 = 正向有流量
                who[i] = v - n
    return total, who


def bf_assignment(cost: Sequence[Sequence[int]]) -> int:
    """暴力枚舉所有排列，僅供對拍（n ≤ 7）。"""
    n = len(cost)
    if n == 0:
        return 0
    return min(sum(cost[i][p[i]] for i in range(n)) for p in itertools.permutations(range(n)))


# ---------------------------------------------------------------- 輸出


def solve(n: int, edges: Sequence[Tuple[int, int, int, int]], s: int, t: int) -> str:
    g1 = MinCostFlow(n)
    g2 = MinCostFlow(n)
    idx1: List[Optional[int]] = []
    idx2: List[Optional[int]] = []
    for u, v, cap, w in edges:
        idx1.append(g1.add_edge(u, v, cap, w))
        idx2.append(g2.add_edge(u, v, cap, w))
    f1, c1 = g1.min_cost_flow_spfa(s, t)
    f2, c2 = g2.min_cost_flow_dijkstra(s, t)
    flows = []
    for i, (_, _, cap, _) in enumerate(edges):
        e = g1.edges[idx1[i]]
        flows.append(str(cap - e.cap) if cap >= 0 else "0")
    return "\n".join([
        f"{f1} {c1}",
        " ".join(flows),
        f"{f2} {c2}",
    ])


# ---------------------------------------------------------------- 測試


def check_graph(n: int, edges: List[Tuple[int, int, int, int]], s: int, t: int) -> None:
    """兩種 MCMF 一致、流量等於最大流、費用等於暴力最小值。

    初始圖若含負費用環則最小費用無下界，只驗證兩個版本都回傳 (0, 0)。
    """
    g1 = MinCostFlow(n)
    g2 = MinCostFlow(n)
    for u, v, c, w in edges:
        g1.add_edge(u, v, c, w)
        g2.add_edge(u, v, c, w)
    if n > 0 and s != t and g1.has_negative_cycle(s):
        assert g1.min_cost_flow_spfa(s, t) == (0, 0)
        assert g2.min_cost_flow_dijkstra(s, t) == (0, 0)
        return
    f1, c1 = g1.min_cost_flow_spfa(s, t)
    f2, c2 = g2.min_cost_flow_dijkstra(s, t)
    assert f1 == f2, (n, edges, s, t)
    assert c1 == c2, (n, edges, s, t)

    din = Dinic(n)
    for u, v, c, _ in edges:
        din.add_edge(u, v, c)
    assert f1 == din.max_flow(s, t), (n, edges, s, t)

    if s == t:                      # s == t 時流量無上界，本題統一定義為 0，不與暴力比對
        assert f1 == 0 and c1 == 0
        return

    bf_f, bf_c = bf_min_cost_max_flow(n, edges, s, t)
    assert f1 == bf_f and c1 == bf_c, (n, edges, s, t, f1, c1, bf_f, bf_c)


def run_tests() -> None:
    rng = random.Random(20261013)

    # ---- 固定用例（先跑過再寫死期望值）----
    assert solve(4, [(0, 1, 3, 1), (0, 2, 2, 4), (1, 3, 2, 2), (2, 3, 3, 1)], 0, 3) \
        == "4 16\n2 2 2 2\n4 16"
    assert solve(2, [(0, 1, 5, 3)], 0, 1) == "5 15\n5\n5 15"
    assert solve(3, [], 0, 2) == "0 0\n\n0 0"
    assert solve(2, [(1, 0, 4, -2)], 0, 1) == "0 0\n0\n0 0"

    check_graph(4, [(0, 1, 3, 1), (0, 2, 2, 4), (1, 3, 2, 2), (2, 3, 3, 1)], 0, 3)
    check_graph(2, [(0, 1, 5, 3)], 0, 1)
    check_graph(3, [], 0, 2)
    check_graph(2, [(1, 0, 4, -2)], 0, 1)

    # ---- 隨機對拍：小圖（費用非負，不會有負環），與暴力枚舉比對 ----
    for _ in range(400):
        n = rng.randint(2, 5)
        m = rng.randint(0, 5)
        edges = []
        for _ in range(m):
            u = rng.randrange(n)
            v = rng.randrange(n)
            edges.append((u, v, rng.randint(0, 3), rng.randint(0, 6)))
        s, t = rng.randrange(n), rng.randrange(n)
        check_graph(n, edges, s, t)

    # ---- 隨機對拍：負費用邊（邊一律 u < v，圖無環 ⟹ 保證無負環）----
    for _ in range(200):
        n = rng.randint(3, 6)
        m = rng.randint(1, 6)
        edges = []
        for _ in range(m):
            u = rng.randrange(n - 1)
            v = rng.randint(u + 1, n - 1)
            edges.append((u, v, rng.randint(0, 3), rng.choice([-3, -1, 0, 2, 5])))
        s, t = rng.randrange(n), rng.randrange(n)
        check_graph(n, edges, s, t)

    # ---- 隨機對拍：任意邊 + 任意費用（可能含負環，只驗證兩個版本一致）----
    for _ in range(150):
        n = rng.randint(2, 5)
        m = rng.randint(0, 5)
        edges = []
        for _ in range(m):
            u = rng.randrange(n)
            v = rng.randrange(n)
            edges.append((u, v, rng.randint(0, 3), rng.randint(-4, 6)))
        s, t = rng.randrange(n), rng.randrange(n)
        check_graph(n, edges, s, t)

    # ---- 隨機對拍：中圖，只比對兩個版本 + Dinic ----
    for _ in range(120):
        n = rng.randint(2, 20)
        m = rng.randint(0, 40)
        edges = []
        for _ in range(m):
            u = rng.randrange(n)
            v = rng.randrange(n)
            edges.append((u, v, rng.randint(0, 8), rng.randint(-10, 10)))
        s, t = rng.randrange(n), rng.randrange(n)
        g1 = MinCostFlow(n)
        g2 = MinCostFlow(n)
        for u, v, c, w in edges:
            g1.add_edge(u, v, c, w)
            g2.add_edge(u, v, c, w)
        if n > 0 and s != t and g1.has_negative_cycle(s):
            assert g1.min_cost_flow_spfa(s, t) == (0, 0)
            assert g2.min_cost_flow_dijkstra(s, t) == (0, 0)
            continue
        f1, c1 = g1.min_cost_flow_spfa(s, t)
        f2, c2 = g2.min_cost_flow_dijkstra(s, t)
        assert f1 == f2 and c1 == c2, (n, edges, s, t)
        din = Dinic(n)
        for u, v, c, _ in edges:
            din.add_edge(u, v, c)
        assert f1 == din.max_flow(s, t), (n, edges, s, t)
        # 費用恆等式：總費用 = Σ 正向邊流量 × 單位費用
        total = 0
        for i, (_, _, _, w) in enumerate(edges):
            total += (edges[i][2] - g1.edges[2 * i].cap) * w
        assert total == c1, (n, edges, s, t)

    # ---- 應用：指派問題對拍 ----
    for _ in range(200):
        n = rng.randint(1, 6)
        cost = [[rng.randint(-5, 20) for _ in range(n)] for _ in range(n)]
        total, who = min_cost_assignment(cost)
        assert total == bf_assignment(cost), cost
        assert sorted(who) == list(range(n)), (cost, who)
        assert total == sum(cost[i][who[i]] for i in range(n))

    # 經典範例：對角線最便宜
    assert min_cost_assignment([[4, 1, 3], [2, 0, 5], [3, 2, 2]])[0] == 5
    assert min_cost_assignment([[10]])[0] == 10
    assert min_cost_assignment([])[0] == 0


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。"""
    tokens = raw.split()

    def read_int(i: int) -> Optional[int]:
        return int(tokens[i]) if i < len(tokens) and NUM_RE.match(tokens[i]) else None

    def clamp(v: Optional[int], lo: int, hi: int, default: int) -> int:
        if v is None:
            return default
        return max(lo, min(hi, v))

    n = clamp(read_int(0), 0, MAXN, 0)
    m = clamp(read_int(1), 0, MAXM, 0)
    s = clamp(read_int(2), 0, MAXN - 1, 0)
    t = clamp(read_int(3), 0, MAXN - 1, 0)

    edges: List[Tuple[int, int, int, int]] = []
    valid: List[bool] = []
    for i in range(m):
        u = read_int(4 + 4 * i)
        v = read_int(5 + 4 * i)
        cap = read_int(6 + 4 * i)
        w = read_int(7 + 4 * i)
        if u is None or v is None or cap is None or w is None:
            break
        ok = 0 <= u < n and 0 <= v < n
        valid.append(ok)
        if not ok:
            continue
        edges.append((u, v, clamp(cap, 0, CAP_MAX, 0), clamp(w, -COST_ABS, COST_ABS, 0)))

    if n == 0 or s >= n or t >= n:
        out = "0 0\n" + " ".join(["0"] * len(valid)) + "\n0 0"
        emit(out + "\n")
        return

    emit(solve(n, edges, s, t) + "\n")


def emit(text: str) -> None:
    """以二進位輸出（latin-1），避免 Windows 把 \n 換成 \r\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
