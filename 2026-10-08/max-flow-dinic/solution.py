"""Dinic 最大流與二分圖匹配（最小割 / 邊不相交路徑 / 匈牙利算法）

題意：
    給定一張**有向圖**（邊容量為非負整數）、源點 s、匯點 t，以及一張二分圖，求：
      1. **最大流**：從 s 到 t 最多能送多少流量（每條邊不超過其容量、每個中間點流入=流出）；
      2. **最小割**：把點分成含 s 的 S 與含 t 的 T，割容量 = Σ 從 S 跨到 T 的邊容量，求最小者
         （並給出源側點集與具體割邊）；
      3. **二分圖最大匹配**：左部 nl 個點、右部 nr 個點，求最多能配多少對不相交的邊
         （匈牙利 / Kuhn 算法，並用最大流建模交叉驗證）；
      4. **邊不相交路徑數**：無向圖中從 s 到 t 最多能找出多少條「沒有共用邊」的路徑。

思路：
    ### 為什麼 Ford-Fulkerson 不夠
    最樸素的增廣路算法每次隨便找一條路，容量是大整數時可能要增廣 O(flow) 次。
    Dinic 的改進是**分層**：每輪先用 BFS 只在「殘量 > 0」的邊上建出層次圖（level[v] = level[u] + 1），
    再用 DFS 在層次圖上把**阻塞流**一次推完（不再有 s→t 的層次路徑為止），然後重新 BFS。
    每輪至少讓 level[t] 加 1，故輪數 O(V)，每輪 O(VE) → 總計 **O(V²E)**；
    對二分圖匹配這類單位容量圖，實際表現是 O(E√V) 級別，非常快。

    ### 反向邊（演算法的靈魂）
    每條邊 u→v（容量 c）同時建一條反向邊 v→u（容量 0）。推送 f 時：正向邊 -f、反向邊 +f。
    反向邊代表「反悔」——把先前流過去的流量退回去，等價於重新安排路徑。
    少了它，貪心選路的錯誤就無法修正，也就拿不到最大流。

    ### 最小割 = 最大流
    最大流跑完後，在**殘量網絡**中從 s 做一次 BFS，可達點集記為 S。
    此時任何從 S 跨到 T 的邊都必定已滿流（否則對面也可達），
    而任何從 T 跨回 S 的邊必定零流（否則反向有殘量、對面也可達），
    所以割容量恰好等於流量 → **最大流最小割定理**。源側點集也能直接讀出來。

    ### 二分圖匹配的兩種做法
      - **匈牙利（Kuhn）**：逐個左點嘗試增廣，遇到已佔用的右點就遞迴「讓位」
        （對該右點當前匹配到的左點重新找對象），O(VE) 但常數極小，且能直接輸出匹配邊；
      - **最大流建模**：超級源 S 連每個左點（容量 1）、左→右連邊（容量 1）、
        每個右點連超級匯 T（容量 1），最大流即最大匹配數。容量為 1 保證一個點只被用一次。

    ### 邊不相交路徑
    把每條無向邊拆成兩條容量 1 的有向邊，最大流就是邊不相交路徑數（Menger 定理的流版本）。
    若題目要求「點不相交」，則把每個點拆成 in/out 兩點、中間連容量 1 的邊即可。

    ### 確定性（保證 Python 與 C++ 輸出逐字節一致）
      - BFS 按鄰接表順序擴展；DFS 嚴格按 `it[u]` 遊標推進；
      - 割邊按 (u, v) 升序輸出，源側點集升序；
      - 匈牙利按左點升序增廣，匹配邊按左點升序輸出。

應用場景：
    交通 / 頻寬網路的吞吐上限（最大流）、網路可靠性的最小代價切割（最小割）、
    任務與工人的分配 / 課程選課配額（二分圖匹配）、機場跑道與航線規劃、
    影像分割（graph cuts 把像素當點、相鄰懲罰當邊）、專案選擇的最大獲利（最小割對偶）。

複雜度：
    Dinic 最大流      O(V²E) 時間（單位容量圖實務上遠快於此）、O(V + E) 空間
    最小割（一次 BFS） O(E) 時間
    匈牙利匹配        O(nl · E) 時間、O(nr) 空間
    Edmonds-Karp 對拍 O(V E²) 時間

輸入格式（stdin，全部以空白分隔）：
    n m s t
    u1 v1 c1
    ...（共 m 行）
    nl nr k
    a1 b1
    ...（共 k 行，左點編號 0..nl-1，右點編號 0..nr-1）
    讀不到時缺的部分補 0；遇到非整數 token 視為輸入結束。
    n 上限 200、m 上限 5000；越界的邊會被忽略；負容量視為 0。
輸出格式（stdout）：
    第 1 行：最大流值
    第 2 行：最小割容量
    第 3 行：源側點集大小 |S|
    第 4 行：源側點集（升序，空格分隔；空則為空行）
    第 5 行：割邊條數 e
    第 6 .. 5+e 行：每條割邊一行 `u v c`（按 u、v 升序）
    第 6+e 行：二分圖最大匹配數
    第 7+e 行起：每條匹配邊一行 `u v`（按左點升序）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from collections import deque
from typing import List, Optional, Set, Tuple

sys.setrecursionlimit(1000000)

INF = 4 * 10 ** 18                           # 大於任何可能的流量，且不溢出 int64
CAP_MAX = 10 ** 15                           # 容量上限：保證 Σ 流量仍在 int64 內（與 C++ 一致）
MAX_N = 200                                  # 防止 IO 模式下圖太大
MAX_M = 5000
_INT_RE = re.compile(r"^[+-]?[0-9]+$")      # 嚴格整數規則，與 C++ 的 tryLL 完全一致


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；非十進制整數則返回 None（視為輸入到此為止）。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


# ---------------------------------------------------------------- Dinic


class Dinic:
    """Dinic 最大流。每條邊存 [to, cap, rev]，rev 是反向邊在其端點鄰接表中的下標。"""

    def __init__(self, n: int) -> None:
        self.n = n
        self.g: List[List[List[int]]] = [[] for _ in range(n)]
        self.edges: List[Tuple[int, int, int]] = []      # 原始邊 (u, v, c)

    def add_edge(self, u: int, v: int, c: int) -> int:
        """加入一條容量為 c 的有向邊（自動加反向邊），返回邊的編號，越界則返回 -1。"""
        if not (0 <= u < self.n and 0 <= v < self.n):
            return -1
        if c < 0:
            c = 0
        self.g[u].append([v, c, len(self.g[v])])
        self.g[v].append([u, 0, len(self.g[u]) - 1])
        self.edges.append((u, v, c))
        return len(self.edges) - 1

    def _bfs(self, s: int, t: int) -> List[int]:
        """在殘量網絡上分層，返回 level；level[t] < 0 表示已無增廣路。"""
        level = [-1] * self.n
        level[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for e in self.g[u]:
                v = e[0]
                if e[1] > 0 and level[v] < 0:
                    level[v] = level[u] + 1
                    q.append(v)
        return level

    def _blocking(self, s: int, t: int, level: List[int], it: List[int]) -> int:
        """在層次圖上找一條增廣路並推送；返回推送量，0 表示本輪阻塞流已推完。"""
        stack = [s]
        path: List[Tuple[int, int]] = []
        while stack:
            u = stack[-1]
            if u == t:
                f = INF
                for (pu, pi) in path:
                    cap = self.g[pu][pi][1]
                    if cap < f:
                        f = cap
                for (pu, pi) in path:
                    e = self.g[pu][pi]
                    e[1] -= f
                    self.g[e[0]][e[2]][1] += f       # 反向邊 +f，保留「反悔」的餘地
                return f
            advanced = False
            while it[u] < len(self.g[u]):
                i = it[u]
                e = self.g[u][i]
                v = e[0]
                if e[1] > 0 and level[v] == level[u] + 1:
                    path.append((u, i))
                    stack.append(v)
                    advanced = True
                    break
                it[u] += 1
            if not advanced:
                level[u] = -1                        # 該點在層次圖上走不到 t，剪掉
                stack.pop()
                if path:
                    pu, pi = path.pop()
                    it[pu] += 1
        return 0

    def max_flow(self, s: int, t: int) -> int:
        """求 s → t 的最大流；s == t 或端點越界時返回 0。"""
        if not (0 <= s < self.n and 0 <= t < self.n) or s == t:
            return 0
        total = 0
        while True:
            level = self._bfs(s, t)
            if level[t] < 0:
                break
            it = [0] * self.n
            while True:
                f = self._blocking(s, t, level, it)
                if f == 0:
                    break
                total += f
        return total

    def reachable_from(self, s: int) -> List[int]:
        """殘量網絡中從 s 可達的點（升序），即最小割的源側 S。"""
        seen = [False] * self.n
        if 0 <= s < self.n:
            seen[s] = True
            q = deque([s])
            while q:
                u = q.popleft()
                for e in self.g[u]:
                    v = e[0]
                    if e[1] > 0 and not seen[v]:
                        seen[v] = True
                        q.append(v)
        return [i for i in range(self.n) if seen[i]]

    def min_cut(self, s: int) -> Tuple[int, List[int], List[Tuple[int, int, int]]]:
        """最大流跑完後求最小割，返回 (割容量, 源側點集升冪, 割邊[(u,v,c)] 按 (u,v) 升冪)。"""
        side: Set[int] = set(self.reachable_from(s))
        cap = 0
        cuts: List[Tuple[int, int, int]] = []
        for (u, v, c) in self.edges:
            if u in side and v not in side:
                cap += c
                cuts.append((u, v, c))
        cuts.sort(key=lambda x: (x[0], x[1]))
        return cap, sorted(side), cuts


def build_flow(n: int, edges: List[Tuple[int, int, int]]) -> Dinic:
    """把邊表建成流網絡（過濾越界邊、負容量歸零）。"""
    net = Dinic(n)
    for (u, v, c) in edges:
        net.add_edge(u, v, c)
    return net


# ---------------------------------------------------------------- 對拍用的其它算法


def edmonds_karp(n: int, edges: List[Tuple[int, int, int]], s: int, t: int) -> int:
    """Edmonds-Karp：每次用 BFS 找最短增廣路，O(V E²)，作為 Dinic 的獨立對拍。"""
    if n <= 0 or not (0 <= s < n and 0 <= t < n) or s == t:
        return 0
    cap = [[0] * n for _ in range(n)]
    for (u, v, c) in edges:
        if 0 <= u < n and 0 <= v < n and c > 0:
            cap[u][v] += c
    flow = 0
    while True:
        pre = [-1] * n
        pre[s] = -2
        q = deque([s])
        while q and pre[t] == -1:
            u = q.popleft()
            for v in range(n):
                if pre[v] == -1 and cap[u][v] > 0:
                    pre[v] = u
                    q.append(v)
        if pre[t] == -1:
            break
        f = INF
        v = t
        while v != s:
            f = min(f, cap[pre[v]][v])
            v = pre[v]
        v = t
        while v != s:
            p = pre[v]
            cap[p][v] -= f
            cap[v][p] += f
            v = p
        flow += f
    return flow


def min_cut_bruteforce(n: int, edges: List[Tuple[int, int, int]], s: int, t: int) -> int:
    """枚舉所有「含 s 不含 t」的點集求最小割容量（n ≤ 14），驗證最大流最小割定理。"""
    if n <= 0 or not (0 <= s < n and 0 <= t < n) or s == t:
        return 0
    best = INF
    for mask in range(1 << n):
        if not (mask >> s & 1) or (mask >> t & 1):
            continue
        total = 0
        for (u, v, c) in edges:
            if (mask >> u & 1) and not (mask >> v & 1):
                total += c
        if total < best:
            best = total
    return 0 if best == INF else best


def edge_disjoint_paths(n: int, edges: List[Tuple[int, int, int]], s: int, t: int) -> int:
    """無向圖邊不相交路徑數：每條邊拆成兩條容量 1 的有向邊後求最大流。"""
    if n <= 0 or not (0 <= s < n and 0 <= t < n) or s == t:
        return 0
    net = Dinic(n)
    for (u, v, _) in edges:
        net.add_edge(u, v, 1)
        net.add_edge(v, u, 1)
    return net.max_flow(s, t)


# ---------------------------------------------------------------- 二分圖匹配


def bipartite_matching(adj: List[List[int]], nl: int, nr: int) -> Tuple[int, List[Tuple[int, int]]]:
    """匈牙利 / Kuhn 算法。adj[u] 為左點 u 可連的右點列表。

    返回 (最大匹配數, 匹配邊 [(左點, 右點)]，按左點升冪)。
    增廣時對每個左點 DFS 嘗試「讓位」：若右點 v 已被佔，就遞迴替它現在的對象另找一個。
    """
    if nl <= 0 or nr <= 0:
        return 0, []
    match_r = [-1] * nr

    def try_kuhn(u: int, seen: List[bool]) -> bool:
        for v in adj[u]:
            if seen[v]:
                continue
            seen[v] = True
            if match_r[v] == -1 or try_kuhn(match_r[v], seen):
                match_r[v] = u
                return True
        return False

    count = 0
    for u in range(nl):
        if try_kuhn(u, [False] * nr):
            count += 1
    pairs = sorted((match_r[v], v) for v in range(nr) if match_r[v] != -1)
    return count, pairs


def matching_via_dinic(adj: List[List[int]], nl: int, nr: int) -> int:
    """把二分圖建模成流網絡求最大匹配：S→左(1)、左→右(1)、右→T(1)。"""
    if nl <= 0 or nr <= 0:
        return 0
    total = nl + nr
    src = total
    dst = total + 1
    net = Dinic(total + 2)
    for u in range(nl):
        net.add_edge(src, u, 1)
    for v in range(nr):
        net.add_edge(nl + v, dst, 1)
    for u in range(nl):
        for v in adj[u]:
            net.add_edge(u, nl + v, 1)
    return net.max_flow(src, dst)


def matching_bruteforce(adj: List[List[int]], nl: int, nr: int) -> int:
    """暴力枚舉匹配（nl ≤ 8），作為匈牙利算法的對拍。"""
    best = 0

    def rec(u: int, used: int) -> None:
        nonlocal best
        if u == nl:
            best = max(best, bin(used).count("1"))
            return
        rec(u + 1, used)                            # 左點 u 不配
        for v in adj[u]:                            # 左點 u 配到右點 v
            if not (used >> v & 1):
                rec(u + 1, used | (1 << v))

    rec(0, 0)
    return best


def normalize_adj(adj: List[List[int]], nl: int, nr: int) -> List[List[int]]:
    """過濾越界與重複的邊，讓兩個語言的行為一致。"""
    out: List[List[int]] = []
    for u in range(nl):
        row = adj[u] if u < len(adj) else []
        seen: Set[int] = set()
        clean: List[int] = []
        for v in row:
            if 0 <= v < nr and v not in seen:
                seen.add(v)
                clean.append(v)
        out.append(clean)
    return out


# ---------------------------------------------------------------- IO 與測試


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> int:
        nonlocal pos
        v = parse_int(toks[pos]) if pos < len(toks) else None
        pos += 1
        if v is None:
            return 0
        if v > CAP_MAX:
            return CAP_MAX
        if v < -CAP_MAX:
            return -CAP_MAX
        return v

    n = nxt()
    m = nxt()
    s = nxt()
    t = nxt()
    if n < 0:
        n = 0
    if n > MAX_N:
        n = MAX_N
    if m < 0:
        m = 0
    if m > MAX_M:
        m = MAX_M
    edges: List[Tuple[int, int, int]] = []
    for _ in range(m):
        u = nxt()
        v = nxt()
        c = nxt()
        if c < 0:
            c = 0
        edges.append((u, v, c))

    net = build_flow(n, edges)
    flow = net.max_flow(s, t)
    cut_cap, side, cuts = net.min_cut(s)

    nl = nxt()
    nr = nxt()
    k = nxt()
    if nl < 0:
        nl = 0
    if nr < 0:
        nr = 0
    if nl > MAX_N:
        nl = MAX_N
    if nr > MAX_N:
        nr = MAX_N
    if k < 0:
        k = 0
    if k > MAX_M:
        k = MAX_M
    adj: List[List[int]] = [[] for _ in range(nl)]
    for _ in range(k):
        u = nxt()
        v = nxt()
        if 0 <= u < nl and 0 <= v < nr:
            adj[u].append(v)
    adj = normalize_adj(adj, nl, nr)
    cnt, pairs = bipartite_matching(adj, nl, nr)

    out: List[str] = [
        str(flow),
        str(cut_cap),
        str(len(side)),
        " ".join(str(x) for x in side),
        str(len(cuts)),
    ]
    for (u, v, c) in cuts:
        out.append("%d %d %d" % (u, v, c))
    out.append(str(cnt))
    for (u, v) in pairs:
        out.append("%d %d" % (u, v))
    sys.stdout.write("\n".join(out) + "\n")


def run_tests() -> None:
    # ---- 經典範例：4 點圖，最大流 = 2（對應最小割 {0,1} / {2,3}）----
    e0 = [(0, 1, 3), (0, 2, 2), (1, 2, 1), (1, 3, 2), (2, 3, 4)]
    net0 = build_flow(4, e0)
    assert net0.max_flow(0, 3) == 5
    assert net0.max_flow(0, 3) == 0          # 殘量網絡已滿，再跑一次為 0
    net0b = build_flow(4, e0)
    assert net0b.max_flow(0, 3) == 5
    cap0, side0, cuts0 = net0b.min_cut(0)
    assert cap0 == 5
    assert side0 == [0]                      # 只有源點可達（其餘邊全滿流）
    assert cuts0 == [(0, 1, 3), (0, 2, 2)]

    # ---- 二分圖：左 3 右 3，最大匹配 = 3（完美匹配）----
    adj1 = [[0, 1], [0, 2], [1]]
    assert bipartite_matching(adj1, 3, 3) == (3, [(0, 0), (1, 2), (2, 1)])
    assert matching_via_dinic(adj1, 3, 3) == 3
    assert matching_bruteforce(adj1, 3, 3) == 3

    # ---- 非完美匹配：左 3 右 2，最多 2 ----
    adj2 = [[0, 1], [0], [1]]
    cnt2, pairs2 = bipartite_matching(adj2, 3, 2)
    assert cnt2 == 2
    assert matching_via_dinic(adj2, 3, 2) == 2
    assert matching_bruteforce(adj2, 3, 2) == 2
    assert len(pairs2) == 2 and len(set(p[0] for p in pairs2)) == 2

    # ---- 空圖 / 退化情形 ----
    assert build_flow(0, []).max_flow(0, 0) == 0
    assert bipartite_matching([], 0, 0) == (0, [])
    assert bipartite_matching([[]], 1, 0) == (0, [])
    assert matching_via_dinic([], 0, 0) == 0
    assert edmonds_karp(0, [], 0, 0) == 0
    assert min_cut_bruteforce(0, [], 0, 0) == 0
    assert edge_disjoint_paths(3, [(0, 1, 1)], 0, 2) == 0

    # ---- 邊不相交路徑：兩條 0→3 的路 ----
    eu = [(0, 1, 1), (0, 2, 1), (1, 3, 1), (2, 3, 1)]
    assert edge_disjoint_paths(4, eu, 0, 3) == 2

    # ---- 隨機對拍：Dinic vs Edmonds-Karp vs 暴力最小割 ----
    rnd = random.Random(20261008)
    for _ in range(120):
        n = rnd.randint(2, 8)
        m = rnd.randint(0, 12)
        es: List[Tuple[int, int, int]] = []
        for _ in range(m):
            u = rnd.randrange(n)
            v = rnd.randrange(n)
            c = rnd.randint(0, 9)
            es.append((u, v, c))
        s = rnd.randrange(n)
        t = rnd.randrange(n)
        if s == t:
            continue
        f1 = build_flow(n, es).max_flow(s, t)
        f2 = edmonds_karp(n, es, s, t)
        f3 = min_cut_bruteforce(n, es, s, t)
        assert f1 == f2 == f3, (n, es, s, t, f1, f2, f3)
        # 最大流 ≤ 源點出邊總容量
        assert f1 <= sum(c for (u, v, c) in es if u == s)
        # 跑完最大流後，Dinic 自己求出的割容量也必須等於流量
        net = build_flow(n, es)
        fl = net.max_flow(s, t)
        cc, sd, ct = net.min_cut(s)
        assert fl == cc
        assert s in sd and t not in sd
        assert sum(c for (u, v, c) in ct) == cc

    # ---- 隨機對拍：匈牙利 vs Dinic vs 暴力匹配 ----
    for _ in range(120):
        nl = rnd.randint(1, 7)
        nr = rnd.randint(1, 7)
        k = rnd.randint(0, 12)
        adj: List[List[int]] = [[] for _ in range(nl)]
        for _ in range(k):
            u = rnd.randrange(nl)
            v = rnd.randrange(nr)
            adj[u].append(v)
        adj = normalize_adj(adj, nl, nr)
        c1, pr = bipartite_matching(adj, nl, nr)
        c2 = matching_via_dinic(adj, nl, nr)
        c3 = matching_bruteforce(adj, nl, nr)
        assert c1 == c2 == c3, (adj, nl, nr, c1, c2, c3)
        assert len(pr) == c1
        assert len(set(p[0] for p in pr)) == c1         # 同一左點只出現一次
        assert len(set(p[1] for p in pr)) == c1         # 同一右點只出現一次
        for (u, v) in pr:
            assert v in adj[u]                          # 每條匹配邊確實存在


if __name__ == "__main__":
    data = sys.stdin.read()
    if data.strip():
        run_io(data)
    else:
        run_tests()
        print("all tests passed")
