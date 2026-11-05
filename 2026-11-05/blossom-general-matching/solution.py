"""帶花樹 Blossom：一般圖最大匹配（Edmonds 花收縮）

題意：
    給定 n 個點、m 條邊的無向簡單圖（可能有奇環，不限二分圖），
    求最大匹配（Maximum Matching）：最多的兩兩不共點的邊集合。
    輸出匹配數，以及每個點的配對點（未匹配點輸出 0）。

    輸入：第一行 n m；接下來 m 行每行兩個整數 u v（1-indexed）。
    輸出：第一行一個整數 k；第二行 n 個整數，第 i 個是點 i 的配對點
          （0 表示未匹配）。

思路：
    ### 1. 二分圖匹配為什麼不夠
    匈牙利算法靠「交錯路徑」（交替出現非匹配邊 / 匹配邊）增廣。
    一般圖裡交錯路徑可能經過**奇環**：奇環上的一點既可以按「內點」
    也可以按「外點」的規則延伸，二分圖的交替染色在奇環上失效。

    ### 2. 花與收縮（Edmonds 1965）
    把「花（blossom）」定義為：一條從樹根出發的交錯路徑 + 一條
    閉合奇環的交錯邊。花的關鍵性質：**整朵花可以像一個點一樣對待**
    ——從花上任意一點都能以正確的奇偶性延伸增廣路。
    找到花時，把花上所有點收縮成一個新節點（並查集維護 base），
    繼續 BFS；找到增廣路後再「展開」花，沿 parent 回溯翻轉匹配邊。

    ### 3. 實作要點
    - base[i]：i 所在花的代表點（並查集，帶路徑壓縮）。
    - p[i]：交錯森林中 i 的父親（i 是「外點」時有意義）。
    - 每輪從某個未匹配點開始 BFS：外點把所有鄰居分類——
      未訪問 → 若鄰居未匹配則找到增廣路，否則鄰居的配偶入隊（內點）；
      遇到另一個外點且 base 不同 → 找到花，lca 後整朵收縮。
    - 找到增廣路後，倒序翻轉：match[u] = v, match[v] = u 沿 parent 鏈。

    ### 4. 對拍
    內建測試用「枚舉邊子集」的暴力（n ≤ 10, m ≤ 15）：
    子集內每個點度數 <= 1 即為合法匹配，取最大值比對；
    另含三角形（奇環）、雙三角形、Petersen 圖等固定用例。

複雜度：
    每個未匹配點最多增廣一次；單次找增廣路 O(n^2)～O(n^3)（含花的 lca），
    總時間 O(n^3)，空間 O(n + m)。

應用場景（Applications）：
    - **一般配對問題**：舞伴配對、宿舍分配、導師-學生互選（圖不是二分圖）
    - **網路設計**：衛星鏈路冗餘規劃中的邊不相交配對
    - 二分圖場景用匈牙利 / Dinic 建模即可（見本倉庫 2026-10-08），
      Blossom 是「圖中存在奇環」時的通用解
"""

import sys
from collections import deque
from typing import Dict, List, Optional, Sequence, Tuple

NMAX = 300
VAL_LIM = 10 ** 18


class Blossom:
    """一般圖最大匹配（帶花樹）。點編號 0..n-1。"""

    def __init__(self, n: int, adj: Sequence[Sequence[int]]) -> None:
        self.n = n
        self.adj: List[List[int]] = [list(a) for a in adj]
        self.match = [-1] * n      # match[u] = 與 u 配對的點（-1 未匹配）
        self.p = [0] * n           # 交錯森林父親（外點）
        self.base = [0] * n        # 所在花的代表（並查集根）
        self.used = [False] * n

    # ---- 並查集（帶路徑壓縮）----
    def _find_base(self, x: int) -> int:
        while self.base[x] != x:
            self.base[x] = self.base[self.base[x]]
            x = self.base[x]
        return x

    # ---- 求兩個外點在交錯森林中的最近公共花代表 ----
    def _lca(self, a: int, b: int) -> int:
        mark = [False] * self.n
        while True:
            a = self._find_base(a)
            mark[a] = True
            if self.match[a] == -1:
                break
            a = self.p[self.match[a]]
        while True:
            b = self._find_base(b)
            if mark[b]:
                return b
            b = self.p[self.match[b]]

    # ---- 把 a 到花代表 b 的路徑全標成花，child 是閉合邊的另一端 ----
    def _mark_path(self, v: int, b: int, child: int, blossom: List[bool]) -> None:
        while self._find_base(v) != b:
            blossom[self._find_base(v)] = True
            blossom[self._find_base(self.match[v])] = True
            self.p[v] = child
            child = self.match[v]
            v = self.p[self.match[v]]

    # ---- 從 root（未匹配點）BFS 找一條增廣路；找到並翻轉後回 True ----
    def _find_path(self, root: int) -> bool:
        n = self.n
        self.used = [False] * n
        self.p = [-1] * n
        for i in range(n):
            self.base[i] = i
        self.used[root] = True
        q: deque = deque([root])
        while q:
            v = q.popleft()
            for to in self.adj[v]:
                # 同一朵花內、或 v-to 本身是匹配邊：跳過
                if self._find_base(v) == self._find_base(to) or self.match[v] == to:
                    continue
                # 情況 1：to 是已被拜訪過的外點（或樹根）→ 找到花
                if to == root or (self.match[to] != -1 and self.p[self.match[to]] != -1):
                    curbase = self._lca(v, to)
                    blossom = [False] * n
                    self._mark_path(v, curbase, to, blossom)
                    self._mark_path(to, curbase, v, blossom)
                    for i in range(n):
                        if blossom[self._find_base(i)]:
                            self.base[i] = curbase
                            if not self.used[i]:
                                self.used[i] = True
                                q.append(i)
                # 情況 2：to 還沒在交錯森林裡
                elif self.p[to] == -1:
                    self.p[to] = v
                    if self.match[to] == -1:
                        # 找到增廣路：沿 p 鏈倒序翻轉匹配邊
                        u = to
                        while u != -1:
                            pv = self.p[u]
                            ppv = self.match[pv]
                            self.match[u] = pv
                            self.match[pv] = u
                            u = ppv
                        return True
                    else:
                        # to 已匹配 → 其配偶作為內點延伸
                        self.used[self.match[to]] = True
                        q.append(self.match[to])
        return False

    def solve(self) -> Tuple[int, List[int]]:
        res = 0
        for v in range(self.n):
            if self.match[v] == -1:
                if self._find_path(v):
                    res += 1
        return res, self.match


def general_matching(n: int, edges: Sequence[Tuple[int, int]]) -> Tuple[int, List[int]]:
    """回傳 (匹配數, match 陣列；未匹配為 -1 → 對外轉 0）。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    b = Blossom(n, adj)
    res, match = b.solve()
    return res, [m + 1 if m != -1 else 0 for m in match]


# --------------------------------------------------------------- 暴力對拍

def brute_max_matching(n: int, edges: Sequence[Tuple[int, int]]) -> int:
    """枚舉邊子集（m <= 15），取合法匹配的最大值。"""
    m = len(edges)
    best = 0
    for mask in range(1 << m):
        deg = [0] * n
        cnt = 0
        ok = True
        mm = mask
        i = 0
        while mm:
            if mm & 1:
                a, b = edges[i]
                deg[a] += 1
                deg[b] += 1
                if deg[a] > 1 or deg[b] > 1:
                    ok = False
                    break
                cnt += 1
            mm >>= 1
            i += 1
        if ok and cnt > best:
            best = cnt
    return best


# --------------------------------------------------------------- 內建測試

def run_tests() -> None:
    # 固定用例 1：三角形（奇環）→ 最大匹配 1
    res, match = general_matching(3, [(0, 1), (1, 2), (2, 0)])
    assert res == 1
    nz = [x for x in match if x]
    assert len(nz) == 2 and match[nz[0] - 1] == nz[1]

    # 固定用例 2：兩個三角形由一條邊相連 → 匹配數 3
    edges2 = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)]
    res, _ = general_matching(6, edges2)
    assert res == 3, f"got {res}"

    # 固定用例 3：Petersen 圖（10 點 15 邊）→ 有完美匹配，答案 5
    pet = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
           (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
           (5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]
    res, _ = general_matching(10, pet)
    assert res == 5, f"got {res}"

    # 固定用例 4：無邊
    res, match = general_matching(4, [])
    assert res == 0 and match == [0, 0, 0, 0]

    # 固定用例 5：奇環 + 尾巴（經典 blossom 觸發形態）
    edges5 = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)]
    res, _ = general_matching(5, edges5)
    assert res == 2, f"got {res}"

    # 隨機對拍（n <= 10, m <= 15）
    rng_state = 987654321

    def rnd(mod: int) -> int:
        nonlocal rng_state
        rng_state = (rng_state * 1103515245 + 12345) % (1 << 31)
        return rng_state % mod

    for it in range(300):
        n = rnd(10) + 1
        # 枚举所有点对后 Fisher-Yates 洗牌取前 m 条：
        # 直接 while 循环抽点在 LCG 低比特周期短时会凑不满 m 条而死循环
        all_edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
        for i in range(len(all_edges) - 1, 0, -1):
            j = rnd(i + 1)
            all_edges[i], all_edges[j] = all_edges[j], all_edges[i]
        m = rnd(min(len(all_edges), 15) + 1)
        edges = all_edges[:m]
        got, match = general_matching(n, edges)
        want = brute_max_matching(n, edges)
        assert got == want, f"random#{it} n={n} edges={edges} got={got} want={want}"
        # 自洽：match 陣列必須互相配對且無共用點
        for a, b in enumerate(match):
            if b != 0:
                assert match[b - 1] == a + 1, f"random#{it} match broken"


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

    n = _parse_int(nxt() or "")
    if n is None or not (0 <= n <= NMAX):
        bad()
        return
    m = _parse_int(nxt() or "")
    if m is None or not (0 <= m <= n * NMAX):
        bad()
        return
    edges: List[Tuple[int, int]] = []
    for _ in range(m):
        a = _parse_int(nxt() or "")
        b = _parse_int(nxt() or "")
        if a is None or b is None or not (1 <= a <= n and 1 <= b <= n) or a == b:
            bad()
            return
        edges.append((a - 1, b - 1))
    res, match = general_matching(n, edges)
    emit(str(res) + "\n")
    emit(" ".join(str(x) for x in match) + "\n")


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
