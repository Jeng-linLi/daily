"""LCA（最近公共祖先）：倍增 / 歐拉序 RMQ，並附帶樹上距離與 k 級祖先

題意：
    給定一棵 n 個節點的無根樹（以邊表給出，0-indexed），以 0 號點為根，要求在線回答：
      1. `lca(u, v)` —— u 與 v 的最近公共祖先；
      2. `dist(u, v)` —— u 到 v 的樹上距離（邊數）；
      3. `kth_ancestor(v, k)` —— v 向上走 k 步到達的祖先（超過根則為 -1）；
      4. `is_ancestor(u, v)` —— u 是否為 v 的祖先（含 u == v）。
    預處理 O(n log n)，單次查詢 O(log n)（倍增）；另外提供歐拉序 + Sparse Table 的
    O(n log n) 預處理 / O(1) 查詢版本作交叉驗證，以及一個 O(n) 爬 parent 的樸素版當基準。

思路：
    ### 倍增（binary lifting）
    `up[k][v]` = v 的第 `2^k` 級祖先（不存在則 -1）。遞推式：
        up[0][v] = parent[v]
        up[k][v] = up[k-1][ up[k-1][v] ]
    它利用了「二進制拆分」：任何跳數都能拆成若干個 2 的冪，因此
      - 把深的點往上拉 `depth[u] - depth[v]` 步：只需對 diff 的每個為 1 的位跳一次；
      - 兩點同深度後**一起**往上跳：從大到小枚舉 k，若 `up[k][u] != up[k][v]`
        就同時跳（相等代表跳上去會「跳過頭」或已經是同一個點，不能跳）。
        最後 `up[0][u]` 就是 LCA。

    ### 歐拉序 + RMQ（另一條完全不同的路線）
    對樹做 DFS，進入節點時記一次、每訪問完一個子節點返回時再記一次父節點，
    得到長度 `2n-1` 的歐拉序列。u 與 v 的 LCA，就是歐拉序中
    `[first[u], first[v]]` 區間內 **depth 最小** 的那個節點。
    區間最小值用 Sparse Table 預處理（O(n log n) 空間、O(1) 查詢）。
    這條路線不需要「對齊深度」，是 LCA 的另一種經典解法，兩者結果必然一致 → 拿來對拍。

    ### 樹上距離
    `dist(u, v) = depth[u] + depth[v] - 2 * depth[lca(u, v)]`。

    ### 祖先判定
    DFS 進出時間 `tin / tout`：u 是 v 的祖先 ⟺ `tin[u] <= tin[v] 且 tout[v] <= tout[u]`。

輸入格式（stdin，全部以空白分隔）：
    n
    n-1 行：u v        （無向邊；n == 1 時沒有邊）
    q
    q 行：lca u v  |  anc v k
輸出格式（stdout）：
    每個查詢一行：
      `lca u v` -> `<lca> <dist>`
      `anc v k` -> `<ancestor>`（不存在為 -1）
    節點編號越界、或該點不在根 0 所在的連通分量內時，輸出 `-1 -1` / `-1`。
輸入被截斷時，缺的部分按 0 個查詢處理。無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import List, Optional, Tuple

# 只接受 [+-]?digits，與 C++ 版本的 tryLL 完全一致；
# 其它 token（例如 op 名稱）一律視為「輸入到此為止」，避免兩語言一個崩一個不崩。
_INT_RE = re.compile(r"^[+-]?[0-9]+$")


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；不是十進制整數則返回 None。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


class LCA:
    """倍增 LCA + 歐拉序 RMQ，兩條獨立路線互相對拍。"""

    def __init__(self, n: int, edges: List[Tuple[int, int]], root: int = 0) -> None:
        self.n = n
        self.root = root
        self.adj: List[List[int]] = [[] for _ in range(n)]
        for u, v in edges:
            if 0 <= u < n and 0 <= v < n and u != v:
                self.adj[u].append(v)
                self.adj[v].append(u)
        for lst in self.adj:
            lst.sort()                      # 排序保證 DFS 順序確定（兩語言一致）

        # LOG = 足夠覆蓋 n 層的位數
        self.LOG = 1
        while (1 << self.LOG) <= max(n, 1):
            self.LOG += 1

        self.depth = [-1] * n
        self.tin = [-1] * n
        self.tout = [-1] * n
        self.up: List[List[int]] = [[-1] * n for _ in range(self.LOG)]
        self.euler: List[int] = []
        self.euler_depth: List[int] = []

        self._build_tree()
        self._build_lifting()
        self._build_euler()
        self._build_sparse()

    # ------------------------------------------------------------ 建樹（迭代 DFS）
    def _build_tree(self) -> None:
        if self.n == 0 or not (0 <= self.root < self.n):
            return
        stack = [(self.root, -1)]
        self.depth[self.root] = 0
        self.up[0][self.root] = -1
        while stack:
            u, p = stack.pop()
            for v in reversed(self.adj[u]):
                if v == p or self.depth[v] >= 0:
                    continue
                self.depth[v] = self.depth[u] + 1
                self.up[0][v] = u
                stack.append((v, u))

    # ------------------------------------------------------------ 倍增表
    def _build_lifting(self) -> None:
        for k in range(1, self.LOG):
            prev, cur = self.up[k - 1], self.up[k]
            for v in range(self.n):
                p = prev[v]
                cur[v] = -1 if p < 0 else prev[p]

    # ------------------------------------------------------------ 歐拉序（迭代 DFS）
    def _build_euler(self) -> None:
        if self.n == 0 or self.depth[self.root] < 0:
            return
        # 棧元素 = (節點, 父節點, 下一個要訪問的孩子下標)
        stack: List[List[int]] = [[self.root, -1, 0]]
        while stack:
            frame = stack[-1]
            u, p, ci = frame[0], frame[1], frame[2]
            if ci == 0:                     # 第一次進入 u
                self.tin[u] = len(self.euler)
                self.euler.append(u)
                self.euler_depth.append(self.depth[u])
            if ci < len(self.adj[u]):
                frame[2] = ci + 1
                v = self.adj[u][ci]
                if v != p:
                    stack.append([v, u, 0])
            else:
                stack.pop()
                self.tout[u] = len(self.euler) - 1
                if p >= 0:                  # 從子節點返回，再記一次父節點
                    self.euler.append(p)
                    self.euler_depth.append(self.depth[p])

    # ------------------------------------------------------------ Sparse Table
    def _build_sparse(self) -> None:
        m = len(self.euler_depth)
        self.st_k = 1
        while (1 << self.st_k) <= max(m, 1):
            self.st_k += 1
        self.log_tbl = [0] * (m + 1)
        for i in range(2, m + 1):
            self.log_tbl[i] = self.log_tbl[i // 2] + 1
        self.st: List[List[int]] = [list(range(m))]
        for k in range(1, self.st_k):
            prev = self.st[k - 1]
            half = 1 << (k - 1)
            nxt = []
            for i in range(0, m - (1 << k) + 1):
                a = prev[i]
                b = prev[i + half]
                # 深度小者勝；深度相同取歐拉序下標小者（保證確定）
                nxt.append(a if self.euler_depth[a] <= self.euler_depth[b] else b)
            self.st.append(nxt)

    def _rmq(self, l: int, r: int) -> int:
        """返回 [l, r] 內 depth 最小的歐拉序下標。"""
        k = self.log_tbl[r - l + 1]
        a = self.st[k][l]
        b = self.st[k][r - (1 << k) + 1]
        return a if self.euler_depth[a] <= self.euler_depth[b] else b

    # ------------------------------------------------------------ 查詢
    def _valid(self, v: int) -> bool:
        return 0 <= v < self.n and self.depth[v] >= 0

    def lca(self, u: int, v: int) -> int:
        """倍增版 LCA，O(log n)。"""
        if not self._valid(u) or not self._valid(v):
            return -1
        if self.depth[u] < self.depth[v]:
            u, v = v, u
        diff = self.depth[u] - self.depth[v]
        k = 0
        while diff:
            if diff & 1:
                u = self.up[k][u]
                if u < 0:
                    return -1
            diff >>= 1
            k += 1
        if u == v:
            return u
        for k in range(self.LOG - 1, -1, -1):
            pu, pv = self.up[k][u], self.up[k][v]
            if pu != pv:                    # 相等代表跳上去會跳過頭
                u = pu if pu >= 0 else u
                v = pv if pv >= 0 else v
        return self.up[0][u]

    def lca_rmq(self, u: int, v: int) -> int:
        """歐拉序 + Sparse Table 版 LCA，O(1) 查詢。"""
        if not self._valid(u) or not self._valid(v):
            return -1
        l, r = self.tin[u], self.tin[v]
        if l > r:
            l, r = r, l
        return self.euler[self._rmq(l, r)]

    def lca_naive(self, u: int, v: int) -> int:
        """逐層往上爬的樸素基準，O(n)。"""
        if not self._valid(u) or not self._valid(v):
            return -1
        a, b = u, v
        while self.depth[a] > self.depth[b]:
            a = self.up[0][a]
        while self.depth[b] > self.depth[a]:
            b = self.up[0][b]
        while a != b:
            a = self.up[0][a]
            b = self.up[0][b]
        return a

    def dist(self, u: int, v: int) -> int:
        """樹上距離（邊數）；不在同一棵樹時返回 -1。"""
        w = self.lca(u, v)
        if w < 0:
            return -1
        return self.depth[u] + self.depth[v] - 2 * self.depth[w]

    def kth_ancestor(self, v: int, k: int) -> int:
        """v 的第 k 級祖先（k == 0 為自己；超過根則 -1）。"""
        if not self._valid(v) or k < 0:
            return -1
        if k > self.depth[v]:
            return -1
        i = 0
        while k:
            if k & 1:
                v = self.up[i][v]
                if v < 0:
                    return -1
            k >>= 1
            i += 1
        return v

    def kth_ancestor_naive(self, v: int, k: int) -> int:
        """逐層往上爬的樸素基準。"""
        if not self._valid(v) or k < 0:
            return -1
        if k > self.depth[v]:
            return -1
        cur = v
        for _ in range(k):
            cur = self.up[0][cur]
        return cur

    def is_ancestor(self, u: int, v: int) -> bool:
        """u 是否為 v 的祖先（含 u == v）。"""
        if not self._valid(u) or not self._valid(v):
            return False
        return self.tin[u] <= self.tin[v] and self.tout[v] <= self.tout[u]


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
    edges: List[Tuple[int, int]] = []
    if n > 1:
        for _ in range(n - 1):
            a, b = parse_int(nxt()), parse_int(nxt())
            if a is None or b is None:
                break                        # 輸入截斷或出現非數字 token
            edges.append((a, b))
    q = next_int(0)

    solver = LCA(n, edges, 0)
    for _ in range(q):
        op = nxt()
        if op is None:
            break
        if op == "lca":
            u, v = parse_int(nxt()), parse_int(nxt())
            if u is None or v is None:
                break
            print(solver.lca(u, v), solver.dist(u, v))
        elif op == "anc":
            v, k = parse_int(nxt()), parse_int(nxt())
            if v is None or k is None:
                break
            print(solver.kth_ancestor(v, k))


# ---------------------------------------------------------------- 測試
def build_random_tree(n: int, rng: random.Random) -> List[Tuple[int, int]]:
    """隨機生成一棵 n 個點的樹（用隨機 parent 的方式保證連通）。"""
    edges = []
    for v in range(1, n):
        p = rng.randint(0, v - 1)
        edges.append((p, v))
    return edges


def run_tests() -> None:
    # 單點樹
    s = LCA(1, [], 0)
    assert s.depth == [0]
    assert s.lca(0, 0) == 0
    assert s.dist(0, 0) == 0
    assert s.kth_ancestor(0, 0) == 0
    assert s.kth_ancestor(0, 1) == -1
    assert s.is_ancestor(0, 0) is True
    assert s.lca_rmq(0, 0) == 0
    assert s.lca_naive(0, 0) == 0

    # README 示例樹：0-1, 0-2, 1-3, 1-4, 2-5, 2-6
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
    s = LCA(7, edges, 0)
    assert s.depth == [0, 1, 1, 2, 2, 2, 2]
    assert s.lca(3, 4) == 1
    assert s.dist(3, 4) == 2
    assert s.lca(3, 5) == 0
    assert s.dist(3, 5) == 4
    assert s.lca(5, 6) == 2
    assert s.dist(5, 6) == 2
    assert s.lca(0, 3) == 0
    assert s.dist(0, 3) == 2
    assert s.lca(3, 3) == 3
    assert s.dist(3, 3) == 0
    assert s.kth_ancestor(3, 0) == 3
    assert s.kth_ancestor(3, 1) == 1
    assert s.kth_ancestor(3, 2) == 0
    assert s.kth_ancestor(3, 3) == -1
    assert s.is_ancestor(0, 6) is True
    assert s.is_ancestor(1, 6) is False
    assert s.is_ancestor(2, 2) is True
    # 三種路線一致
    for u in range(7):
        for v in range(7):
            assert s.lca(u, v) == s.lca_rmq(u, v) == s.lca_naive(u, v)

    # 鏈狀樹（深度最大）：0-1-2-...-9
    chain = [(i, i + 1) for i in range(9)]
    s = LCA(10, chain, 0)
    assert s.depth == list(range(10))
    assert s.lca(0, 9) == 0
    assert s.dist(0, 9) == 9
    assert s.dist(3, 7) == 4
    assert s.kth_ancestor(9, 9) == 0
    assert s.kth_ancestor(9, 10) == -1
    for u in range(10):
        for v in range(10):
            assert s.lca(u, v) == min(u, v)
            assert s.dist(u, v) == abs(u - v)
            assert s.lca(u, v) == s.lca_rmq(u, v) == s.lca_naive(u, v)

    # 星狀樹：中心 0，其餘都是葉子
    star = [(0, i) for i in range(1, 8)]
    s = LCA(8, star, 0)
    for u in range(1, 8):
        for v in range(1, 8):
            assert s.lca(u, v) == (0 if u != v else u)
            assert s.dist(u, v) == (2 if u != v else 0)

    # 森林 / 不連通：節點 4, 5 不在根的分量內
    forest_edges = [(0, 1), (1, 2), (4, 5)]
    s = LCA(6, forest_edges, 0)
    assert s.depth[3] == -1
    assert s.depth[4] == -1
    assert s.lca(0, 4) == -1
    assert s.dist(0, 4) == -1
    assert s.lca(0, 2) == 0
    assert s.dist(0, 2) == 2
    assert s.kth_ancestor(4, 0) == -1
    assert s.is_ancestor(0, 4) is False

    # 自環與越界邊被忽略
    s = LCA(3, [(0, 0), (0, 1), (1, 2), (1, 99)], 0)
    assert s.depth == [0, 1, 2]
    assert s.lca(1, 2) == 1

    # 空樹
    s = LCA(0, [], 0)
    assert s.lca(0, 0) == -1
    assert s.dist(0, 0) == -1

    # 越界查詢
    s = LCA(3, [(0, 1), (1, 2)], 0)
    assert s.lca(0, 7) == -1
    assert s.lca(-1, 1) == -1
    assert s.dist(0, 9) == -1
    assert s.kth_ancestor(9, 1) == -1
    assert s.kth_ancestor(1, -1) == -1
    assert s.is_ancestor(0, 9) is False

    # 隨機對拍：倍增 vs RMQ vs 樸素；k 級祖先 vs 樸素；祖先判定 vs 暴力上爬
    rng = random.Random(20261005)
    for _ in range(120):
        n = rng.randint(1, 14)
        edges = build_random_tree(n, rng)
        s = LCA(n, edges, 0)
        assert s.depth[0] == 0
        assert all(d >= 0 for d in s.depth)          # 隨機 parent 建出來必連通
        for u in range(n):
            for v in range(n):
                w = s.lca(u, v)
                assert w == s.lca_rmq(u, v) == s.lca_naive(u, v)
                assert s.dist(u, v) == s.depth[u] + s.depth[v] - 2 * s.depth[w]
            for k in range(0, n + 2):
                assert s.kth_ancestor(u, k) == s.kth_ancestor_naive(u, k)
            # 祖先判定對拍：暴力沿 parent 上爬
            for v in range(n):
                brute = False
                cur = v
                while cur >= 0:
                    if cur == u:
                        brute = True
                        break
                    cur = s.up[0][cur]
                assert s.is_ancestor(u, v) == brute

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
