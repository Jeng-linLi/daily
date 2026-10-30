"""虛擬樹（Virtual Tree）：把 O(n) 的原樹壓縮成 O(k) 的小樹

題意：
    給定一棵有 n 個節點的根樹（根為 1），節點 u 有權重 w[u]。
    有 m 次詢問，每次詢問給出 k 個「關鍵節點」（key nodes）。
    對每次詢問，要建立這 k 個關鍵節點的**虛擬樹**並回答：
        (1) 虛擬樹共有幾個節點
        (2) 虛擬樹所有邊的長度總和（邊長 = 原樹上兩端點的距離）
        (3) 對每個關鍵節點 u（依節點編號遞增輸出）：
            cnt = 位於 u 的虛擬樹子樹中的關鍵節點個數
            sum = 位於 u 的虛擬樹子樹中的關鍵節點的權重總和

    限制：1 <= n <= 100000，0 <= m <= 100000，0 <= k <= n，
          |w[u]| <= 10^9。

思路：
    ### 1. 為什麼需要虛擬樹
    原樹有 n 個節點，但每次詢問只關心 k 個關鍵節點（k << n）。
    若每次都在原樹上做 DFS，複雜度是 O(nm)，完全不可行。
    虛擬樹的想法是：**只保留 k 個關鍵節點，以及它們彼此的 LCA**，
    把中間沒有關鍵節點的部分壓縮成一條邊。壓縮後的樹最多 2k-1 個節點。

    ### 2. 虛擬樹的節點集合
    令 V = {關鍵節點} ∪ {LCA(相鄰的兩個關鍵節點)}。
    這裡「相鄰」指的是把關鍵節點依 DFS 序排序後相鄰的那一對。
    定理：這 k-1 個 LCA 恰好涵蓋了所有需要的內部節點 ——
    因為一個節點若是兩個關鍵節點的 LCA，那麼在 DFS 序上它們之間
    必然存在一對「相鄰且分居兩側」的關鍵節點。
    所以 |V| <= 2k-1。

    ### 3. 建構虛擬樹
    把 V 依 DFS 序排序，用一個堆疊線性建樹：
    對每個節點 v，把堆疊頂端「不是 v 的祖先」的節點彈掉，
    剩下的堆疊頂就是 v 在虛擬樹上的父親。
    因為 V 已依 DFS 序排序，任何時刻堆疊中的節點都是一條祖先鏈。
    邊長 = depth[child] - depth[parent]。

    ### 4. 樹上統計
    虛擬樹建好後，cnt / sum 就是普通的樹上 DP：
    先把關鍵節點的自身權重放進去，再由子往父累加。
    因為虛擬樹的節點依 DFS 序排序時「父一定排在子前面」，
    倒序掃描一次即可完成累加。

    ### 5. 為什麼 cnt / sum 跟原樹子樹的答案一樣
    虛擬樹保留了關鍵節點之間的祖先關係：若 u 是 v 的祖先
    （u、v 都是關鍵節點），則 u 在虛擬樹上依然是 v 的祖先。
    所以「u 的虛擬樹子樹中的關鍵節點」與「u 的原樹子樹中的關鍵節點」
    是同一個集合。

複雜度：
    預處理（LCA 倍增表 + DFS 序）  時間 O(n log n)，空間 O(n log n)
    單次詢問                      O(k log k + k log n)
    總計                          O((n + sum(k)) log n)

應用場景：
    - 樹上批量路徑查詢（Virtual Tree / 虛擬樹優化）：把 O(n) 的一次遍歷
      壓縮成 O(k)，是「多次少量查詢」的標準加速手段
    - 靜態 LCA 問題（ZJOI / SNOI 的「接近的祖先」類題目）
    - 集合運算型子樹查詢：合併子集合、統計子樹內不同元素的個數
    - 動態規劃換根 / 樹形 DP 中的「壓縮樹」技巧
    - 莫隊 + 樹上啟發式合併的基礎設施
"""

import sys
from typing import List, Optional, Tuple

NMAX = 100000
MMAX = 100000
WMAX = 10 ** 9


class Tree:
    """根樹 + DFS 序 + 倍增 LCA。"""

    def __init__(self, n: int, parent: List[int], weight: List[int]) -> None:
        self.n = n
        self.parent = parent
        self.weight = weight
        self.children: List[List[int]] = [[] for _ in range(n + 1)]
        for u in range(2, n + 1):
            self.children[parent[u]].append(u)
        self.depth: List[int] = [0] * (n + 1)
        self.tin: List[int] = [0] * (n + 1)
        self.tout: List[int] = [0] * (n + 1)
        self.order: List[int] = []
        self.up: List[List[int]] = []
        self._build()

    def _build(self) -> None:
        n = self.n
        # 迭代 DFS：避免深鏈造成 Python 遞迴爆堆
        stack = [1]
        order = self.order
        while stack:
            u = stack.pop()
            self.tin[u] = len(order)
            order.append(u)
            for c in self.children[u]:
                self.depth[c] = self.depth[u] + 1
                stack.append(c)
        # 子樹大小 -> tout（DFS 序中子樹是一段連續區間）
        size = [1] * (n + 1)
        for u in reversed(order):
            for c in self.children[u]:
                size[u] += size[c]
            self.tout[u] = self.tin[u] + size[u]
        # 倍增表
        log = max(1, n.bit_length())
        up = [self.parent]
        for _ in range(1, log):
            prev = up[-1]
            up.append([prev[prev[u]] for u in range(n + 1)])
        self.up = up

    def is_ancestor(self, u: int, v: int) -> bool:
        return self.tin[u] <= self.tin[v] < self.tout[u]

    def lca(self, a: int, b: int) -> int:
        if self.depth[a] < self.depth[b]:
            a, b = b, a
        diff = self.depth[a] - self.depth[b]
        k = 0
        while diff:
            if diff & 1:
                a = self.up[k][a]
            diff >>= 1
            k += 1
        if a == b:
            return a
        for k in range(len(self.up) - 1, -1, -1):
            if self.up[k][a] != self.up[k][b]:
                a = self.up[k][a]
                b = self.up[k][b]
        return self.parent[a]


def build_virtual(tree: Tree, keys: List[int]) -> Tuple[List[int], List[int]]:
    """建立虛擬樹。

    回傳 (vnodes, vparent)：
        vnodes[i]  —— 虛擬樹第 i 個節點的原編號（依 DFS 序遞增）
        vparent[i] —— 第 i 個節點的虛擬父節點在 vnodes 中的索引，根為 -1
    """
    if not keys:
        return [], []
    skeys = sorted(keys, key=lambda x: tree.tin[x])
    nodes = list(skeys)
    for i in range(len(skeys) - 1):
        nodes.append(tree.lca(skeys[i], skeys[i + 1]))
    uniq = sorted(set(nodes), key=lambda x: tree.tin[x])
    vparent = [-1] * len(uniq)
    stack: List[int] = []
    for idx, v in enumerate(uniq):
        while stack and not tree.is_ancestor(uniq[stack[-1]], v):
            stack.pop()
        if stack:
            vparent[idx] = stack[-1]
        stack.append(idx)
    return uniq, vparent


def solve_query(tree: Tree, keys: List[int]) -> Tuple[int, int, List[Tuple[int, int]]]:
    """回傳 (虛擬樹節點數, 邊長總和, [(cnt, sum), ...])。

    第三個列表依關鍵節點編號遞增排列。
    """
    uniq, vparent = build_virtual(tree, keys)
    k = len(uniq)
    if k == 0:
        return 0, 0, []
    keyset = set(keys)
    cnt = [0] * k
    wsum = [0] * k
    for i, v in enumerate(uniq):
        if v in keyset:
            cnt[i] = 1
            wsum[i] = tree.weight[v]
    edge_total = 0
    # 倒序掃描：子節點的 DFS 序必定大於父節點，故子先於父被處理
    for i in range(k - 1, -1, -1):
        p = vparent[i]
        if p >= 0:
            cnt[p] += cnt[i]
            wsum[p] += wsum[i]
            edge_total += tree.depth[uniq[i]] - tree.depth[uniq[p]]
    pos = {v: i for i, v in enumerate(uniq)}
    res = [(cnt[pos[u]], wsum[pos[u]]) for u in sorted(keys)]
    return k, edge_total, res


def brute_query(tree: Tree, keys: List[int]) -> Tuple[int, int, List[Tuple[int, int]]]:
    """暴力參考實作：O(k^2) 求虛擬樹節點集合與每個關鍵節點的統計。"""
    if not keys:
        return 0, 0, []
    vset = set(keys)
    for a in keys:
        for b in keys:
            vset.add(tree.lca(a, b))
    edge_total = 0
    for v in vset:
        best_anc = None
        for u in vset:
            if u != v and tree.is_ancestor(u, v):
                if best_anc is None or tree.depth[u] > tree.depth[best_anc]:
                    best_anc = u
        if best_anc is not None:
            edge_total += tree.depth[v] - tree.depth[best_anc]
    res: List[Tuple[int, int]] = []
    for u in sorted(keys):
        c = 0
        s = 0
        for v in keys:
            if tree.is_ancestor(u, v):
                c += 1
                s += tree.weight[v]
        res.append((c, s))
    return len(vset), edge_total, res


def run_tests() -> None:
    # ---- 固定用例 ----
    # 樹：1-2, 1-3, 2-4, 2-5, 3-6（根 = 1）
    parent = [0, 0, 1, 1, 2, 2, 3]
    weight = [0, 10, 20, 30, 40, 50, 60]
    tree = Tree(6, parent, weight)

    # 單一關鍵節點：虛擬樹只有它自己
    n, e, res = solve_query(tree, [4])
    assert (n, e, res) == (1, 0, [(1, 40)]), (n, e, res)
    assert brute_query(tree, [4]) == (1, 0, [(1, 40)])

    # 兩個兄弟節點 4、5：LCA = 2，虛擬樹是 2-4、2-5，邊長 1+1 = 2
    #   注意 cnt 是「自己子樹內的關鍵節點數」，4 的子樹只有 4 自己
    n, e, res = solve_query(tree, [4, 5])
    assert (n, e) == (3, 2), (n, e)
    assert res == [(1, 40), (1, 50)], res
    assert brute_query(tree, [4, 5]) == (3, 2, [(1, 40), (1, 50)])

    # 祖先與後代：3 與 6，LCA = 3，虛擬樹只有 3-6 一條邊
    #   3 的虛擬樹子樹含 {3, 6}，故 cnt = 2、sum = 30+60 = 90
    n, e, res = solve_query(tree, [3, 6])
    assert (n, e) == (2, 1), (n, e)
    assert res == [(2, 90), (1, 60)], res
    assert brute_query(tree, [3, 6]) == (2, 1, [(2, 90), (1, 60)])

    # 跨分支：4 與 6，LCA = 1，虛擬樹節點 {1, 4, 6}，
    #   邊 1-4 長 2、1-6 長 2，總和 4
    n, e, res = solve_query(tree, [4, 6])
    assert (n, e) == (3, 4), (n, e)
    assert res == [(1, 40), (1, 60)], res
    assert brute_query(tree, [4, 6]) == (3, 4, [(1, 40), (1, 60)])

    # 四個關鍵節點：V = {1, 2, 4, 5, 6}，共 5 個節點
    #   1 的子樹含全部 4 個關鍵節點 → cnt = 4、sum = 40+50+60+10 = 160
    n, e, res = solve_query(tree, [1, 4, 5, 6])
    assert (n, e) == (5, 5), (n, e)
    assert res == [(4, 160), (1, 40), (1, 50), (1, 60)], res
    assert brute_query(tree, [1, 4, 5, 6]) == (5, 5, [(4, 160), (1, 40), (1, 50), (1, 60)])

    # 空詢問
    assert solve_query(tree, []) == (0, 0, [])
    assert brute_query(tree, []) == (0, 0, [])

    # 負權重：鏈 1-2-3，權重 -5 / 7 / -1
    parent2 = [0, 0, 1, 2]
    weight2 = [0, -5, 7, -1]
    tree2 = Tree(3, parent2, weight2)
    n, e, res = solve_query(tree2, [1, 3])
    assert (n, e) == (2, 2), (n, e)
    assert res == [(2, -6), (1, -1)], res
    assert brute_query(tree2, [1, 3]) == (2, 2, [(2, -6), (1, -1)])

    # 單節點樹
    tree3 = Tree(1, [0, 0], [0, 42])
    assert solve_query(tree3, [1]) == (1, 0, [(1, 42)])
    assert brute_query(tree3, [1]) == (1, 0, [(1, 42)])

    # ---- 隨機對拍 ----
    import random

    rng = random.Random(20261030)
    for _ in range(200):
        n = rng.randint(1, 14)
        par = [0] * (n + 1)
        for u in range(2, n + 1):
            par[u] = rng.randint(1, u - 1)
        w = [0] + [rng.randint(-20, 20) for _ in range(n)]
        t = Tree(n, par, w)
        k = rng.randint(0, n)
        keys = rng.sample(range(1, n + 1), k)
        assert solve_query(t, keys) == brute_query(t, keys), (n, par, keys)

    # 較大規模：隨機樹 + 隨機關鍵節點集合
    rng = random.Random(99887766)
    for _ in range(20):
        n = rng.randint(50, 400)
        par = [0] * (n + 1)
        for u in range(2, n + 1):
            par[u] = rng.randint(1, u - 1)
        w = [0] + [rng.randint(-1000, 1000) for _ in range(n)]
        t = Tree(n, par, w)
        keys = sorted(rng.sample(range(1, n + 1), rng.randint(1, n)))
        assert solve_query(t, keys) == brute_query(t, keys)

    # 深鏈：驗證迭代 DFS 不爆堆
    n = 5000
    par = [0] * (n + 1)
    for u in range(2, n + 1):
        par[u] = u - 1
    t = Tree(n, par, [0] + [1] * n)
    n_, e_, res_ = solve_query(t, [1, n // 2, n])
    assert n_ == 3 and e_ == n - 1, (n_, e_)
    assert res_[0] == (3, 3) and res_[2] == (1, 1), res_


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

    n = _parse_int(nxt() or "")
    if n is None or not (1 <= n <= NMAX):
        bad()
        return
    weight = [0] * (n + 1)
    for u in range(1, n + 1):
        v = _parse_int(nxt() or "")
        if v is None or abs(v) > WMAX:
            bad()
            return
        weight[u] = v
    parent = [0] * (n + 1)
    for u in range(2, n + 1):
        p = _parse_int(nxt() or "")
        if p is None or not (1 <= p < u):
            bad()
            return
        parent[u] = p
    m = _parse_int(nxt() or "")
    if m is None or not (0 <= m <= MMAX):
        bad()
        return
    tree = Tree(n, parent, weight)
    for _ in range(m):
        k = _parse_int(nxt() or "")
        if k is None or not (0 <= k <= n):
            bad()
            return
        keys: List[int] = []
        seen = set()
        ok = True
        for _ in range(k):
            x = _parse_int(nxt() or "")
            if x is None or not (1 <= x <= n) or x in seen:
                ok = False
                break
            seen.add(x)
            keys.append(x)
        if not ok:
            bad()
            return
        cnt_nodes, edge_total, res = solve_query(tree, keys)
        parts = [str(cnt_nodes), str(edge_total)]
        for c, s in res:
            parts.append(str(c))
            parts.append(str(s))
        out.append(" ".join(parts))
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
