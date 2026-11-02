"""線段樹合併（Segment Tree Merge）：樹上子樹的相異值統計

題意：
    給定一棵有 n 個節點的**根樹**（根為 1），每個節點 u 有一個整數標籤 c[u]。
    對每個節點 u，求其子樹中出現過的**相異**（distinct）標籤數量。
    限制：1 <= n <= 100000，|c[u]| <= 10^9。

思路：
    ### 1. 為什麼不能用樸素做法
    最直接的想法是「對每個節點都遍歷一次子樹，把標籤塞進 set」。
    但一個鏈狀樹上，子樹大小總和是 O(n^2)，直接超時。
    改進一點：把小集合併進大集合（启发式合并 / small-to-large），
    可以做到 O(n log^2 n)，但仍然要維護 hash 集合，常數不小。

    ### 2. 線段樹合併的核心想法
    換一個視角：**自底向上，把子樹的資訊「合併」到父節點**。
    對每個節點 u 維護一棵「值域線段樹」，代表 u 子樹裡出現過的所有標籤：
        葉子（某個具體標籤值）= 該標籤是否出現（0/1）
        內部節點 = 左右兩邊的出現數之和
    於是 u 子樹的相異標籤數 = u 那棵線段樹的根節點的 cnt。

    合併兩棵線段樹 merge(a, b) 的規則非常簡單：
        若 a 為空 → 直接回傳 b
        若 b 為空 → 直接回傳 a
        否則（兩者都非空，這是最後一種情形也是唯一需要遞迴的情形）：
            a 的左子樹 = merge(a.left, b.left)
            a 的右子樹 = merge(a.right, b.right)
            a.cnt = a.left.cnt + a.right.cnt
    底層標籤是「出現/不出現」，所以葉子合併時取 or（1）而不是相加 ——
    這正是「相異」的語意。

    ### 3. 為什麼總複雜度是 O(n log V) 而不是 O(n^2)
    merge 裡唯一會「消滅」節點的情形是「a、b 都非空」；
    此時 b 的根節點被丟棄（不再被任何線段樹引用）。
    每次 merge 的「有消耗」呼叫都會消滅一個節點，
    而全程建立的節點總數不超過 O(n log V)（每個節點插入一個標籤，
    最多開 log V 個節點）。
    所以所有 merge 呼叫加總是 O(節點總數) = O(n log V)。
    這是線段樹合併最漂亮的地方：**它不需要任何平衡技巧**。

    ### 4. 實作順序（後序遍歷）
    用顯式堆疊做 DFS 取得「父在子前」的順序，再倒過來就是「子在父前」。
    對每個節點 u（依此順序）：
        1. 把 u 自己的標籤插入 root[u]（此時 root[u] 已含所有子樹的貢獻）
        2. ans[u] = cnt[root[u]]
        3. 把 root[u] 合併進 root[parent[u]]
    注意 ans[u] 必須在「合併進父親」之前記錄，
    因為合併會改寫（甚至換掉）那些節點。

複雜度：
    離散化            O(n log n)
    建立 + 插入       O(n log V)
    所有 merge 合計    O(n log V)（節點消耗 amortize）
    總計              時間 O(n log n)，空間 O(n log V)

應用場景（Applications）：
    - **專案成員技能盤點**：每個節點是一個人，子樹是團隊，
      「團隊涵蓋了幾種技能」就是子樹的相異值數
    - **物種普查**：把觀測點依地理層級組成樹，
      每個區域的物種豐富度（species richness）就是這個問題
    - **授權 / 權限盤點**：群組樹上的「有效權限種類數」
    - **程式碼重複偵測**：以 AST 節點組成樹，統計子樹內的相異節點型別數
"""

import sys
from array import array
from typing import List, Optional, Sequence, Tuple

NMAX = 100000
VAL_MAX = 10 ** 9


class SegmentTreePool:
    """動態開點的線段樹節點池（用平行陣列存，避免大量 Python 物件）。

    節點 0 代表「空節點」（cnt 恆為 0）。
    """

    __slots__ = ("left", "right", "cnt")

    def __init__(self) -> None:
        self.left = array("i", [0])
        self.right = array("i", [0])
        self.cnt = array("i", [0])

    def new_node(self) -> int:
        self.left.append(0)
        self.right.append(0)
        self.cnt.append(0)
        return len(self.cnt) - 1

    def insert(self, node: int, pos: int, lo: int, hi: int) -> int:
        """在 pos 處插入標籤（存在則維持 1，即集合語意）。回傳新根。"""
        left = self.left
        right = self.right
        cnt = self.cnt
        if node == 0:
            node = self.new_node()
        if lo == hi:
            cnt[node] = 1
            return node
        mid = (lo + hi) >> 1
        if pos <= mid:
            left[node] = self.insert(left[node], pos, lo, mid)
        else:
            right[node] = self.insert(right[node], pos, mid + 1, hi)
        cnt[node] = cnt[left[node]] + cnt[right[node]]
        return node

    def merge(self, a: int, b: int, lo: int, hi: int) -> int:
        """合併兩棵線段樹（集合聯集）。回傳合併後的根。"""
        if a == 0:
            return b
        if b == 0:
            return a
        left = self.left
        right = self.right
        cnt = self.cnt
        if lo == hi:
            # 同一個標籤值：出現即為 1（相異，不是計數）
            cnt[a] = 1
            return a
        mid = (lo + hi) >> 1
        left[a] = self.merge(left[a], left[b], lo, mid)
        right[a] = self.merge(right[a], right[b], mid + 1, hi)
        cnt[a] = cnt[left[a]] + cnt[right[a]]
        return a

    def count(self, node: int) -> int:
        """根節點的 cnt = 該樹中的相異標籤數。"""
        return self.cnt[node] if node else 0


def subtree_distinct_counts(n: int, labels: Sequence[int],
                            edges: Sequence[Tuple[int, int]]) -> List[int]:
    """求每個節點子樹內的相異標籤數。

    參數：
        n      —— 節點數（1..n）
        labels —— 每個節點的標籤（索引 0 起算，長度 n）
        edges  —— n-1 條無向邊 ((u, v), ...)
    回傳：長度 n 的列表，第 i 項是節點 i+1 的答案。
    """
    # ---- 1. 標籤離散化 ----
    uniq = sorted(set(labels))
    m = len(uniq)
    rank_of = {v: i + 1 for i, v in enumerate(uniq)}  # 1-based 值域

    # ---- 2. 建鄰接表 ----
    adj: List[List[int]] = [[] for _ in range(n + 1)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # ---- 3. 取得「父在子前」的遍歷順序（顯式堆疊，避免遞迴爆栈）----
    parent = [0] * (n + 1)
    seen = bytearray(n + 1)
    order: List[int] = []
    stack = [1]
    seen[1] = 1
    parent[1] = -1
    while stack:
        u = stack.pop()
        order.append(u)
        for v in adj[u]:
            if not seen[v]:
                seen[v] = 1
                parent[v] = u
                stack.append(v)

    # ---- 4. 自底向上：插入 + 記錄答案 + 合併到父親 ----
    pool = SegmentTreePool()
    root = [0] * (n + 1)
    ans = [0] * (n + 1)
    for u in reversed(order):
        r = pool.insert(root[u], rank_of[labels[u - 1]], 1, m)
        root[u] = r
        ans[u] = pool.count(r)
        p = parent[u]
        if p > 0:
            root[p] = pool.merge(root[p], r, 1, m)
    return ans[1:]


def brute_force(n: int, labels: Sequence[int],
                edges: Sequence[Tuple[int, int]]) -> List[int]:
    """參考實作：對每個節點遍歷子樹收集標籤去重（O(n^2)，僅供對拍）。"""
    adj: List[List[int]] = [[] for _ in range(n + 1)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    # 先求出以 1 為根的父子關係，樸素遍歷只能往下走（不能回頭走到父親）
    parent = [0] * (n + 1)
    seen = bytearray(n + 1)
    seen[1] = 1
    stack = [1]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if not seen[v]:
                seen[v] = 1
                parent[v] = u
                stack.append(v)
    children: List[List[int]] = [[] for _ in range(n + 1)]
    for v in range(2, n + 1):
        children[parent[v]].append(v)

    ans = [0] * (n + 1)
    for start in range(1, n + 1):
        vals = set()
        stack = [start]
        while stack:
            u = stack.pop()
            vals.add(labels[u - 1])
            stack.extend(children[u])
        ans[start] = len(vals)
    return ans[1:]


def demo_resource_pool_merge() -> None:
    """應用程式示範：把「資源池合併」抽象成同一套 merge 機制。

    場景：三個團隊各自持有一批資源編號，合併後要問「這個聯集有幾種資源」。
    這在權限盤點、標籤聯集、特徵庫合併等場景都會用到。
    """
    pool = SegmentTreePool()
    # 池 A = {1, 5, 9}，池 B = {5, 7}，池 C = {9, 11}
    sets = [(1, 5, 9), (5, 7), (9, 11)]
    roots = []
    for s in sets:
        r = 0
        for v in s:
            r = pool.insert(r, v, 1, 16)
        roots.append(r)
    assert pool.count(roots[0]) == 3
    assert pool.count(roots[1]) == 2
    assert pool.count(roots[2]) == 2

    # A ← A ∪ B
    merged_ab = pool.merge(roots[0], roots[1], 1, 16)
    assert pool.count(merged_ab) == 4          # {1,5,9,7}
    # 再併 C → {1,5,7,9,11}
    merged_abc = pool.merge(merged_ab, roots[2], 1, 16)
    assert pool.count(merged_abc) == 5
    # 合併同一個池子兩次不會重複計數（集合語意）
    again = pool.merge(merged_abc, merged_abc, 1, 16)
    assert pool.count(again) == 5


def run_tests() -> None:
    # ---- 固定用例 ----
    # 1) 單節點
    assert subtree_distinct_counts(1, [7], []) == [1]
    # 2) 鏈狀樹，全部不同 → 每個節點答案都是子樹大小
    #    1 - 2 - 3
    assert subtree_distinct_counts(3, [1, 2, 3], [(1, 2), (2, 3)]) == [3, 2, 1]
    # 3) 鏈狀樹，標籤全部相同 → 每個答案都是 1
    assert subtree_distinct_counts(3, [5, 5, 5], [(1, 2), (2, 3)]) == [1, 1, 1]
    # 4) 星狀樹：根 + 3 個葉子，其中兩個葉子同標籤
    #    1 = 根(9)，2/3/4 = 葉(1, 2, 1)
    #    根子樹 = {9,1,2} → 3；葉子各自 = 1
    assert subtree_distinct_counts(4, [9, 1, 2, 1],
                                   [(1, 2), (1, 3), (1, 4)]) == [3, 1, 1, 1]
    # 5) 兩個分支各貢獻不同標籤
    #    1(根) -> 2(子) -> 4(孫)，1 -> 3(子) -> 5(孫)
    #    標籤：1:'a', 2:'b', 3:'c', 4:'a', 5:'d'
    #    ans[1] = {a,b,c,d} = 4；ans[2] = {b,a} = 2；ans[3] = {c,d} = 2
    assert subtree_distinct_counts(5, [1, 2, 3, 1, 4],
                                   [(1, 2), (2, 4), (1, 3), (3, 5)]) == [4, 2, 2, 1, 1]
    # 6) 負數與 0 標籤
    assert subtree_distinct_counts(3, [-1, 0, -1],
                                   [(1, 2), (1, 3)]) == [2, 1, 1]
    # 7) 單一標籤的所有節點
    assert subtree_distinct_counts(4, [8, 8, 8, 8],
                                   [(1, 2), (2, 3), (3, 4)]) == [1, 1, 1, 1]
    # 8) 完全平衡的二叉樹，標籤對稱
    #      1
    #    2   3
    #   4 5  6 7
    #    標籤：1:1 2:2 3:3 4:4 5:2 6:3 7:4
    #    ans[1]={1,2,3,4}=4, ans[2]={2,4}=2, ans[3]={3,4}=2, 葉子=1
    assert subtree_distinct_counts(7, [1, 2, 3, 4, 2, 3, 4],
                                   [(1, 2), (1, 3), (2, 4), (2, 5),
                                    (3, 6), (3, 7)]) == [4, 2, 2, 1, 1, 1, 1]
    # 9) 應用程式示範
    demo_resource_pool_merge()

    # ---- 與參考實作對拍 ----
    import random

    rng = random.Random(20261102)
    for _ in range(200):
        n = rng.randint(1, 10)
        # 小值域 → 經常有重複標籤，能壓到 merge 的「葉子合併」分支
        labels = [rng.randint(-3, 3) for _ in range(n)]
        edges = []
        for v in range(2, n + 1):
            u = rng.randint(1, v - 1)
            edges.append((u, v))
        got = subtree_distinct_counts(n, labels, edges)
        want = brute_force(n, labels, edges)
        assert got == want, (n, labels, edges, got, want)

    # 大值域 → 幾乎沒有重複，壓到 merge 的「內部節點合併」分支
    for _ in range(150):
        n = rng.randint(1, 12)
        labels = [rng.randint(-1000, 1000) for _ in range(n)]
        edges = []
        for v in range(2, n + 1):
            u = rng.randint(1, v - 1)
            edges.append((u, v))
        got = subtree_distinct_counts(n, labels, edges)
        want = brute_force(n, labels, edges)
        assert got == want, (n, labels, edges, got, want)

    # 隨機樹（不保證是鏈，節點編號與父子關係隨機）
    for _ in range(150):
        n = rng.randint(1, 12)
        labels = [rng.randint(-2, 2) for _ in range(n)]
        edges = []
        for v in range(2, n + 1):
            u = rng.randint(1, v - 1)
            edges.append((u, v))
        rng.shuffle(edges)
        got = subtree_distinct_counts(n, labels, edges)
        want = brute_force(n, labels, edges)
        assert got == want, (n, labels, edges, got, want)

    # 較大規模：深鏈 + 大量重複
    n = 400
    labels = [rng.randint(0, 5) for _ in range(n)]
    edges = [(i, i + 1) for i in range(1, n)]
    got = subtree_distinct_counts(n, labels, edges)
    want = brute_force(n, labels, edges)
    assert got == want


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
    n = n_opt

    labels: List[int] = []
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None or abs(v) > VAL_MAX:
            bad()
            return
        labels.append(v)

    edges: List[Tuple[int, int]] = []
    for _ in range(n - 1):
        u = _parse_int(nxt() or "")
        v = _parse_int(nxt() or "")
        if u is None or v is None or not (1 <= u <= n) or not (1 <= v <= n):
            bad()
            return
        edges.append((u, v))

    for a in subtree_distinct_counts(n, labels, edges):
        out.append(str(a))
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
