"""矩陣樹定理（Matrix-Tree Theorem / Kirchhoff 定理）：帶權生成樹計數

題意：
    給定一個 n 個點、m 條邊的無向帶權圖（允許重邊、不允許自環；
    權值可為 0，權值為 0 的邊等價於不存在）。
    定義一棵生成樹的權值為樹上所有邊權的**乘積**，
    求所有生成樹的權值**總和**，結果對 998244353 取模。
    圖不連通（不存在生成樹）時答案為 0。
    輸入：第一行 n m，接下來 m 行每行 u v w（1 起算，u != v）。
    輸出：一個整數。

思路：
    ### 1. Kirchhoff 拉普拉斯矩陣
    設度數矩陣 D（對角線 D[u][u] = 與 u 相連的所有邊權之和），
    鄰接矩陣 A（A[u][v] = u、v 之間的邊權之和，重邊累加），
    拉普拉斯矩陣 L = D - A。
    去掉 L 的任意一行一列（約定去掉第 0 行第 0 列）得到 n-1 階子式 M，
    帶權矩陣樹定理告訴我們：
        det(M) = 所有生成樹的「邊權乘積」之和
    無權圖取所有 w = 1，就退化成生成樹**個數**（Cayley 公式即其特例）。

    ### 2. 模質數意義下的行列式
    在模 p（p = 998244353，質數）意義下做高斯消元求行列式：
        * 對每一列，在當前行以下找一個非零主元；找不到則 det = 0
          （恰好對應「圖不連通 → 無生成樹」）。
        * 交換兩行，行列式變號。
        * det *= 主元；並用主元的逆元把下方各行的該列消成 0。
    全程只在整數（模 p）上運算，沒有任何捨入誤差——
    這正是不用實數消元的原因：生成樹計數本質上是整數，
    實數高斯消元在 n 稍大時就會因捨入完全失真。

    ### 3. 對拍
    內建測試用 O(2^m * m * alpha) 的暴力：
    枚舉所有 (n-1) 條邊的子集，用並查集判是否構成生成樹，
    是則把邊權乘積累加。與 Kirchhoff 行列式逐一比對。

複雜度：
    建矩陣     O(n + m)
    行列式消元 O(n^3)
    空間       O(n^2)

應用場景（Applications）：
    - **網路可靠性**：通訊網路的生成樹數量 = 連通骨架的冗餘度（可靠性多項式）
    - **電路分析**：Kirchhoff 矩陣的行列式 / 主子式給出等效電阻（電網、晶片電源網）
    - **圖論計數**：驗證 Cayley 公式（K_n 有 n^(n-2) 棵生成樹）、
      計數「拆掉哪些邊會讓網路不再有備援」
    - **馬爾可夫鏈**：生成樹計數與隨機遊走的平穩分佈（Markov chain tree theorem）
"""

import sys
from typing import List, Optional, Sequence, Tuple

MOD = 998244353
NMAX = 200
MMAX = 2000
WMAX = 10 ** 9


def det_mod(mat: List[List[int]]) -> int:
    """模 p 高斯消元求行列式（mat 為 k x k 方陣，會被就地修改）。"""
    k = len(mat)
    det = 1
    for col in range(k):
        # 找非零主元
        piv = -1
        for r in range(col, k):
            if mat[r][col] != 0:
                piv = r
                break
        if piv == -1:
            return 0
        if piv != col:
            mat[col], mat[piv] = mat[piv], mat[col]
            det = MOD - det  # 換行變號（MOD - 1 ≡ -1）
        det = det * mat[col][col] % MOD
        inv = pow(mat[col][col], MOD - 2, MOD)
        for r in range(col + 1, k):
            f = mat[r][col] * inv % MOD
            if f:
                row_r = mat[r]
                row_c = mat[col]
                for c in range(col, k):
                    row_r[c] = (row_r[c] - f * row_c[c]) % MOD
    return det % MOD


def spanning_tree_count(n: int, edges: Sequence[Tuple[int, int, int]]) -> int:
    """矩陣樹定理：回傳帶權生成樹總權值和（mod MOD）。

    edges 為 (u, v, w) 列表，u / v 為 0 起算且 u != v。
    """
    if n == 1:
        return 1  # 空積：只有一個點的圖恰有一棵「空生成樹」
    # 拉普拉斯矩陣
    L = [[0] * n for _ in range(n)]
    for u, v, w in edges:
        if w == 0:
            continue
        L[u][u] = (L[u][u] + w) % MOD
        L[v][v] = (L[v][v] + w) % MOD
        L[u][v] = (L[u][v] - w) % MOD
        L[v][u] = (L[v][u] - w) % MOD
    # 去掉第 0 行第 0 列
    minor = [row[1:] for row in L[1:]]
    return det_mod(minor)


class DSU:
    """並查集（僅供暴力對拍使用）。"""

    def __init__(self, n: int) -> None:
        self.fa = list(range(n))

    def find(self, x: int) -> int:
        while self.fa[x] != x:
            self.fa[x] = self.fa[self.fa[x]]
            x = self.fa[x]
        return x

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        self.fa[ra] = rb
        return True


def brute_force(n: int, edges: Sequence[Tuple[int, int, int]]) -> int:
    """暴力：枚舉所有 (n-1) 條邊的子集（m 很小時才可用）。"""
    m = len(edges)
    total = 0
    for mask in range(1 << m):
        # popcount 檢查
        if bin(mask).count("1") != n - 1:
            continue
        dsu = DSU(n)
        prod = 1
        ok = True
        for i in range(m):
            if mask >> i & 1:
                u, v, w = edges[i]
                prod = prod * w % MOD
                if not dsu.union(u, v):
                    ok = False
                    break
        if ok and prod:
            total = (total + prod) % MOD
    return total


# 自製 LCG（兩語言一致的隨機序列）
def lcg_next(state: int) -> Tuple[int, int]:
    state = (state * 1103515245 + 12345) & 0x7FFFFFFF
    return state, state


def run_tests() -> None:
    E = Tuple[int, int, int]
    # ---- 固定用例 ----
    # 1) 兩點一邊：帶權計數 = 邊權乘積 = 5
    assert spanning_tree_count(2, [(0, 1, 5)]) == 5
    # 2) 三角形帶權：去掉一條邊剩兩條，2*3 + 2*5 + 3*5 = 31
    assert spanning_tree_count(3, [(0, 1, 2), (1, 2, 3), (0, 2, 5)]) == 31
    # 3) 不連通
    assert spanning_tree_count(3, [(0, 1, 1)]) == 0
    # 4) 重邊：兩點兩條平行邊，權 1 和 3，各有 1 棵 → 4
    assert spanning_tree_count(2, [(0, 1, 1), (0, 1, 3)]) == 4
    # 5) 權為 0 的邊等價於不存在
    assert spanning_tree_count(2, [(0, 1, 0)]) == 0
    # 6) 單點圖：空生成樹，答案 1
    assert spanning_tree_count(1, []) == 1
    # 7) 四階完全圖 K4：Cayley 公式 4^2 = 16 棵生成樹
    k4: List[E] = []
    for u in range(4):
        for v in range(u + 1, 4):
            k4.append((u, v, 1))
    assert spanning_tree_count(4, k4) == 16

    # ---- 與暴力對拍（自製 LCG，兩語言同一序列）----
    state = 20261103

    def nxt(mod: int) -> Tuple[int, int]:
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state, state % mod

    for _ in range(300):
        n = nxt(4)[1] + 2  # 2..5
        m = nxt(n * 2)[1]  # 0..2n-2 條邊
        edges: List[E] = []
        for _ in range(m):
            u = nxt(n)[1]
            v = nxt(n)[1]
            if u == v:
                v = (v + 1) % n
            w = nxt(4)[1]
            edges.append((u, v, w))
        got = spanning_tree_count(n, edges)
        want = brute_force(n, edges)
        assert got == want, (n, edges, got, want)


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


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


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
    m_opt = _parse_int(nxt() or "")
    if m_opt is None or not (0 <= m_opt <= MMAX):
        bad()
        return
    n, m = n_opt, m_opt

    edges: List[Tuple[int, int, int]] = []
    for _ in range(m):
        u = _parse_int(nxt() or "")
        v = _parse_int(nxt() or "")
        w = _parse_int(nxt() or "")
        if u is None or v is None or w is None:
            bad()
            return
        if not (1 <= u <= n and 1 <= v <= n) or u == v or not (0 <= w <= WMAX):
            bad()
            return
        edges.append((u - 1, v - 1, w))

    emit(str(spanning_tree_count(n, edges)) + "\n")


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
