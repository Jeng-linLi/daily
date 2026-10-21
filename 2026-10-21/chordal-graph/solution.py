"""弦圖與完美消除序列：MCS / Lex-BFS 辨識，並在弦圖上線性時間求 ω、χ、α

題意：
    給一張 n 點 m 邊的**無向簡單圖**，回答若干詢問：

        mode 1   這張圖是不是**弦圖**（chordal / triangulated）？輸出 1 / 0
        mode 2   輸出 MCS（Maximum Cardinality Search）的訪問順序
        mode 3   輸出**完美消除序列**（PEO，= MCS 順序的反序）
        mode 4   ω：最大團的大小
        mode 5   χ：最小染色數，並給一組具體染色（每個點的顏色，0 起算）
        mode 6   α：最大獨立集的大小，並給一組具體的最大獨立集（點號升序）
        mode 7   輸出 PEO 中每個位置的「後續鄰居數」d_i（ω = 1 + max d_i）

思路：
    ### 1. 什麼是弦圖
    「長度 ≥ 4 的迴圈都**有弦**（有一條邊直接連起迴圈上兩個不相鄰的點）」，
    等價說法是「不存在 induced cycle of length ≥ 4」。樹、區間圖、分裂圖都是弦圖。

    弦圖最強的性質是它**有完美消除序列（PEO）**：
    一個排列 v1..vn，使得每個 vi 在「後面」（v_{i+1}..v_n）的鄰居**構成一個團**。

    ### 2. 兩個產生 PEO 的演算法（互相對拍）
    - **MCS**：每次選「已訪問鄰居數最多」的未訪問點（平手取編號最小）。
    - **Lex-BFS**：維護一個「集合序列」，每次取第一個集合的第一個點，
      再把每個集合依「是否與該點相鄰」切成兩半（不相鄰的在前）。

    兩者的**反序**都是 PEO —— 但前提是圖真的是弦圖。
    所以「把候選順序拿去驗 PEO」本身就是一個**弦圖辨識演算法**：
    G 是弦圖 ⟺ reverse(MCS order) 是 PEO。

    ### 3. 為什麼要用兩套
    MCS 與 Lex-BFS 的思路完全不同（一個數鄰居、一個維護字典序分區），
    卻必須對所有圖給出**相同**的「是否弦圖」答案。這是本題最主要的交叉驗證。
    另外在 n ≤ 11 時再用「枚舉所有點子集、找 induced cycle」的 2^n 暴力法做第三方確認。

    ### 4. 弦圖上三個 NP-hard 問題全部變線性
    一般圖上求 ω（最大團）、χ（最小染色）、α（最大獨立集）都是 NP-hard，
    但拿到 PEO 之後：

        ω = 1 + max_i |{v_i 在 PEO 中後面的鄰居}|         （那個集合本身就是團）
        χ：依 PEO 的**反序**做貪心染色（每次選最小的可用顏色）→ 恰好 ω 種顏色
        α：依 PEO 的**正序**做貪心（與已選集合不相鄰就選）→ 最大獨立集

    因為弦圖是**完美圖**（perfect graph），χ = ω；而 α = θ（最小團覆蓋數）。
    本題會把 χ 與 ω 的相等性當成一條不變式來檢查。

    ### 5. 兩語言一致的細節
    - 平手規則全部寫死：MCS 取編號最小；Lex-BFS 的每個集合內部維持編號升序、取第一個。
    - 非弦圖時 ω / α 的 PEO 公式只是**下界**（染色仍是合法但不一定最少），
      這時只檢查「≤ 暴力值」與「染色合法」，不檢查相等。
    - 輸入的邊先正規化：自環、越界、重邊一律丟棄，兩語言才會看到同一張圖。

應用場景：
    - 稀疏矩陣 Cholesky 分解：弦圖 = 完美消除順序 = 分解時不產生新的非零元（fill-in）
    - 貝氏網路 / 概率圖模型：junction tree 演算法必須先把圖「三角化」成弦圖
    - 資料庫：acyclic hypergraph（α-acyclic）與 join tree，等價於弦圖性質
    - 約束滿足與排程：區間圖（弦圖的子類）上的教室排課、頻道分配是多項式時間
    - 電腦圖學：三角化多邊形（polygon triangulation）本質上就是弦圖

複雜度（n 點 m 邊）：
    MCS / Lex-BFS（本題的樸素實作）     O(n²)
    PEO 驗證                            O(n + m)
    ω / χ / α                           O(n + m)
    暴力對拍（n ≤ 14 / n ≤ 11）          O(2^n · n)
    空間                                O(n + m)
"""

import random
import re
import sys
from collections import deque
from typing import Deque, Dict, List, Set, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 60           # IO 模式下的點數上限
MMAX = 1200         # IO 模式下的邊數上限
QMAX = 20           # IO 模式下的詢問數上限
BRUTE_N = 14        # n 不超過這個值才做 2^n 最大團 / 獨立集對拍
BRUTE_CYCLE_N = 11  # n 不超過這個值才做 2^n induced cycle 對拍


# ------------------------------------------------------------------ 圖


class Graph:
    def __init__(self, n: int) -> None:
        self.n = n
        self.adj: List[Set[int]] = [set() for _ in range(n)]
        self.amask: List[int] = [0] * n       # 鄰接位元遮罩（0-based）

    def add_edge(self, u: int, v: int) -> None:
        if u == v or not (0 <= u < self.n) or not (0 <= v < self.n):
            return
        self.adj[u].add(v)
        self.adj[v].add(u)

    def finish(self) -> None:
        for i in range(self.n):
            m = 0
            for j in self.adj[i]:
                m |= 1 << j
            self.amask[i] = m


# ------------------------------------------------------------ 順序產生


def mcs_order(g: Graph) -> List[int]:
    """Maximum Cardinality Search：每次選已訪問鄰居最多的未訪問點。"""
    n = g.n
    visited = [False] * n
    cnt = [0] * n
    order: List[int] = []
    for _ in range(n):
        best = -1
        for v in range(n):
            if not visited[v] and (best < 0 or cnt[v] > cnt[best]):
                best = v               # 嚴格大於才換人 ⟹ 平手取編號最小
        visited[best] = True
        order.append(best)
        for w in sorted(g.adj[best]):
            cnt[w] += 1
    return order


def lexbfs_order(g: Graph) -> List[int]:
    """Lexicographic BFS：集合序列分區，每次取第一個集合的第一個點。"""
    labels: List[List[int]] = [sorted(range(g.n))] if g.n else []
    order: List[int] = []
    picked = [False] * g.n
    while labels:
        s = labels.pop(0)
        v = s.pop(0)
        if s:
            labels.insert(0, s)
        picked[v] = True
        order.append(v)
        new_labels: List[List[int]] = []
        nb = g.adj[v]
        for t in labels:
            t = [x for x in t if not picked[x]]
            if not t:
                continue
            # 與 v 相鄰者標籤較大（字典序較大），必須排在**前面**
            a = [x for x in t if x in nb]
            b = [x for x in t if x not in nb]
            if a:
                new_labels.append(a)
            if b:
                new_labels.append(b)
        labels = new_labels
    return order


# ------------------------------------------------------------------ PEO


def peo_later_counts(g: Graph, peo: List[int]) -> List[int]:
    """PEO 中每個位置的「後續鄰居數」。"""
    n = g.n
    pos = [0] * n
    for i, v in enumerate(peo):
        pos[v] = i
    out = [0] * n
    for i, v in enumerate(peo):
        c = 0
        for w in g.adj[v]:
            if pos[w] > i:
                c += 1
        out[i] = c
    return out


def is_peo(g: Graph, peo: List[int]) -> bool:
    """每個點在序列後面的鄰居必須構成團。"""
    n = g.n
    pos = [0] * n
    for i, v in enumerate(peo):
        pos[v] = i
    for i, v in enumerate(peo):
        later = sorted((w for w in g.adj[v] if pos[w] > i), key=lambda w: pos[w])
        for a in range(len(later)):
            for b in range(a + 1, len(later)):
                if later[b] not in g.adj[later[a]]:
                    return False
    return True


def max_clique_peo(g: Graph, peo: List[int]) -> int:
    d = peo_later_counts(g, peo)
    return (1 + max(d)) if d else 0


def color_by_peo(g: Graph, peo: List[int]) -> List[int]:
    """依 PEO 的**反序**貪心染色，每次選最小的可用顏色。"""
    n = g.n
    color = [-1] * n
    for i in range(n - 1, -1, -1):
        v = peo[i]
        used = set()
        for w in g.adj[v]:
            if color[w] >= 0:
                used.add(color[w])
        c = 0
        while c in used:
            c += 1
        color[v] = c
    return color


def mis_by_peo(g: Graph, peo: List[int]) -> List[int]:
    """依 PEO 的**正序**貪心取點：與已選集合不相鄰就選。"""
    chosen: List[int] = []
    taken = set()
    for v in peo:
        if any(w in taken for w in g.adj[v]):
            continue
        chosen.append(v)
        taken.add(v)
    return chosen


# ------------------------------------------------------------- 暴力對照


def brute_max_clique(g: Graph) -> int:
    n = g.n
    best = 0
    for mask in range(1, 1 << n):
        ok = True
        m = mask
        while m:
            b = m & -m
            v = b.bit_length() - 1
            if (g.amask[v] | b) & mask != mask:
                ok = False
                break
            m -= b
        if ok:
            c = bin(mask).count("1")
            if c > best:
                best = c
    return best


def brute_mis(g: Graph) -> int:
    """最大獨立集 = 補圖的最大團。"""
    n = g.n
    full = (1 << n) - 1
    best = 0
    for mask in range(1, 1 << n):
        ok = True
        m = mask
        while m:
            b = m & -m
            v = b.bit_length() - 1
            if (g.amask[v] & mask) != 0:
                ok = False
                break
            m -= b
        if ok:
            c = bin(mask).count("1")
            if c > best:
                best = c
    return best


def brute_is_chordal(g: Graph) -> bool:
    """枚舉所有點子集，檢查是否存在 induced cycle（長度 ≥ 4）。"""
    n = g.n
    for mask in range(1, 1 << n):
        c = bin(mask).count("1")
        if c < 4:
            continue
        ok = True
        m = mask
        while m:
            b = m & -m
            v = b.bit_length() - 1
            if bin(g.amask[v] & mask).count("1") != 2:
                ok = False
                break
            m -= b
        if not ok:
            continue
        # 所有點度數都是 2，再確認連通 ⟹ 是一個 induced cycle
        start = (mask & -mask).bit_length() - 1
        seen = 1 << start
        q: Deque[int] = deque([start])
        while q:
            x = q.popleft()
            for w in sorted(g.adj[x]):
                if (mask >> w) & 1 and not (seen >> w) & 1:
                    seen |= 1 << w
                    q.append(w)
        if seen == mask:
            return False
    return True


def coloring_proper(g: Graph, color: List[int]) -> bool:
    for v in range(g.n):
        for w in g.adj[v]:
            if color[v] == color[w]:
                return False
    return True


def independent_valid(g: Graph, s: List[int]) -> bool:
    st = set(s)
    for v in s:
        for w in g.adj[v]:
            if w in st:
                return False
    return True


# ------------------------------------------------------------------ 驅動


def run(n: int, edges: List[Tuple[int, int]], modes: List[int]) -> str:
    g = Graph(n)
    for u, v in edges:
        g.add_edge(u, v)
    g.finish()

    mcs = mcs_order(g)
    lex = lexbfs_order(g)
    peo = list(reversed(mcs))
    peo2 = list(reversed(lex))

    chordal = is_peo(g, peo)
    chordal2 = is_peo(g, peo2)
    assert chordal == chordal2, ("chordal mismatch", n, edges, chordal, chordal2)

    if n <= BRUTE_CYCLE_N:
        bc = brute_is_chordal(g)
        assert chordal == bc, ("brute chordal mismatch", n, edges, chordal, bc)

    d = peo_later_counts(g, peo)
    omega = (1 + max(d)) if d else 0
    color = color_by_peo(g, peo)
    chi = (max(color) + 1) if n else 0
    mis = mis_by_peo(g, peo)
    alpha = len(mis)

    assert coloring_proper(g, color), ("coloring invalid", n, edges)
    assert independent_valid(g, mis), ("MIS invalid", n, edges)

    if n <= BRUTE_N:
        bo = brute_max_clique(g)
        ba = brute_mis(g)
        # α 一定是合法獨立集，任何圖都 ≤ 暴力值；ω 的 PEO 公式只在弦圖才成立
        assert alpha <= ba, ("alpha > brute", n, edges, alpha, ba)
        assert chi >= bo, ("chi < clique", n, edges, chi, bo)
        if chordal:
            assert omega == bo, ("chordal omega", n, edges, omega, bo)
            assert alpha == ba, ("chordal alpha", n, edges, alpha, ba)
            assert chi == bo, ("chordal chi", n, edges, chi, bo)

    out: List[str] = []
    for m in modes:
        if m == 1:
            out.append("1" if chordal else "0")
        elif m == 2:
            out.append(" ".join(str(v + 1) for v in mcs))
        elif m == 3:
            out.append(" ".join(str(v + 1) for v in peo))
        elif m == 4:
            out.append(str(omega))
        elif m == 5:
            s = str(chi)
            for v in range(n):
                s += " %d" % color[v]
            out.append(s)
        elif m == 6:
            s = str(alpha)
            for v in sorted(mis):
                s += " %d" % (v + 1)
            out.append(s)
        elif m == 7:
            out.append(" ".join(str(x) for x in d))
        else:
            out.append("error")
    return "".join(s + "\n" for s in out)


def _check(n: int, edges: List[Tuple[int, int]], modes: List[int]) -> None:
    run(n, edges, modes)


ALL_MODES = [1, 2, 3, 4, 5, 6, 7]


def run_tests() -> None:
    # 固定用例 1：路徑 1-2-3-4（樹 ⟹ 弦圖）
    _check(4, [(0, 1), (1, 2), (2, 3)], ALL_MODES)

    # 固定用例 2：星形（中心 1，葉子 2..5）
    _check(5, [(0, 1), (0, 2), (0, 3), (0, 4)], ALL_MODES)

    # 固定用例 3：4-迴圈 1-2-3-4-1（無弦 ⟹ 不是弦圖）
    _check(4, [(0, 1), (1, 2), (2, 3), (3, 0)], ALL_MODES)

    # 固定用例 4：K4（完全圖 ⟹ 弦圖）
    _check(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], ALL_MODES)

    # 固定用例 5：兩個三角形共用一個點
    _check(5, [(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (2, 4)], ALL_MODES)

    # 固定用例 6：空圖 / 單點 / 自環 / 重邊 / 越界
    _check(0, [], ALL_MODES)
    _check(1, [], ALL_MODES)
    _check(3, [(0, 0), (0, 1), (0, 1), (0, 5)], ALL_MODES)

    # 固定用例 7：5-迴圈加一條弦（變成弦圖）
    _check(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 2)], ALL_MODES)

    rnd = random.Random(20261021)

    # 隨機樹（一定弦圖）
    for _ in range(60):
        n = rnd.randint(1, 12)
        edges: List[Tuple[int, int]] = []
        for v in range(1, n):
            edges.append((v, rnd.randint(0, v - 1)))
        _check(n, edges, ALL_MODES)

    # 隨機一般圖
    for _ in range(120):
        n = rnd.randint(0, 12)
        p = rnd.choice([0.1, 0.25, 0.4, 0.6])
        edges = []
        for u in range(n):
            for v in range(u + 1, n):
                if rnd.random() < p:
                    edges.append((u, v))
        _check(n, edges, ALL_MODES)

    # 隨機「分裂圖」：團 + 獨立集（一定弦圖）
    for _ in range(40):
        n = rnd.randint(2, 12)
        k = rnd.randint(1, n)
        edges = []
        for u in range(k):
            for v in range(u + 1, k):
                edges.append((u, v))
        for v in range(k, n):
            for u in range(k):
                if rnd.random() < 0.5:
                    edges.append((u, v))
        _check(n, edges, ALL_MODES)

    # 隨機含雜訊邊的輸入（自環 / 越界 / 重邊）
    for _ in range(40):
        n = rnd.randint(2, 10)
        edges = []
        for _ in range(rnd.randint(1, 12)):
            edges.append((rnd.randint(0, n + 1), rnd.randint(0, n + 1)))
        _check(n, edges, ALL_MODES)


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt(default: int = 0) -> int:
        nonlocal pos
        if pos >= len(toks):
            return default
        t = toks[pos]
        pos += 1
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-(10 ** 15), min(10 ** 15, v))

    n = nxt(0)
    if n < 0:
        n = 0
    if n > NMAX:
        n = NMAX
    m = nxt(0)
    if m < 0:
        m = 0
    if m > MMAX:
        m = MMAX
    edges: List[Tuple[int, int]] = []
    for _ in range(m):
        # 輸入是 1-based，轉成 0-based；越界 / 自環由 add_edge 過濾
        a = nxt(1) - 1
        b = nxt(1) - 1
        edges.append((a, b))
    q = nxt(0)
    if q < 0:
        q = 0
    if q > QMAX:
        q = QMAX
    modes = [nxt(0) for _ in range(q)]
    emit(run(n, edges, modes))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
