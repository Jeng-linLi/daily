"""最大團與 Bron–Kerbosch 演算法（含貪心著色上界與補圖獨立集）

題意：
    給一張 n 點的無向簡單圖（鄰接矩陣），求：

    1. **最大團大小**（maximum clique）：最大的「兩兩相鄰」的點集有多大
    2. **極大團的數量**（maximal cliques）：不能再加點的團有幾個
    3. **字典序最小的那個最大團**的點集合
    4. **最大獨立集**與**最小點覆蓋**
    5. 全部結果都與「枚舉 2^n 個子集」的暴力版對拍

思路：
    ### 1. Bron–Kerbosch：列舉所有極大團
    維護三個集合：
        R = 目前這條搜尋路徑上已經選進團裡的點
        P = 還能加進來的候選點（與 R 中所有點都相鄰）
        X = 已經處理過、用來避免重複列舉同一個極大團的點
    當 P 與 X 同時為空，R 就是一個**極大團**（加不進任何點了）。

    ### 2. 樞軸（pivot）剪枝
    沒有樞軸的版本會對同一個極大團重複搜尋很多次。
    取 u ∈ P ∪ X 中「|P ∩ N(u)| 最大」的那個當樞軸（平手取編號最小的）：
    任何極大團要嘛包含 u、要嘛包含某個「與 u 不相鄰」的點，
    所以只需遞迴 **P \\ N(u)** 裡的點，候選數目明顯變小，且不漏解。

    ### 3. 貪心著色上界（branch and bound）
    只求「大小」不求「列舉」時可以剪得更兇。對候選集合 P 做一次貪心著色：
    依序給每個點塗上「鄰居沒用過的最小顏色」，得到 color[v]。
    因為同一個團裡任兩點相鄰，顏色必互不相同，所以
        |R| + color[v] ≤ 目前最佳解  →  從 v 出發不可能更好，直接剪掉
    這就是 Tomita 的 MCQ 演算法，實務上比純樞軸版快非常多。

    ### 4. 最大獨立集 = 補圖的最大團
    原圖的獨立集（兩兩**不**相鄰）正好是補圖的團（兩兩相鄰）。
    而對任意圖都有「α + τ = n」（最大獨立集 + 最小點覆蓋 = 點數），
    因為點覆蓋的補集一定是獨立集，反之亦然。

    ### 5. 跨語言一致的排序
    著色的順序、樞軸的平手規則、遞迴枚舉順序全部寫死成「編號由小到大」，
    否則 Python 與 C++ 的集合迭代順序不同，剪枝量就不同，輸出對不起來。

應用場景：
    - **社群網路**：找出「一群互相認識的人」（ clique 就是最緊密的社群結構）
    - **生物資訊**：蛋白質交互作用網路裡找功能模組；基因共表現聚類
    - **編譯器**：暫存器分配時的干涉圖、指令排程的衝突圖
    - **化學資訊學**：用最大團做分子結構比對（兩個分子圖的最大公共子結構）
    - 排班 / 會議安排：獨立集 = 可以同時進行的任務集合

複雜度：
    樞軸版列舉所有極大團   O(3^{n/3}) 時間（最壞情況，且是緊的）
    著色上界版求最大團     最壞同階，實務上遠快
    空間                   O(n²) 鄰接矩陣 + O(n) 遞迴深度
    暴力版（n ≤ 14）       O(2^n · n²) 時間
"""

import random
import re
import sys
from typing import List, Set, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 20            # IO 模式下的最大點數
MAXM = 190           # 最大邊數（簡單圖上限 n(n-1)/2）
BRUTE_MAX = 14       # 頂點數不超過這個值才跑暴力對拍


# ---------------------------------------------------------------- 工具


def build_adj(n: int, edges: List[Tuple[int, int]]) -> List[List[int]]:
    adj = [[0] * n for _ in range(n)]
    for u, v in edges:
        if 0 <= u < n and 0 <= v < n and u != v:
            adj[u][v] = 1
            adj[v][u] = 1
    return adj


def complement(adj: List[List[int]]) -> List[List[int]]:
    n = len(adj)
    return [[(1 if i != j and not adj[i][j] else 0) for j in range(n)] for i in range(n)]


# ---------------------------------------------------------------- 演算法一：樞軸版（列舉極大團）


def bron_kerbosch(adj: List[List[int]]) -> Tuple[int, int, List[int]]:
    """回傳 (最大團大小, 極大團數量, 字典序最小的最大團)。"""
    n = len(adj)
    best = 0
    count = 0
    best_set: List[int] = []

    def rec(r: List[int], p: List[int], x: List[int]) -> None:
        nonlocal best, count, best_set
        if not p and not x:
            count += 1
            if len(r) > best:
                best = len(r)
                best_set = sorted(r)
            elif len(r) == best and best > 0:
                cand = sorted(r)
                if not best_set or cand < best_set:
                    best_set = cand
            return
        # 樞軸：P ∪ X 中與 P 交集最大者，平手取編號最小
        pu = -1
        best_cover = -1
        for u in p:
            c = sum(1 for w in p if adj[u][w])
            if c > best_cover or (c == best_cover and (pu < 0 or u < pu)):
                best_cover = c
                pu = u
        for u in x:
            c = sum(1 for w in p if adj[u][w])
            if c > best_cover or (c == best_cover and (pu < 0 or u < pu)):
                best_cover = c
                pu = u
        if pu < 0:
            return
        cand = [v for v in p if not adj[pu][v]]
        for v in cand:
            np_ = [w for w in p if adj[v][w]]
            nx_ = [w for w in x if adj[v][w]]
            rec(r + [v], np_, nx_)
            p = [w for w in p if w != v]
            x = x + [v]

    rec([], list(range(n)), [])
    return (best, count, best_set)


# ---------------------------------------------------------------- 演算法二：著色上界版


def greedy_color(adj: List[List[int]], p: List[int]) -> List[int]:
    """對候選集 p 依編號順序貪心著色，回傳每個位置的顏色編號（由 1 起算）。"""
    n = len(adj)
    color_of = [-1] * n
    out = [0] * len(p)
    for idx, v in enumerate(p):
        used = [False] * (len(p) + 2)
        for u in p:
            if color_of[u] > 0 and adj[v][u]:
                used[color_of[u]] = True
        c = 1
        while c <= len(p) + 1 and used[c]:
            c += 1
        color_of[v] = c
        out[idx] = c
    return out


def max_clique_color(adj: List[List[int]]) -> int:
    """Tomita 風格：貪心著色上界 + 分支限界。回傳最大團大小。"""
    best = 0

    def expand(cur: int, p: List[int], colors: List[int]) -> None:
        nonlocal best
        if not p:
            if cur > best:
                best = cur
            return
        # 貪心著色不保證 colors 隨下標遞增，所以上界要用「前綴最大值」
        pmax = [0] * len(p)
        acc = 0
        for i, c in enumerate(colors):
            if c > acc:
                acc = c
            pmax[i] = acc
        for i in range(len(p) - 1, -1, -1):
            if cur + pmax[i] <= best:
                return
            if cur + len(p) <= best:
                return
            v = p[i]
            np_ = [p[j] for j in range(i) if adj[v][p[j]]]
            ncolors = greedy_color(adj, np_)
            expand(cur + 1, np_, ncolors)
            p.pop(i)

    expand(0, list(range(len(adj))), greedy_color(adj, list(range(len(adj)))))
    return best


# ---------------------------------------------------------------- 演算法三：暴力枚舉子集


def brute_all(adj: List[List[int]]) -> Tuple[int, int, List[int], int]:
    """枚舉 2^n 個子集。回傳 (最大團, 極大團數, 字典序最小的最大團, 最大獨立集)。
    n > BRUTE_MAX 時回傳 (-1, -1, [], -1)。"""
    n = len(adj)
    if n > BRUTE_MAX:
        return (-1, -1, [], -1)
    best = 0
    best_set: List[int] = []
    cnt_maximal = 0
    best_ind = 0

    for mask in range(1 << n):
        verts = [i for i in range(n) if (mask >> i) & 1]
        is_clique = True
        for a in range(len(verts)):
            for b in range(a + 1, len(verts)):
                if not adj[verts[a]][verts[b]]:
                    is_clique = False
                    break
            if not is_clique:
                break
        if is_clique:
            # 字典序最小的最大團（子集枚舉的順序不是字典序，所以要顯式比較）
            if len(verts) > best:
                best = len(verts)
                best_set = verts
            elif len(verts) == best and verts < best_set:
                best_set = verts
            # 是否極大：沒有任何外面的點可以加進來
            maximal = True
            for u in range(n):
                if (mask >> u) & 1:
                    continue
                ok = True
                for w in verts:
                    if not adj[u][w]:
                        ok = False
                        break
                if ok:
                    maximal = False
                    break
            if maximal:
                cnt_maximal += 1
        is_ind = True
        for a in range(len(verts)):
            for b in range(a + 1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    is_ind = False
                    break
            if not is_ind:
                break
        if is_ind and len(verts) > best_ind:
            best_ind = len(verts)
    return (best, cnt_maximal, best_set, best_ind)


# ---------------------------------------------------------------- 主流程


def run(n: int, edges: List[Tuple[int, int]]) -> str:
    adj = build_adj(n, edges)
    ok = True

    mc1, cnt1, best_set = bron_kerbosch(adj)
    mc2 = max_clique_color(adj)
    comp = complement(adj)
    mis = max_clique_color(comp)
    mvc = n - mis

    bmc, bcnt, bset, bmis = brute_all(adj)
    if bmc != -1:
        if mc1 != bmc or mc2 != bmc or cnt1 != bcnt:
            ok = False
        if sorted(best_set) != sorted(bset):
            ok = False
        if mis != bmis:
            ok = False

    if mc1 != mc2:
        ok = False
    if not (0 <= mis <= n) or not (0 <= mvc <= n):
        ok = False
    if mis + mvc != n:
        ok = False
    # 獨立集與點覆蓋的互補關係：補圖的最大團必須等於原圖的最大獨立集
    if bron_kerbosch(comp)[0] != mis:
        ok = False

    lines = []
    lines.append(f"{mc1} {mc2} {bmc}")
    lines.append(f"{cnt1} {bcnt}")
    line3 = [str(len(best_set))]
    line3.extend(str(v) for v in best_set)
    lines.append(" ".join(line3))
    lines.append(f"{mis} {mvc} {bmis}")
    lines.append("1" if ok else "0")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 空圖 ----
    out = run(0, [])
    assert out == "0 0 0\n1 1\n0\n0 0 0\n1\n", out

    # ---- 單點 ----
    out = run(1, [])
    assert out == "1 1 1\n1 1\n1 0\n1 0 1\n1\n", out

    # ---- 三角形 ----
    out = run(3, [(0, 1), (1, 2), (0, 2)])
    assert out == "3 3 3\n1 1\n3 0 1 2\n1 2 1\n1\n", out

    # ---- 路徑 0-1-2（最大團 2）----
    out = run(3, [(0, 1), (1, 2)])
    assert out == "2 2 2\n2 2\n2 0 1\n2 1 2\n1\n", out

    # ---- 三點全空（獨立集 3）----
    out = run(3, [])
    assert out == "1 1 1\n3 3\n1 0\n3 0 3\n1\n", out

    # ---- 四點：0-1-2 三角形 + 2-3 邊 ----
    out = run(4, [(0, 1), (1, 2), (0, 2), (2, 3)])
    assert out.split("\n")[0] == "3 3 3", out
    assert out.split("\n")[-2] == "1", out

    # ---- 隨機對拍 ----
    rng = random.Random(20261019)
    for _ in range(400):
        n = rng.randint(0, 14)
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < rng.choice([0.2, 0.5, 0.8]):
                    edges.append((i, j))
        out = run(n, edges)
        assert out.split("\n")[-2] == "1", (n, edges, out)

    # ---- 較大的圖（不跑暴力）----
    for _ in range(120):
        n = rng.randint(15, MAXN)
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < 0.3:
                    edges.append((i, j))
        out = run(n, edges)
        assert out.split("\n")[-2] == "1", (n, out)

    # ---- 完全圖 Kn ----
    for n in range(0, 9):
        edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
        out = run(n, edges)
        assert out.split("\n")[0] == f"{max(n, 0)} {max(n, 0)} {max(n, 0)}", (n, out)
        assert out.split("\n")[1] == "1 1", (n, out)


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    toks = raw.split()
    idx = 0

    def next_int(default: int) -> int:
        nonlocal idx
        if idx >= len(toks):
            return default
        t = toks[idx]
        idx += 1
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-(10 ** 18), min(10 ** 18, v))

    n = next_int(0)
    if n < 0:
        n = 0
    if n > MAXN:
        n = MAXN
    m = next_int(0)
    if m < 0:
        m = 0
    if m > MAXM:
        m = MAXM
    edges: List[Tuple[int, int]] = []
    for _ in range(m):
        u = next_int(0)
        v = next_int(0)
        if 0 <= u < n and 0 <= v < n and u != v:
            edges.append((u, v))
    emit(run(n, edges))


def main() -> None:
    sys.setrecursionlimit(1000000)
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
