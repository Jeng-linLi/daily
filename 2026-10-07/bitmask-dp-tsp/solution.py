"""狀態壓縮 DP（TSP 旅行商 / 最短哈密頓路徑 / 集合劃分 / SOS 子集和 DP）

題意：
    給定 n 個點（n ≤ 14）與一張**有向帶權圖**的距離矩陣 `dist`，以及一個長度 2^n 的成本陣列 `cost`，
    用狀態壓縮動態規劃解決四個問題：
      1. **TSP 最短迴路**：從 0 號點出發，訪問每個點恰好一次，最後回到 0 號點的最小總代價與具體路徑；
      2. **最短哈密頓路徑**：從 0 號點出發訪問全部點，但**不要求**回到 0 號點；
      3. **集合劃分**：把 {0..n-1} 劃分成若干非空子集，總成本 = Σ cost[子集]，求最小總成本與具體劃分；
      4. **SOS DP（Sum over Subsets）**：對每個 mask 求 f[mask] = Σ_{sub ⊆ mask} cost[sub]。

思路：
    ### 狀態壓縮的核心
    當 n 不大（≤ 20 左右）時，把「一個子集」編碼成一個整數 mask（第 i 位為 1 表示 i 在集合內），
    於是以「子集」為狀態的 DP 就變成一維陣列下標，位運算代替集合運算：
      - 加入元素 i：`mask | (1 << i)`
      - 移除元素 i：`mask & ~(1 << i)`
      - 枚舉 mask 的所有子集：`for (sub = mask; ; sub = (sub - 1) & mask)`（含空集，O(2^popcount)）

    ### TSP：dp[mask][u]
    `dp[mask][u]` = 從 0 出發、已訪問集合恰為 mask、當前停在 u 的最小代價。
    轉移：枚舉下一個未訪問的點 v，
        dp[mask | 1<<v][v] = min(dp[mask | 1<<v][v], dp[mask][u] + dist[u][v])
    共 2^n · n 個狀態、每狀態 n 次轉移 → **O(2^n · n^2)** 時間、**O(2^n · n)** 空間。
    相比枚舉所有排列的 O(n!)，n = 14 時從 8.7e10 降到約 3.2e6，這就是狀態壓縮的威力。
    答案 = `min_v dp[ALL][v] + dist[v][0]`（補上回到起點那一跳）。
    路徑用 `parent[mask][v]` 記錄前驅，回溯即可得到。

    ### 最短哈密頓路徑
    同一張 dp 表，只是答案不補最後一跳：`min_v dp[ALL][v]`。

    ### 集合劃分：dp[mask]
    `dp[mask]` = 把 mask 這個集合劃分成若干非空子集的最小總成本。
    枚舉包含 mask 最低位的那個子集 s（這樣每個劃分只被算一次，避免重複計數）：
        dp[mask] = min_{s ⊆ mask, s 含 lowestbit(mask)} ( dp[mask ^ s] + cost[s] )
    枚舉量是 3^n 的子集和級別 → **O(3^n)**，比枚舉所有劃分（Bell 數，B(14) ≈ 1.9e8）小得多。

    ### SOS DP
    `f[mask] = Σ_{sub ⊆ mask} a[sub]` 若暴力枚舉子集是 O(3^n)；
    SOS DP 按位處理：對每一位 i，把「去掉第 i 位」的結果累加進來：
        for i in 0..n-1: for mask: if mask 的第 i 位為 1: f[mask] += f[mask ^ (1<<i)]
    → **O(n · 2^n)**。常用於「帶位掩碼的數位 DP / 子集計數 / 相容性統計」。

    ### 平手規則（保證 Python 與 C++ 輸出逐字節一致）
      - 狀態枚舉順序、v 的枚舉順序皆為升序；
      - 轉移只在**嚴格更優**時更新（先來先佔）；
      - 集合劃分用「降序枚舉子集 + 嚴格更優才更新」。

應用場景：
    物流配送路線規劃、晶片鑽孔 / 焊接路徑優化（TSP 的直接應用）、
    任務排程與分批（集合劃分）、競賽中的位掩碼計數（SOS DP）、基因組組裝的 overlap 圖。

複雜度：
    TSP / 哈密頓路徑   O(2^n · n^2) 時間、O(2^n · n) 空間
    集合劃分           O(3^n) 時間、O(2^n) 空間
    SOS DP             O(n · 2^n) 時間、O(2^n) 空間

輸入格式（stdin，全部以空白分隔）：
    n
    n × n 個整數：dist[i][j]（i 行 j 列；可為負，但不要求對稱）
    2^n 個整數：cost[0] .. cost[2^n - 1]
    讀不到那麼多時，缺的部分補 0；遇到非整數 token 視為輸入結束。
    為避免狀態數爆炸，n 會被截斷到 14。
輸出格式（stdout）：
    第 1 行：TSP 最短迴路代價
    第 2 行：TSP 路徑（0-indexed，空格分隔，首尾皆為 0）
    第 3 行：最短哈密頓路徑代價
    第 4 行：哈密頓路徑（起點為 0）
    第 5 行：集合劃分最小成本
    第 6 行：劃分組數
    第 7 行起：每組一行（組內元素升序，組間按組內最小元素升序），共「組數」行
    最後一行：SOS DP 的 f[mask]（0 ≤ mask < 2^n，空格分隔）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from itertools import permutations
from typing import List, Optional, Tuple

INF = 10 ** 18
MAX_N = 14                      # 狀態數上限：2^14 · 14 ≈ 2.3e5，避免 IO 模式下爆炸
_INT_RE = re.compile(r"^[+-]?[0-9]+$")   # 嚴格整數規則，與 C++ 的 tryLL 完全一致


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；非十進制整數則返回 None（視為輸入到此為止）。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


def popcount(x: int) -> int:
    """二進制中 1 的個數。"""
    return bin(x).count("1")


def bits_of(mask: int) -> List[int]:
    """把 mask 拆成升序的元素列表。"""
    out = []
    i = 0
    while mask:
        if mask & 1:
            out.append(i)
        mask >>= 1
        i += 1
    return out


# ---------------------------------------------------------------- TSP


def tsp(dist: List[List[int]], n: int) -> Tuple[int, List[int]]:
    """TSP 最短迴路：從 0 出發訪問全部點再回到 0，返回 (代價, 路徑)。"""
    if n == 0:
        return 0, []
    size = 1 << n
    dp = [INF] * (size * n)
    par = [-1] * (size * n)
    dp[0 * n + 0] = 0                                   # dp[1][0]，mask=1 即下標 1
    dp[1 * n + 0] = 0
    for mask in range(size):
        if not (mask & 1):
            continue
        base = mask * n
        for u in range(n):
            if dp[base + u] == INF:
                continue
            cur = dp[base + u]
            row = dist[u]
            for v in range(n):
                if mask >> v & 1:
                    continue
                nmask = mask | (1 << v)
                cand = cur + row[v]
                idx = nmask * n + v
                if cand < dp[idx]:
                    dp[idx] = cand
                    par[idx] = u
    full = size - 1
    best_v, best_c = -1, INF
    for v in range(n):
        if dp[full * n + v] == INF:
            continue
        c = dp[full * n + v] + dist[v][0]
        if c < best_c:                                  # 平手取較小的 v
            best_c, best_v = c, v
    if best_v < 0:
        return 0, []
    path = []
    mask, u = full, best_v
    while mask != -1 and u != -1:
        path.append(u)
        p = par[mask * n + u]
        if p == -1:
            break
        mask ^= (1 << u)
        u = p
    path.reverse()
    path.append(0)                                      # 補上回到起點
    return best_c, path


def tsp_bruteforce(dist: List[List[int]], n: int) -> int:
    """TSP 暴力版（枚舉排列），n ≤ 8，作為對拍基準。"""
    if n == 0:
        return 0
    if n == 1:
        return dist[0][0]
    best = INF
    for perm in permutations(range(1, n)):
        c = dist[0][perm[0]]
        for i in range(len(perm) - 1):
            c += dist[perm[i]][perm[i + 1]]
        c += dist[perm[-1]][0]
        if c < best:
            best = c
    return best


def hamiltonian_path(dist: List[List[int]], n: int) -> Tuple[int, List[int]]:
    """最短哈密頓路徑：從 0 出發訪問全部點，不必回到 0，返回 (代價, 路徑)。"""
    if n == 0:
        return 0, []
    size = 1 << n
    dp = [INF] * (size * n)
    par = [-1] * (size * n)
    dp[1 * n + 0] = 0
    for mask in range(size):
        if not (mask & 1):
            continue
        base = mask * n
        for u in range(n):
            if dp[base + u] == INF:
                continue
            cur = dp[base + u]
            row = dist[u]
            for v in range(n):
                if mask >> v & 1:
                    continue
                nmask = mask | (1 << v)
                cand = cur + row[v]
                idx = nmask * n + v
                if cand < dp[idx]:
                    dp[idx] = cand
                    par[idx] = u
    full = size - 1
    best_v, best_c = -1, INF
    for v in range(n):
        if dp[full * n + v] == INF:
            continue
        if dp[full * n + v] < best_c:
            best_c, best_v = dp[full * n + v], v
    if best_v < 0:
        return 0, []
    path = []
    mask, u = full, best_v
    while True:
        path.append(u)
        p = par[mask * n + u]
        if p == -1:
            break
        mask ^= (1 << u)
        u = p
    path.reverse()
    return best_c, path


def hamiltonian_bruteforce(dist: List[List[int]], n: int) -> int:
    """哈密頓路徑暴力版（枚舉排列），n ≤ 8。"""
    if n == 0:
        return 0
    if n == 1:
        return 0
    best = INF
    for perm in permutations(range(1, n)):
        c = dist[0][perm[0]]
        for i in range(len(perm) - 1):
            c += dist[perm[i]][perm[i + 1]]
        if c < best:
            best = c
    return best


# ---------------------------------------------------------------- 集合劃分


def min_partition_cost(cost: List[int], n: int) -> Tuple[int, List[List[int]]]:
    """集合劃分最小成本：返回 (最小總成本, 劃分方案)。"""
    size = 1 << n
    dp = [INF] * size
    choice = [0] * size
    dp[0] = 0
    for mask in range(1, size):
        low = mask & -mask
        rest = mask ^ low
        sub = rest
        while True:
            s = sub | low
            cand = dp[mask ^ s] + cost[s]
            if cand < dp[mask]:                         # 嚴格更優才更新 → 確定性
                dp[mask] = cand
                choice[mask] = s
            if sub == 0:
                break
            sub = (sub - 1) & rest
    groups: List[List[int]] = []
    mask = size - 1
    while mask:
        s = choice[mask]
        groups.append(bits_of(s))
        mask ^= s
    groups.sort(key=lambda g: g[0])
    return dp[size - 1], groups


def partition_bruteforce(cost: List[int], n: int) -> int:
    """集合劃分暴力版：枚舉所有集合劃分（n ≤ 8），Bell(8) = 4140。"""
    if n == 0:
        return 0
    best = [INF]

    def rec(i: int, groups: List[int]) -> None:
        if i == n:
            total = sum(cost[g] for g in groups)
            if total < best[0]:
                best[0] = total
            return
        for k in range(len(groups)):
            groups[k] |= (1 << i)
            rec(i + 1, groups)
            groups[k] ^= (1 << i)
        groups.append(1 << i)
        rec(i + 1, groups)
        groups.pop()

    rec(0, [])
    return best[0]


# ---------------------------------------------------------------- SOS DP


def sos_dp(a: List[int], n: int) -> List[int]:
    """SOS DP：f[mask] = Σ_{sub ⊆ mask} a[sub]，O(n · 2^n)。"""
    size = 1 << n
    f = a[:size]
    while len(f) < size:
        f.append(0)
    for i in range(n):
        bit = 1 << i
        for mask in range(size):
            if mask & bit:
                f[mask] += f[mask ^ bit]
    return f


def sos_bruteforce(a: List[int], n: int) -> List[int]:
    """SOS 暴力版：直接枚舉子集求和，O(3^n)，n ≤ 10。"""
    size = 1 << n
    out = []
    for mask in range(size):
        total = 0
        sub = mask
        while True:
            total += a[sub] if sub < len(a) else 0
            if sub == 0:
                break
            sub = (sub - 1) & mask
        out.append(total)
    return out


# ---------------------------------------------------------------- IO 與測試


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> int:
        nonlocal pos
        v = parse_int(toks[pos]) if pos < len(toks) else None
        pos += 1
        return 0 if v is None else v

    n = nxt()
    if n < 0:
        n = 0
    if n > MAX_N:
        n = MAX_N
    dist = [[nxt() for _ in range(n)] for _ in range(n)]
    size = 1 << n
    cost = [nxt() for _ in range(size)]

    tc, tpath = tsp(dist, n)
    hc, hpath = hamiltonian_path(dist, n)
    pc, groups = min_partition_cost(cost, n)
    f = sos_dp(cost, n)

    out = [
        str(tc),
        " ".join(str(x) for x in tpath),
        str(hc),
        " ".join(str(x) for x in hpath),
        str(pc),
        str(len(groups)),
    ]
    for g in groups:
        out.append(" ".join(str(x) for x in g))
    out.append(" ".join(str(x) for x in f))
    sys.stdout.write("\n".join(out) + "\n")


def run_tests() -> None:
    assert popcount(0) == 0
    assert popcount(0b101101) == 4
    assert bits_of(0b10110) == [1, 2, 4]

    # ---- SOS DP 對拍 ----
    a = [1, 2, 3, 4]
    assert sos_dp(a, 2) == [1, 3, 4, 10]
    assert sos_dp(a, 2) == sos_bruteforce(a, 2)

    # ---- 固定用例：對稱 TSP ----
    d4 = [
        [0, 10, 15, 20],
        [10, 0, 35, 25],
        [15, 35, 0, 30],
        [20, 25, 30, 0],
    ]
    assert tsp(d4, 4)[0] == 80                       # 0→1→3→2→0 = 10+25+30+15
    assert tsp(d4, 4)[0] == tsp_bruteforce(d4, 4)
    assert hamiltonian_path(d4, 4)[0] == 65          # 0→1→3→2 = 10+25+30
    assert hamiltonian_path(d4, 4)[0] == hamiltonian_bruteforce(d4, 4)

    tp = tsp(d4, 4)[1]
    assert tp[0] == 0 and tp[-1] == 0 and sorted(tp[1:-1]) == [1, 2, 3]
    hp = hamiltonian_path(d4, 4)[1]
    assert hp[0] == 0 and sorted(hp) == [0, 1, 2, 3]

    # ---- 邊界 ----
    assert tsp([], 0) == (0, [])
    assert hamiltonian_path([], 0) == (0, [])
    assert min_partition_cost([0], 0) == (0, [])
    assert sos_dp([7], 0) == [7]
    assert tsp([[5]], 1)[0] == 5                     # 只有一個點：0 → 0
    assert tsp([[5]], 1)[1] == [0, 0]
    assert hamiltonian_path([[5]], 1) == (0, [0])

    # ---- 隨機對拍：TSP / 哈密頓路徑 vs 排列枚舉（n ≤ 8）----
    random.seed(20261007)
    for _ in range(120):
        n = random.randint(1, 7)
        dist = [[random.randint(0, 30) for _ in range(n)] for _ in range(n)]
        tc, tp = tsp(dist, n)
        assert tc == tsp_bruteforce(dist, n)
        assert tp[0] == 0 and tp[-1] == 0 and sorted(tp[1:-1]) == list(range(1, n))
        hc, hp = hamiltonian_path(dist, n)
        assert hc == hamiltonian_bruteforce(dist, n)
        assert hp[0] == 0 and sorted(hp) == list(range(n))

    # ---- 隨機對拍：含負權邊 ----
    for _ in range(120):
        n = random.randint(1, 7)
        dist = [[random.randint(-20, 20) for _ in range(n)] for _ in range(n)]
        assert tsp(dist, n)[0] == tsp_bruteforce(dist, n)
        assert hamiltonian_path(dist, n)[0] == hamiltonian_bruteforce(dist, n)

    # ---- 隨機對拍：集合劃分 vs 枚舉所有劃分（n ≤ 7）----
    for _ in range(120):
        n = random.randint(1, 6)
        size = 1 << n
        cost = [random.randint(-10, 30) for _ in range(size)]
        pc, groups = min_partition_cost(cost, n)
        assert pc == partition_bruteforce(cost, n)
        covered = sorted(x for g in groups for x in g)
        assert covered == list(range(n))             # 劃分必須恰好覆蓋每個元素一次
        assert sum(cost[sum(1 << x for x in g)] for g in groups) == pc

    # ---- 隨機對拍：SOS DP vs 暴力（n ≤ 8）----
    for _ in range(120):
        n = random.randint(0, 8)
        size = 1 << n
        a = [random.randint(-10, 10) for _ in range(size)]
        assert sos_dp(a, n) == sos_bruteforce(a, n)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
