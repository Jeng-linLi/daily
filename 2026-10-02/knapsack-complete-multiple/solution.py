"""背包變體：多重背包（二進制拆分）/ 完全背包 / 湊齊容量的最少件數

題意：
    給定 n 種物品（重量 w、價值 v、數量上限 c）與背包容量 C，依次回答三問：
      1. **多重背包**：每種物品最多取 c_i 個，容量 C 內的最大價值，並輸出一組達到最優的選取個數；
      2. **完全背包**：每種物品可以無限取，容量 C 內的最大價值，並輸出一組達到最優的選取個數；
      3. **湊齊容量的最少件數**：每種物品無限取，湊出「恰好總重 C」最少要幾件（-1 表示無解）。

思路：
    0-1 背包的 1D 寫法是「容量**倒序**遍歷」，目的是讓 `dp[cap - w]` 讀到的是**上一輪**（還沒放過本物品）的值，
    從而保證每件物品只被用一次。把遍歷方向反過來，就能得到兩種變體：

    1. **完全背包**：容量**正序**遍歷。`dp[cap - w]` 讀到的是**本輪已更新過**的值，
       相當於「本物品可以被反覆加入」，一次正序掃描就等價於枚舉取 0,1,2,… 個。O(n·C)。
    2. **多重背包**：每種物品有數量上限，既不能倒序（那會變成 0-1），也不能正序（那會變成無限）。
       標準做法是**二進制拆分**：把 c 拆成 1, 2, 4, …, 剩餘 的若干「捆」，
       每一捆當成一個獨立的 0-1 物品（重量 k·w、價值 k·v）。
       由於 1..c 的任意整數都能由這組 2 的冪唯一表示，這些捆的 0-1 組合恰好覆蓋「取 0..c 個」的所有情況，
       複雜度從 O(n·C·c) 降到 O(n·C·log c)。
    3. **最少件數湊齊容量**：把 max 換成 min、初值換成 INF，`dp[0] = 0`，
       `dp[cap] = min(dp[cap], dp[cap − w] + 1)`，容量正序（無限取）。
       最後 `dp[C]` 若仍是 INF 就是無解，輸出 −1。

    還原選取方案：
      - 多重背包用**二維 DP**（`dp[i][cap]` 表示只考慮前 i 種、容量不超過 cap 的最優值），
        從 `i = n` 往回走，對每種物品枚舉取了幾個 t，找第一個滿足
        `dp[i][cap] == dp[i-1][cap − t·w] + t·v` 的 t 即可（取最小的 t，保證輸出唯一）。
      - 完全背包用 1D DP 加上 `par[cap]` 記錄「這個容量最後是被哪種物品更新的」，
        然後從 cap = C 一路往回減。`par` 回溯對 0-1 / 多重**不成立**（會把同一捆重複使用而超出數量上限），
        只對完全背包這種「可重複使用」的情形成立。

    約定：重量必須 ≥ 1；`w = 0` 的物品直接忽略（否則價值無界，問題無意義）。

輸入格式（stdin）：
    n C
    w1 v1 c1
    w2 v2 c2
    …… （共 n 行）
輸出格式（stdout）：
    第 1 行：多重背包的最大價值
    第 2 行：多重背包的選取個數（n 個非負整數，空格分隔；n = 0 時輸出空行）
    第 3 行：完全背包的最大價值
    第 4 行：完全背包的選取個數（n 個非負整數）
    第 5 行：湊出恰好總重 C 的最少件數（無解輸出 -1）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import List, Tuple


# ---------------------------------------------------------------- 多重背包
def multiple_knapsack_2d(nf: int, C: int, w: List[int], v: List[int],
                         c: List[int]) -> Tuple[int, List[int]]:
    """多重背包：二維 DP（O(n·C·c)），回傳 (最大價值, 每種物品取幾個)。

    dp[i][cap] = 只考慮前 i 種物品、總重不超過 cap 時的最大價值。
    轉移：dp[i][cap] = max over t in [0, c_i] of dp[i-1][cap − t·w_i] + t·v_i
    """
    dp = [[0] * (C + 1) for _ in range(nf + 1)]
    for i in range(1, nf + 1):
        wi, vi, ci = w[i - 1], v[i - 1], c[i - 1]
        prev = dp[i - 1]
        cur = dp[i]
        for cap in range(C + 1):
            best = prev[cap]                      # 取 0 個
            t = 1
            while t <= ci and t * wi <= cap:
                cand = prev[cap - t * wi] + t * vi
                if cand > best:
                    best = cand
                t += 1
            cur[cap] = best

    # 回溯：從第 nf 種往回，取滿足轉移式的最小 t
    take = [0] * nf
    cap = C
    for i in range(nf, 0, -1):
        wi, vi, ci = w[i - 1], v[i - 1], c[i - 1]
        chosen = 0
        for t in range(0, min(ci, cap // wi) + 1):
            if dp[i][cap] == dp[i - 1][cap - t * wi] + t * vi:
                chosen = t
                break
        take[i - 1] = chosen
        cap -= chosen * wi
    return dp[nf][C], take


def multiple_knapsack_split(nf: int, C: int, w: List[int], v: List[int],
                            c: List[int]) -> int:
    """多重背包：二進制拆分 + 0-1 背包（O(n·C·log c)），只求最大價值（測試裡與 2D 版對拍）。"""
    dp = [0] * (C + 1)
    for i in range(nf):
        wi, vi = w[i], v[i]
        rest = c[i]
        k = 1
        while rest > 0:
            take = min(k, rest)                   # 拆成 1, 2, 4, …, 剩餘
            bw, bv = take * wi, take * vi
            for cap in range(C, bw - 1, -1):      # 0-1 背包：容量倒序
                if dp[cap - bw] + bv > dp[cap]:
                    dp[cap] = dp[cap - bw] + bv
            rest -= take
            k <<= 1
    return dp[C]


# ---------------------------------------------------------------- 完全背包
def complete_knapsack(nf: int, C: int, w: List[int], v: List[int]) -> Tuple[int, List[int]]:
    """完全背包：1D DP 容量正序（O(n·C)），回傳 (最大價值, 每種物品取幾個)。"""
    dp = [0] * (C + 1)
    par = [-1] * (C + 1)                          # par[cap] = 該容量最後被哪種物品更新
    for i in range(nf):
        wi, vi = w[i], v[i]
        for cap in range(wi, C + 1):              # 正序：允許同一種物品被重複取用
            cand = dp[cap - wi] + vi
            if cand > dp[cap]:
                dp[cap] = cand
                par[cap] = i
    take = [0] * nf
    cap = C
    while cap > 0:
        i = par[cap]
        if i < 0:
            break
        take[i] += 1
        cap -= w[i]
    return dp[C], take


def min_items_exact(nf: int, C: int, w: List[int]) -> int:
    """每種物品無限取，湊出恰好總重 C 的最少件數；無解回傳 -1。"""
    INF = 10 ** 9
    dp = [INF] * (C + 1)
    dp[0] = 0
    for i in range(nf):
        wi = w[i]
        for cap in range(wi, C + 1):
            if dp[cap - wi] + 1 < dp[cap]:
                dp[cap] = dp[cap - wi] + 1
    return dp[C] if dp[C] < INF else -1


# ---------------------------------------------------------------- 暴力基準（只用於測試）
def brute_best(nf: int, C: int, w: List[int], v: List[int], c: List[int]) -> int:
    """枚舉所有 (t1, …, tn) 組合求最大價值，只在測試裡當基準。"""
    best = 0

    def dfs(i: int, weight: int, value: int) -> None:
        nonlocal best
        if i == nf:
            if weight <= C and value > best:
                best = value
            return
        wi, vi, ci = w[i], v[i], c[i]
        t = 0
        while t <= ci and weight + t * wi <= C:
            dfs(i + 1, weight + t * wi, value + t * vi)
            t += 1

    dfs(0, 0, 0)
    return best


def brute_min_items(nf: int, C: int, w: List[int]) -> int:
    """Bellman-Ford 式鬆弛求最少件數，只在測試裡當基準。"""
    INF = 10 ** 9
    dp = [INF] * (C + 1)
    dp[0] = 0
    for _ in range(C):
        changed = False
        for cap in range(1, C + 1):
            for i in range(nf):
                if cap - w[i] >= 0 and dp[cap - w[i]] + 1 < dp[cap]:
                    dp[cap] = dp[cap - w[i]] + 1
                    changed = True
        if not changed:
            break
    return dp[C] if dp[C] < INF else -1


# ---------------------------------------------------------------- IO 模式
def run_io(data: str) -> None:
    tokens = data.split()
    if not tokens:
        return
    pos = 0

    def nxt() -> int:
        nonlocal pos
        if pos < len(tokens):
            val = int(tokens[pos])
            pos += 1
            return val
        return 0

    n = nxt()
    C = nxt()
    w_all: List[int] = []
    v_all: List[int] = []
    c_all: List[int] = []
    for _ in range(max(0, n)):
        w_all.append(nxt())
        v_all.append(nxt())
        c_all.append(nxt())

    n = max(0, n)
    C = max(0, C)
    idx = [i for i in range(n) if w_all[i] >= 1]      # 忽略 w = 0 的物品（價值無界）
    w = [w_all[i] for i in idx]
    v = [v_all[i] for i in idx]
    c = [max(0, c_all[i]) for i in idx]
    nf = len(idx)

    best_multi, take_multi = multiple_knapsack_2d(nf, C, w, v, c)
    cnt_multi = [0] * n
    for j, i in enumerate(idx):
        cnt_multi[i] = take_multi[j]
    print(best_multi)
    print(" ".join(map(str, cnt_multi)))

    best_full, take_full = complete_knapsack(nf, C, w, v)
    cnt_full = [0] * n
    for j, i in enumerate(idx):
        cnt_full[i] = take_full[j]
    print(best_full)
    print(" ".join(map(str, cnt_full)))

    print(min_items_exact(nf, C, w))


# ---------------------------------------------------------------- 測試
def run_tests() -> None:
    # README 示例：(2,3,2) (3,4,1) (4,5,3)，C = 10
    w = [2, 3, 4]
    v = [3, 4, 5]
    c = [2, 1, 3]
    C = 10
    best, take = multiple_knapsack_2d(3, C, w, v, c)
    assert best == 13                                  # 1 個 w=2 + 2 個 w=4 → 重 10、值 13
    assert take == [1, 0, 2]
    assert best == multiple_knapsack_split(3, C, w, v, c)
    assert best == brute_best(3, C, w, v, c)
    b2, t2 = complete_knapsack(3, C, w, v)
    assert b2 == 15                                    # 5 個 w=2 → 重 10、值 15
    assert sum(t2[i] * v[i] for i in range(3)) == b2
    assert sum(t2[i] * w[i] for i in range(3)) <= C
    assert min_items_exact(3, C, w) == 3               # 4 + 4 + 2 = 10，三件

    # 空輸入
    assert multiple_knapsack_2d(0, 0, [], [], []) == (0, [])
    assert multiple_knapsack_split(0, 0, [], [], []) == 0
    assert complete_knapsack(0, 0, [], []) == (0, [])
    assert min_items_exact(0, 0, []) == 0
    assert min_items_exact(1, 5, [3]) == -1            # 3 湊不出 5

    # 容量 0：價值 0、件數 0
    assert multiple_knapsack_2d(2, 0, [1, 2], [5, 9], [3, 3]) == (0, [0, 0])
    assert complete_knapsack(2, 0, [1, 2], [5, 9]) == (0, [0, 0])
    assert min_items_exact(2, 0, [1, 2]) == 0

    # 數量上限為 0 時該物品完全不能用
    assert multiple_knapsack_2d(2, 5, [1, 2], [5, 9], [0, 3])[0] == 18   # 只能用第二種：2 個 w=2

    # 裝不下任何東西
    assert multiple_knapsack_2d(1, 3, [5], [100], [1]) == (0, [0])
    assert complete_knapsack(1, 3, [5], [100]) == (0, [0])

    # 隨機對拍：2D 版 vs 二進制拆分版 vs 暴力枚舉
    random.seed(20261002)
    for _ in range(400):
        nf = random.randint(0, 4)
        C = random.randint(0, 12)
        w = [random.randint(1, 5) for _ in range(nf)]
        v = [random.randint(0, 9) for _ in range(nf)]
        c = [random.randint(0, 3) for _ in range(nf)]

        best, take = multiple_knapsack_2d(nf, C, w, v, c)
        assert best == multiple_knapsack_split(nf, C, w, v, c)
        assert best == brute_best(nf, C, w, v, c)
        # 回溯出的方案必須真的可行且剛好達到最優值
        assert all(0 <= take[i] <= c[i] for i in range(nf))
        assert sum(take[i] * w[i] for i in range(nf)) <= C
        assert sum(take[i] * v[i] for i in range(nf)) == best

        # 完全背包 = 把數量上限設成「最多能裝幾個」之後的多重背包
        cap_c = [C // w[i] for i in range(nf)]
        bf, tf = complete_knapsack(nf, C, w, v)
        assert bf == brute_best(nf, C, w, v, cap_c)
        assert sum(tf[i] * w[i] for i in range(nf)) <= C
        assert sum(tf[i] * v[i] for i in range(nf)) == bf

        assert min_items_exact(nf, C, w) == brute_min_items(nf, C, w)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
