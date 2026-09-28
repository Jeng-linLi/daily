"""0-1 背包（Knapsack 0-1，動態規劃）

題意：有 n 個物品，第 i 個重量 w[i]、價值 v[i]，背包容量爲 W。
    每個物品只能選或不選（不能切分、不能重複選），求能裝下的最大總價值，
    並給出一種達到該價值的選取方案。

思路：
    定義 dp[i][c] = 只考慮前 i 個物品、容量爲 c 時能取得的最大價值。
    對第 i 個物品只有兩種決策：
      不選：dp[i][c] = dp[i-1][c]
      選它：dp[i][c] = dp[i-1][c-w[i]] + v[i]   （要求 c >= w[i]）
    取兩者最大，這就是最優子結構；物品只能選一次，所以轉移只依賴 i-1 層，
    不會產生「完全背包」那種同層自我疊加。

    空間可壓到一維：dp[c] = 前 i 個物品下的最優值，但**容量 c 必須倒序遍歷**
    （從 W 到 w[i]）。倒序保證 dp[c-w[i]] 讀到的仍是 i-1 層的舊值；
    若正序遍歷，同一個物品會被重複放入，退化成完全背包——這是最經典的坑。

    要還原選取方案則必須保留二維表：從 dp[n][W] 反推，
    若 dp[i][c] != dp[i-1][c] 說明第 i 個物品被選了，跳到 dp[i-1][c-w[i]] 繼續。

輸入格式（stdin）：
    第一行 n W
    接下來 n 行，每行 w_i v_i
輸出格式（stdout）：
    第一行：最大總價值
    第二行：被選中的物品下標（0-based，升序，空格分隔；一個都沒選則輸出空行）
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from typing import List, Tuple


def knapsack_max(weights: List[int], values: List[int], capacity: int) -> int:
    """只求最大價值：一維滾動數組，時間 O(n*W)，空間 O(W)。"""
    dp = [0] * (capacity + 1)
    for w, v in zip(weights, values):
        # 容量必須倒序，正序會讓同一物品被重複選取（變成完全背包）
        for c in range(capacity, w - 1, -1):
            cand = dp[c - w] + v
            if cand > dp[c]:
                dp[c] = cand
    return dp[capacity]


def knapsack_with_items(
    weights: List[int], values: List[int], capacity: int
) -> Tuple[int, List[int]]:
    """求最大價值並還原一種選取方案：保留二維表，空間 O(n*W)。"""
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        w, v = weights[i - 1], values[i - 1]
        for c in range(capacity + 1):
            best = dp[i - 1][c]                       # 不選第 i-1 個物品
            if c >= w:
                take = dp[i - 1][c - w] + v           # 選它
                if take > best:
                    best = take
            dp[i][c] = best

    # 反向還原：dp[i][c] 比 dp[i-1][c] 大，說明第 i-1 個物品被選中了
    chosen: List[int] = []
    c = capacity
    for i in range(n, 0, -1):
        if dp[i][c] != dp[i - 1][c]:
            chosen.append(i - 1)
            c -= weights[i - 1]
    chosen.reverse()
    return dp[n][capacity], chosen


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
    tokens = data.split()
    if not tokens:
        return
    n, W = int(tokens[0]), int(tokens[1])
    weights: List[int] = []
    values: List[int] = []
    idx = 2
    for _ in range(n):
        weights.append(int(tokens[idx]))
        values.append(int(tokens[idx + 1]))
        idx += 2
    _, chosen = knapsack_with_items(weights, values, W)
    total = sum(values[i] for i in chosen)
    print(total)
    print(" ".join(str(i) for i in chosen))


def brute_force(weights: List[int], values: List[int], capacity: int) -> int:
    """對照用的指數級枚舉，僅用於小規模測試驗證。"""
    n = len(weights)
    best = 0
    for mask in range(1 << n):
        tw = tv = 0
        for i in range(n):
            if mask >> i & 1:
                tw += weights[i]
                tv += values[i]
        if tw <= capacity and tv > best:
            best = tv
    return best


def run_tests() -> None:
    # 經典用例：容量 10 的最優值爲 12，但存在多解
    #   (a) 2 號 + 3 號：w=4+6=10, v=5+7=12
    #   (b) 0 號 + 1 號 + 2 號：w=2+3+4=9, v=3+4+5=12
    # 反推時「dp[i][c] == dp[i-1][c] 視爲未選」，會優先得到 (b)，故只斷言值最優、方案合法
    w = [2, 3, 4, 6]
    v = [3, 4, 5, 7]
    best, chosen = knapsack_with_items(w, v, 10)
    assert best == 12
    assert sum(w[i] for i in chosen) <= 10
    assert sum(v[i] for i in chosen) == 12
    assert knapsack_max(w, v, 10) == 12

    # 容量爲 0 / 物品爲空
    assert knapsack_max([], [], 10) == 0
    assert knapsack_max([5], [9], 0) == 0
    assert knapsack_with_items([5], [9], 0) == (0, [])

    # 單件裝不下
    assert knapsack_max([5], [9], 4) == 0

    # 全部都能裝下
    assert knapsack_max([1, 2, 3], [1, 2, 3], 10) == 6

    # 按價值密度貪心會選錯：密度最高的是 0 號(2.0)，但最優解是 1+2 號 w=10 v=18
    assert knapsack_max([5, 4, 6], [10, 7, 11], 10) == 18
    assert knapsack_max([5, 4, 6], [10, 7, 11], 10) == brute_force([5, 4, 6], [10, 7, 11], 10)

    # 零重量物品：價值白拿，且不會造成死循環
    assert knapsack_max([0, 3], [5, 4], 3) == 9

    # 與暴力枚舉隨機對照：同時校驗「最優值一致」與「方案確實可行且達到最優值」
    import random
    random.seed(20260922)
    for _ in range(300):
        n = random.randint(1, 10)
        cap = random.randint(0, 20)
        ws = [random.randint(0, 8) for _ in range(n)]
        vs = [random.randint(0, 20) for _ in range(n)]
        best2, chosen2 = knapsack_with_items(ws, vs, cap)
        assert best2 == brute_force(ws, vs, cap)      # 最優值與暴力一致
        assert knapsack_max(ws, vs, cap) == best2     # 一維版與二維版一致
        assert sum(ws[i] for i in chosen2) <= cap     # 方案不超容量
        assert sum(vs[i] for i in chosen2) == best2   # 方案確實達到最優值

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
