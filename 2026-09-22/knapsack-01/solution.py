"""0-1 背包（Knapsack 0-1，动态规划）

题意：有 n 个物品，第 i 个重量 w[i]、价值 v[i]，背包容量为 W。
    每个物品只能选或不选（不能切分、不能重复选），求能装下的最大总价值，
    并给出一种达到该价值的选取方案。

思路：
    定义 dp[i][c] = 只考虑前 i 个物品、容量为 c 时能取得的最大价值。
    对第 i 个物品只有两种决策：
      不选：dp[i][c] = dp[i-1][c]
      选它：dp[i][c] = dp[i-1][c-w[i]] + v[i]   （要求 c >= w[i]）
    取两者最大，这就是最优子结构；物品只能选一次，所以转移只依赖 i-1 层，
    不会产生「完全背包」那种同层自我叠加。

    空间可压到一维：dp[c] = 前 i 个物品下的最优值，但**容量 c 必须倒序遍历**
    （从 W 到 w[i]）。倒序保证 dp[c-w[i]] 读到的仍是 i-1 层的旧值；
    若正序遍历，同一个物品会被重复放入，退化成完全背包——这是最经典的坑。

    要还原选取方案则必须保留二维表：从 dp[n][W] 反推，
    若 dp[i][c] != dp[i-1][c] 说明第 i 个物品被选了，跳到 dp[i-1][c-w[i]] 继续。

输入格式（stdin）：
    第一行 n W
    接下来 n 行，每行 w_i v_i
输出格式（stdout）：
    第一行：最大总价值
    第二行：被选中的物品下标（0-based，升序，空格分隔；一个都没选则输出空行）
无 stdin 输入时运行内置断言测试。
"""

import sys
from typing import List, Tuple


def knapsack_max(weights: List[int], values: List[int], capacity: int) -> int:
    """只求最大价值：一维滚动数组，时间 O(n*W)，空间 O(W)。"""
    dp = [0] * (capacity + 1)
    for w, v in zip(weights, values):
        # 容量必须倒序，正序会让同一物品被重复选取（变成完全背包）
        for c in range(capacity, w - 1, -1):
            cand = dp[c - w] + v
            if cand > dp[c]:
                dp[c] = cand
    return dp[capacity]


def knapsack_with_items(
    weights: List[int], values: List[int], capacity: int
) -> Tuple[int, List[int]]:
    """求最大价值并还原一种选取方案：保留二维表，空间 O(n*W)。"""
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        w, v = weights[i - 1], values[i - 1]
        for c in range(capacity + 1):
            best = dp[i - 1][c]                       # 不选第 i-1 个物品
            if c >= w:
                take = dp[i - 1][c - w] + v           # 选它
                if take > best:
                    best = take
            dp[i][c] = best

    # 反向还原：dp[i][c] 比 dp[i-1][c] 大，说明第 i-1 个物品被选中了
    chosen: List[int] = []
    c = capacity
    for i in range(n, 0, -1):
        if dp[i][c] != dp[i - 1][c]:
            chosen.append(i - 1)
            c -= weights[i - 1]
    chosen.reverse()
    return dp[n][capacity], chosen


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
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
    """对照用的指数级枚举，仅用于小规模测试验证。"""
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
    # 经典用例：容量 10 的最优值为 12，但存在多解
    #   (a) 2 号 + 3 号：w=4+6=10, v=5+7=12
    #   (b) 0 号 + 1 号 + 2 号：w=2+3+4=9, v=3+4+5=12
    # 反推时「dp[i][c] == dp[i-1][c] 视为未选」，会优先得到 (b)，故只断言值最优、方案合法
    w = [2, 3, 4, 6]
    v = [3, 4, 5, 7]
    best, chosen = knapsack_with_items(w, v, 10)
    assert best == 12
    assert sum(w[i] for i in chosen) <= 10
    assert sum(v[i] for i in chosen) == 12
    assert knapsack_max(w, v, 10) == 12

    # 容量为 0 / 物品为空
    assert knapsack_max([], [], 10) == 0
    assert knapsack_max([5], [9], 0) == 0
    assert knapsack_with_items([5], [9], 0) == (0, [])

    # 单件装不下
    assert knapsack_max([5], [9], 4) == 0

    # 全部都能装下
    assert knapsack_max([1, 2, 3], [1, 2, 3], 10) == 6

    # 按价值密度贪心会选错：密度最高的是 0 号(2.0)，但最优解是 1+2 号 w=10 v=18
    assert knapsack_max([5, 4, 6], [10, 7, 11], 10) == 18
    assert knapsack_max([5, 4, 6], [10, 7, 11], 10) == brute_force([5, 4, 6], [10, 7, 11], 10)

    # 零重量物品：价值白拿，且不会造成死循环
    assert knapsack_max([0, 3], [5, 4], 3) == 9

    # 与暴力枚举随机对照：同时校验「最优值一致」与「方案确实可行且达到最优值」
    import random
    random.seed(20260922)
    for _ in range(300):
        n = random.randint(1, 10)
        cap = random.randint(0, 20)
        ws = [random.randint(0, 8) for _ in range(n)]
        vs = [random.randint(0, 20) for _ in range(n)]
        best2, chosen2 = knapsack_with_items(ws, vs, cap)
        assert best2 == brute_force(ws, vs, cap)      # 最优值与暴力一致
        assert knapsack_max(ws, vs, cap) == best2     # 一维版与二维版一致
        assert sum(ws[i] for i in chosen2) <= cap     # 方案不超容量
        assert sum(vs[i] for i in chosen2) == best2   # 方案确实达到最优值

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
