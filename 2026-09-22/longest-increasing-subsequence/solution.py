"""最长递增子序列（Longest Increasing Subsequence, LIS）

题意：给定整数数组 nums，求最长的「严格递增」子序列的长度，并还原出一条具体方案。
子序列不要求元素连续，但相对顺序必须与原数组一致。

思路：
  方法一 · 动态规划 O(n^2)
      dp[i] 表示「以 nums[i] 作为结尾」的最长递增子序列长度。
      dp[i] = 1 + max{ dp[j] | j < i 且 nums[j] < nums[i] }，若不存在这样的 j 则 dp[i] = 1。
      用 pre[i] 记录前驱下标，最后从 dp 最大的位置往回跳即可还原序列。

  方法二 · 贪心 + 二分 O(n log n)
      维护 tails 数组：tails[k] = 「长度为 k+1 的递增子序列」的结尾元素的最小可能值。
      tails 严格递增，所以对 x = nums[i] 二分找到第一个 >= x 的位置 pos：
        - pos == len(tails)：x 比所有结尾都大，可以接长，长度 +1；
        - 否则：用 x 覆盖 tails[pos]。结尾越小、后续接长的潜力越大，这是贪心的关键。
      注意 tails 本身未必是一条合法子序列，只有它的长度 len(tails) 是正确答案。

输入格式（stdin）：
    第一行 n
    第二行 n 个整数
输出格式（stdout）：
    第一行 LIS 长度
    第二行 一条 LIS（空格分隔，n == 0 时输出空行）
无 stdin 输入时运行内置断言测试。
"""

import sys
from bisect import bisect_left
from typing import List


def length_of_lis(nums: List[int]) -> int:
    """贪心 + 二分，O(n log n)，只求长度。"""
    tails: List[int] = []  # tails[k] = 长度 k+1 的递增子序列的最小结尾
    for x in nums:
        pos = bisect_left(tails, x)  # 第一个 >= x 的下标（保证严格递增）
        if pos == len(tails):
            tails.append(x)          # x 能接在所有已有序列后面
        else:
            tails[pos] = x           # 用更小的结尾替换，留出增长空间
    return len(tails)


def lis_dp(nums: List[int]) -> List[int]:
    """动态规划，O(n^2)，返回一条具体的最长递增子序列。"""
    n = len(nums)
    if n == 0:
        return []
    dp = [1] * n        # dp[i]：以 nums[i] 结尾的 LIS 长度
    pre = [-1] * n      # pre[i]：最优前驱下标
    best = 0            # dp 最大值的下标
    for i in range(n):
        for j in range(i):
            if nums[j] < nums[i] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                pre[i] = j
        if dp[i] > dp[best]:
            best = i

    seq: List[int] = []
    k = best
    while k != -1:
        seq.append(nums[k])
        k = pre[k]
    seq.reverse()
    return seq


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    tokens = data.split()
    n = int(tokens[0])
    nums = [int(t) for t in tokens[1:1 + n]]
    out = [str(len(lis_dp(nums)))]
    out.append(" ".join(str(x) for x in lis_dp(nums)))
    print("\n".join(out))


def run_tests() -> None:
    # 经典用例
    assert length_of_lis([10, 9, 2, 5, 3, 7, 101, 18]) == 4
    assert lis_dp([10, 9, 2, 5, 3, 7, 101, 18]) == [2, 5, 7, 101]
    # 相等元素不算递增
    assert length_of_lis([7, 7, 7, 7]) == 1
    assert lis_dp([7, 7, 7, 7]) == [7]
    # 含重复但仍能取更长
    assert length_of_lis([0, 1, 0, 3, 2, 3]) == 4
    assert lis_dp([0, 1, 0, 3, 2, 3]) == [0, 1, 2, 3]
    # 边界
    assert length_of_lis([]) == 0
    assert lis_dp([]) == []
    assert length_of_lis([1]) == 1
    assert lis_dp([1]) == [1]
    # 完全递减
    assert length_of_lis([5, 4, 3, 2, 1]) == 1
    # 完全递增
    assert length_of_lis([1, 2, 3, 4, 5]) == 5
    assert lis_dp([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]
    # 两种方法结果必须一致（随机小数据交叉验证）
    import random
    for _ in range(200):
        arr = [random.randint(-5, 5) for _ in range(random.randint(0, 12))]
        assert length_of_lis(arr) == len(lis_dp(arr))
    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
