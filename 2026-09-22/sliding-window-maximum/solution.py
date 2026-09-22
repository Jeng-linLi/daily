"""滑动窗口最大值（单调队列 / Monotonic Queue）

题意：给定数组 nums 和窗口大小 k，窗口从左向右每次滑动一格，
    求每个窗口内的最大值，共 n-k+1 个结果。

思路：
    朴素做法是每个窗口扫一遍取最大，O(n*k)，k 大时会超时。
    单调队列把它优化到 O(n)：
      队列中保存的是"下标"，且对应的值严格递减（队首永远是当前窗口最大值）。
      1) 入队前，从队尾弹出所有 <= 当前值的元素——它们既比当前值小、
         又比当前值早出窗口，永远不可能成为答案，可以安全丢弃；
      2) 入队当前下标；
      3) 从队首弹出所有已滑出窗口的下标（下标 <= i-k）；
      4) 当 i >= k-1 时，队首下标对应的值就是当前窗口最大值。
    每个元素恰好入队一次、出队一次，所以总时间线性。

输入格式（stdin）：
    第一行 n k
    第二行 n 个整数
输出格式（stdout）：
    一行 n-k+1 个整数，空格分隔，为各窗口最大值
无 stdin 输入时运行内置断言测试。
"""

import sys
from collections import deque
from typing import List


def max_sliding_window(nums: List[int], k: int) -> List[int]:
    """返回长度为 k 的滑动窗口在每个位置的最大值，时间 O(n)，空间 O(k)。"""
    if not nums or k <= 0:
        return []
    if k == 1:
        return list(nums)
    if k >= len(nums):
        return [max(nums)]

    q: deque = deque()  # 存下标，保证 nums[q[0]] > nums[q[1]] > ...
    ans: List[int] = []

    for i, val in enumerate(nums):
        # 1) 队尾所有不大于当前值的下标都不可能再成为答案
        while q and nums[q[-1]] <= val:
            q.pop()
        # 2) 当前下标入队
        q.append(i)
        # 3) 队首已滑出窗口的下标出队
        while q and q[0] <= i - k:
            q.popleft()
        # 4) 窗口成型后，队首即最大值
        if i >= k - 1:
            ans.append(nums[q[0]])

    return ans


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    tokens = data.split()
    if not tokens:
        return
    n, k = int(tokens[0]), int(tokens[1])
    nums = [int(x) for x in tokens[2:2 + n]]
    print(" ".join(str(x) for x in max_sliding_window(nums, k)))


def brute_force(nums: List[int], k: int) -> List[int]:
    """对照用的朴素实现，O(n*k)，仅用于测试验证。"""
    return [max(nums[i:i + k]) for i in range(len(nums) - k + 1)]


def run_tests() -> None:
    assert max_sliding_window([1, 3, -1, -3, 5, 3, 6, 7], 3) == [3, 3, 5, 5, 6, 7]
    assert max_sliding_window([1], 1) == [1]
    assert max_sliding_window([1, -1], 1) == [1, -1]
    assert max_sliding_window([9, 8, 7, 6, 5], 3) == [9, 8, 7]   # 递减：队首不断被挤出
    assert max_sliding_window([1, 2, 3, 4, 5], 3) == [3, 4, 5]   # 递增：队尾不断被弹出
    assert max_sliding_window([5, 5, 5, 5], 2) == [5, 5, 5]      # 全相等
    assert max_sliding_window([-7, -8, -7, -6, -5], 3) == [-7, -6, -5]
    assert max_sliding_window([1, 3, 1, 2, 0, 5], 3) == [3, 3, 2, 5]
    assert max_sliding_window([4, 3, 2, 1], 4) == [4]            # k == n
    assert max_sliding_window([4, 3, 2, 1], 5) == [4]            # k > n
    assert max_sliding_window([], 3) == []                       # 空数组

    # 与朴素实现随机对照：确保单调队列没有边界错误
    import random
    random.seed(20260922)
    for _ in range(200):
        n = random.randint(1, 40)
        k = random.randint(1, n)
        arr = [random.randint(-50, 50) for _ in range(n)]
        assert max_sliding_window(arr, k) == brute_force(arr, k)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
