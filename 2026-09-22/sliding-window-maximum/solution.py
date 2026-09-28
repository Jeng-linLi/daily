"""滑動窗口最大值（單調隊列 / Monotonic Queue）

題意：給定數組 nums 和窗口大小 k，窗口從左向右每次滑動一格，
    求每個窗口內的最大值，共 n-k+1 個結果。

思路：
    樸素做法是每個窗口掃一遍取最大，O(n*k)，k 大時會超時。
    單調隊列把它優化到 O(n)：
      隊列中保存的是"下標"，且對應的值嚴格遞減（隊首永遠是當前窗口最大值）。
      1) 入隊前，從隊尾彈出所有 <= 當前值的元素——它們既比當前值小、
         又比當前值早出窗口，永遠不可能成爲答案，可以安全丟棄；
      2) 入隊當前下標；
      3) 從隊首彈出所有已滑出窗口的下標（下標 <= i-k）；
      4) 當 i >= k-1 時，隊首下標對應的值就是當前窗口最大值。
    每個元素恰好入隊一次、出隊一次，所以總時間線性。

輸入格式（stdin）：
    第一行 n k
    第二行 n 個整數
輸出格式（stdout）：
    一行 n-k+1 個整數，空格分隔，爲各窗口最大值
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from collections import deque
from typing import List


def max_sliding_window(nums: List[int], k: int) -> List[int]:
    """返回長度爲 k 的滑動窗口在每個位置的最大值，時間 O(n)，空間 O(k)。"""
    if not nums or k <= 0:
        return []
    if k == 1:
        return list(nums)
    if k >= len(nums):
        return [max(nums)]

    q: deque = deque()  # 存下標，保證 nums[q[0]] > nums[q[1]] > ...
    ans: List[int] = []

    for i, val in enumerate(nums):
        # 1) 隊尾所有不大於當前值的下標都不可能再成爲答案
        while q and nums[q[-1]] <= val:
            q.pop()
        # 2) 當前下標入隊
        q.append(i)
        # 3) 隊首已滑出窗口的下標出隊
        while q and q[0] <= i - k:
            q.popleft()
        # 4) 窗口成型後，隊首即最大值
        if i >= k - 1:
            ans.append(nums[q[0]])

    return ans


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
    tokens = data.split()
    if not tokens:
        return
    n, k = int(tokens[0]), int(tokens[1])
    nums = [int(x) for x in tokens[2:2 + n]]
    print(" ".join(str(x) for x in max_sliding_window(nums, k)))


def brute_force(nums: List[int], k: int) -> List[int]:
    """對照用的樸素實現，O(n*k)，僅用於測試驗證。"""
    return [max(nums[i:i + k]) for i in range(len(nums) - k + 1)]


def run_tests() -> None:
    assert max_sliding_window([1, 3, -1, -3, 5, 3, 6, 7], 3) == [3, 3, 5, 5, 6, 7]
    assert max_sliding_window([1], 1) == [1]
    assert max_sliding_window([1, -1], 1) == [1, -1]
    assert max_sliding_window([9, 8, 7, 6, 5], 3) == [9, 8, 7]   # 遞減：隊首不斷被擠出
    assert max_sliding_window([1, 2, 3, 4, 5], 3) == [3, 4, 5]   # 遞增：隊尾不斷被彈出
    assert max_sliding_window([5, 5, 5, 5], 2) == [5, 5, 5]      # 全相等
    assert max_sliding_window([-7, -8, -7, -6, -5], 3) == [-7, -6, -5]
    assert max_sliding_window([1, 3, 1, 2, 0, 5], 3) == [3, 3, 2, 5]
    assert max_sliding_window([4, 3, 2, 1], 4) == [4]            # k == n
    assert max_sliding_window([4, 3, 2, 1], 5) == [4]            # k > n
    assert max_sliding_window([], 3) == []                       # 空數組

    # 與樸素實現隨機對照：確保單調隊列沒有邊界錯誤
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
