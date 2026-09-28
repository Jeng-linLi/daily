"""最長遞增子序列（Longest Increasing Subsequence, LIS）

題意：給定整數數組 nums，求最長的「嚴格遞增」子序列的長度，並還原出一條具體方案。
子序列不要求元素連續，但相對順序必須與原數組一致。

思路：
  方法一 · 動態規劃 O(n^2)
      dp[i] 表示「以 nums[i] 作爲結尾」的最長遞增子序列長度。
      dp[i] = 1 + max{ dp[j] | j < i 且 nums[j] < nums[i] }，若不存在這樣的 j 則 dp[i] = 1。
      用 pre[i] 記錄前驅下標，最後從 dp 最大的位置往回跳即可還原序列。

  方法二 · 貪心 + 二分 O(n log n)
      維護 tails 數組：tails[k] = 「長度爲 k+1 的遞增子序列」的結尾元素的最小可能值。
      tails 嚴格遞增，所以對 x = nums[i] 二分找到第一個 >= x 的位置 pos：
        - pos == len(tails)：x 比所有結尾都大，可以接長，長度 +1；
        - 否則：用 x 覆蓋 tails[pos]。結尾越小、後續接長的潛力越大，這是貪心的關鍵。
      注意 tails 本身未必是一條合法子序列，只有它的長度 len(tails) 是正確答案。

輸入格式（stdin）：
    第一行 n
    第二行 n 個整數
輸出格式（stdout）：
    第一行 LIS 長度
    第二行 一條 LIS（空格分隔，n == 0 時輸出空行）
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from bisect import bisect_left
from typing import List


def length_of_lis(nums: List[int]) -> int:
    """貪心 + 二分，O(n log n)，只求長度。"""
    tails: List[int] = []  # tails[k] = 長度 k+1 的遞增子序列的最小結尾
    for x in nums:
        pos = bisect_left(tails, x)  # 第一個 >= x 的下標（保證嚴格遞增）
        if pos == len(tails):
            tails.append(x)          # x 能接在所有已有序列後面
        else:
            tails[pos] = x           # 用更小的結尾替換，留出增長空間
    return len(tails)


def lis_dp(nums: List[int]) -> List[int]:
    """動態規劃，O(n^2)，返回一條具體的最長遞增子序列。"""
    n = len(nums)
    if n == 0:
        return []
    dp = [1] * n        # dp[i]：以 nums[i] 結尾的 LIS 長度
    pre = [-1] * n      # pre[i]：最優前驅下標
    best = 0            # dp 最大值的下標
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
    """按統一輸入輸出格式處理 stdin 數據。"""
    tokens = data.split()
    n = int(tokens[0])
    nums = [int(t) for t in tokens[1:1 + n]]
    out = [str(len(lis_dp(nums)))]
    out.append(" ".join(str(x) for x in lis_dp(nums)))
    print("\n".join(out))


def run_tests() -> None:
    # 經典用例
    assert length_of_lis([10, 9, 2, 5, 3, 7, 101, 18]) == 4
    assert lis_dp([10, 9, 2, 5, 3, 7, 101, 18]) == [2, 5, 7, 101]
    # 相等元素不算遞增
    assert length_of_lis([7, 7, 7, 7]) == 1
    assert lis_dp([7, 7, 7, 7]) == [7]
    # 含重複但仍能取更長
    assert length_of_lis([0, 1, 0, 3, 2, 3]) == 4
    assert lis_dp([0, 1, 0, 3, 2, 3]) == [0, 1, 2, 3]
    # 邊界
    assert length_of_lis([]) == 0
    assert lis_dp([]) == []
    assert length_of_lis([1]) == 1
    assert lis_dp([1]) == [1]
    # 完全遞減
    assert length_of_lis([5, 4, 3, 2, 1]) == 1
    # 完全遞增
    assert length_of_lis([1, 2, 3, 4, 5]) == 5
    assert lis_dp([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]
    # 兩種方法結果必須一致（隨機小數據交叉驗證）
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
