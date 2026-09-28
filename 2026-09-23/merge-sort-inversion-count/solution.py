"""歸併排序與逆序對計數（Merge Sort & Inversion Count）

題意：給定長度爲 n 的整數序列 a，求其中「逆序對」的個數，即滿足
    i < j 且 a[i] > a[j] 的二元組 (i, j) 的數量，並輸出升序排序後的序列。

思路：
    暴力做法是雙重循環 O(n^2)。歸併排序之所以能順便數出逆序對，是因爲
    逆序對天生就是「分治三分類」的：對區間 [lo, hi) 以 mid 切分後，
    任何一個逆序對 (i, j) 恰好屬於下面三類之一，且不重不漏：

      1. i, j 都在左半邊  -> 遞歸統計
      2. i, j 都在右半邊  -> 遞歸統計
      3. i 在左半邊、j 在右半邊（跨中線的逆序對）

    關鍵在第 3 類：合併兩個「已各自有序」的子數組時，若 arr[i] > arr[j]，
    由於左半邊 arr[i..mid) 是升序，arr[i] 後面的元素全都 >= arr[i] > arr[j]，
    於是 arr[i], arr[i+1], ..., arr[mid-1] 與 arr[j] 一次性構成 (mid - i) 個
    逆序對。也就是說，一次比較就能批量結算一整段，這正是把 O(n^2) 降到
    O(n log n) 的原因。

    另兩個容易寫錯的點：
      - 比較必須寫成 arr[i] <= arr[j] 才走左半邊（取等號），否則相等元素會被
        誤判成逆序對，破壞「逆序對 = 嚴格大於」的定義。
      - 歸併排序是穩定排序，是否穩定就取決於這個等號的方向。

    補充：逆序對個數恰好等於「只允許交換相鄰元素」時把序列排好序所需的最少
    交換次數（冒泡排序的交換次數）；完全逆序的序列逆序對數爲 n*(n-1)/2，
    是本問題的上界，所以計數變量要用 64 位整數（C++ 側用 long long）。

輸入格式（stdin）：
    第一行：n
    第二行：n 個整數（可跨行書寫）
輸出格式（stdout）：
    第一行：逆序對個數
    第二行：升序排序後的序列（空格分隔；n = 0 時輸出空行）
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from typing import List, Tuple


def _merge_sort_count(arr: List[int], buf: List[int], lo: int, hi: int) -> int:
    """對 arr[lo:hi) 歸併排序，返回其中的逆序對個數（原地寫回 arr）。"""
    if hi - lo <= 1:
        return 0

    mid = (lo + hi) // 2
    # 左右兩半內部的逆序對各自遞歸統計
    inv = _merge_sort_count(arr, buf, lo, mid)
    inv += _merge_sort_count(arr, buf, mid, hi)

    i, j, k = lo, mid, lo
    while i < mid and j < hi:
        if arr[i] <= arr[j]:
            # 取等號：相等元素不構成逆序對，且保證排序穩定
            buf[k] = arr[i]
            i += 1
        else:
            buf[k] = arr[j]
            j += 1
            # 左半邊 arr[i..mid) 全部 > arr[j]，一次性結算 mid - i 個逆序對
            inv += mid - i
        k += 1

    while i < mid:
        buf[k] = arr[i]
        i += 1
        k += 1
    while j < hi:
        buf[k] = arr[j]
        j += 1
        k += 1

    arr[lo:hi] = buf[lo:hi]
    return inv


def sort_and_count(nums: List[int]) -> Tuple[int, List[int]]:
    """返回 (逆序對個數, 升序排序後的新列表)。時間 O(n log n)，空間 O(n)。"""
    arr = list(nums)
    buf = [0] * len(arr)
    inv = _merge_sort_count(arr, buf, 0, len(arr))
    return inv, arr


def count_inversions_brute(nums: List[int]) -> int:
    """對照用的 O(n^2) 暴力枚舉，僅用於小規模測試驗證。"""
    n = len(nums)
    return sum(1 for i in range(n) for j in range(i + 1, n) if nums[i] > nums[j])


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
    tokens = data.split()
    if not tokens:
        return
    n = int(tokens[0])
    nums = [int(x) for x in tokens[1 : 1 + n]]
    inv, sorted_nums = sort_and_count(nums)
    print(inv)
    print(" ".join(str(x) for x in sorted_nums))


def run_tests() -> None:
    # README 中的示例：2 3 8 6 1 -> 逆序對 5 個
    #   (2,1) (3,1) (8,6) (8,1) (6,1)
    inv, sorted_nums = sort_and_count([2, 3, 8, 6, 1])
    assert inv == 5
    assert sorted_nums == [1, 2, 3, 6, 8]
    assert count_inversions_brute([2, 3, 8, 6, 1]) == 5

    # 空序列與單元素
    assert sort_and_count([]) == (0, [])
    assert sort_and_count([42]) == (0, [42])

    # 已升序：0 個逆序對
    assert sort_and_count([1, 2, 3, 4, 5]) == (0, [1, 2, 3, 4, 5])

    # 完全逆序：n*(n-1)/2 個逆序對，驗證 64 位計數不溢出
    assert sort_and_count([5, 4, 3, 2, 1]) == (10, [1, 2, 3, 4, 5])
    assert sort_and_count(list(range(2000, 0, -1)))[0] == 2000 * 1999 // 2

    # 相等元素不算逆序對（這裡最容易把 <= 寫成 < 而數多）
    assert sort_and_count([2, 2, 1]) == (2, [1, 2, 2])
    assert sort_and_count([1, 1, 1]) == (0, [1, 1, 1])
    assert sort_and_count([3, 1, 3, 1]) == (3, [1, 1, 3, 3])

    # 負數與零：逆序對爲 (-1,-3) (-1,-2) (0,-2) (2,-2)，共 4 個
    assert sort_and_count([-1, -3, 0, 2, -2]) == (4, [-3, -2, -1, 0, 2])

    # 與暴力解隨機對拍：同時校驗「逆序對數一致」與「結果確實有序且是原序列的排列」
    import random

    random.seed(20260923)
    for _ in range(300):
        n = random.randint(0, 40)
        nums = [random.randint(-20, 20) for _ in range(n)]
        inv2, sorted2 = sort_and_count(nums)
        assert inv2 == count_inversions_brute(nums)      # 與暴力枚舉一致
        assert sorted2 == sorted(nums)                   # 排序結果正確
        assert len(sorted2) == n                         # 元素一個不多一個不少

    # 排序不應改動調用方傳入的原列表
    original = [3, 1, 2]
    sort_and_count(original)
    assert original == [3, 1, 2]

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
