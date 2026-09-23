"""归并排序与逆序对计数（Merge Sort & Inversion Count）

题意：给定长度为 n 的整数序列 a，求其中「逆序对」的个数，即满足
    i < j 且 a[i] > a[j] 的二元组 (i, j) 的数量，并输出升序排序后的序列。

思路：
    暴力做法是双重循环 O(n^2)。归并排序之所以能顺便数出逆序对，是因为
    逆序对天生就是「分治三分类」的：对区间 [lo, hi) 以 mid 切分后，
    任何一个逆序对 (i, j) 恰好属于下面三类之一，且不重不漏：

      1. i, j 都在左半边  -> 递归统计
      2. i, j 都在右半边  -> 递归统计
      3. i 在左半边、j 在右半边（跨中线的逆序对）

    关键在第 3 类：合并两个「已各自有序」的子数组时，若 arr[i] > arr[j]，
    由于左半边 arr[i..mid) 是升序，arr[i] 后面的元素全都 >= arr[i] > arr[j]，
    于是 arr[i], arr[i+1], ..., arr[mid-1] 与 arr[j] 一次性构成 (mid - i) 个
    逆序对。也就是说，一次比较就能批量结算一整段，这正是把 O(n^2) 降到
    O(n log n) 的原因。

    另两个容易写错的点：
      - 比较必须写成 arr[i] <= arr[j] 才走左半边（取等号），否则相等元素会被
        误判成逆序对，破坏「逆序对 = 严格大于」的定义。
      - 归并排序是稳定排序，是否稳定就取决于这个等号的方向。

    补充：逆序对个数恰好等于「只允许交换相邻元素」时把序列排好序所需的最少
    交换次数（冒泡排序的交换次数）；完全逆序的序列逆序对数为 n*(n-1)/2，
    是本问题的上界，所以计数变量要用 64 位整数（C++ 侧用 long long）。

输入格式（stdin）：
    第一行：n
    第二行：n 个整数（可跨行书写）
输出格式（stdout）：
    第一行：逆序对个数
    第二行：升序排序后的序列（空格分隔；n = 0 时输出空行）
无 stdin 输入时运行内置断言测试。
"""

import sys
from typing import List, Tuple


def _merge_sort_count(arr: List[int], buf: List[int], lo: int, hi: int) -> int:
    """对 arr[lo:hi) 归并排序，返回其中的逆序对个数（原地写回 arr）。"""
    if hi - lo <= 1:
        return 0

    mid = (lo + hi) // 2
    # 左右两半内部的逆序对各自递归统计
    inv = _merge_sort_count(arr, buf, lo, mid)
    inv += _merge_sort_count(arr, buf, mid, hi)

    i, j, k = lo, mid, lo
    while i < mid and j < hi:
        if arr[i] <= arr[j]:
            # 取等号：相等元素不构成逆序对，且保证排序稳定
            buf[k] = arr[i]
            i += 1
        else:
            buf[k] = arr[j]
            j += 1
            # 左半边 arr[i..mid) 全部 > arr[j]，一次性结算 mid - i 个逆序对
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
    """返回 (逆序对个数, 升序排序后的新列表)。时间 O(n log n)，空间 O(n)。"""
    arr = list(nums)
    buf = [0] * len(arr)
    inv = _merge_sort_count(arr, buf, 0, len(arr))
    return inv, arr


def count_inversions_brute(nums: List[int]) -> int:
    """对照用的 O(n^2) 暴力枚举，仅用于小规模测试验证。"""
    n = len(nums)
    return sum(1 for i in range(n) for j in range(i + 1, n) if nums[i] > nums[j])


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    tokens = data.split()
    if not tokens:
        return
    n = int(tokens[0])
    nums = [int(x) for x in tokens[1 : 1 + n]]
    inv, sorted_nums = sort_and_count(nums)
    print(inv)
    print(" ".join(str(x) for x in sorted_nums))


def run_tests() -> None:
    # README 中的示例：2 3 8 6 1 -> 逆序对 5 个
    #   (2,1) (3,1) (8,6) (8,1) (6,1)
    inv, sorted_nums = sort_and_count([2, 3, 8, 6, 1])
    assert inv == 5
    assert sorted_nums == [1, 2, 3, 6, 8]
    assert count_inversions_brute([2, 3, 8, 6, 1]) == 5

    # 空序列与单元素
    assert sort_and_count([]) == (0, [])
    assert sort_and_count([42]) == (0, [42])

    # 已升序：0 个逆序对
    assert sort_and_count([1, 2, 3, 4, 5]) == (0, [1, 2, 3, 4, 5])

    # 完全逆序：n*(n-1)/2 个逆序对，验证 64 位计数不溢出
    assert sort_and_count([5, 4, 3, 2, 1]) == (10, [1, 2, 3, 4, 5])
    assert sort_and_count(list(range(2000, 0, -1)))[0] == 2000 * 1999 // 2

    # 相等元素不算逆序对（这里最容易把 <= 写成 < 而数多）
    assert sort_and_count([2, 2, 1]) == (2, [1, 2, 2])
    assert sort_and_count([1, 1, 1]) == (0, [1, 1, 1])
    assert sort_and_count([3, 1, 3, 1]) == (3, [1, 1, 3, 3])

    # 负数与零：逆序对为 (-1,-3) (-1,-2) (0,-2) (2,-2)，共 4 个
    assert sort_and_count([-1, -3, 0, 2, -2]) == (4, [-3, -2, -1, 0, 2])

    # 与暴力解随机对拍：同时校验「逆序对数一致」与「结果确实有序且是原序列的排列」
    import random

    random.seed(20260923)
    for _ in range(300):
        n = random.randint(0, 40)
        nums = [random.randint(-20, 20) for _ in range(n)]
        inv2, sorted2 = sort_and_count(nums)
        assert inv2 == count_inversions_brute(nums)      # 与暴力枚举一致
        assert sorted2 == sorted(nums)                   # 排序结果正确
        assert len(sorted2) == n                         # 元素一个不多一个不少

    # 排序不应改动调用方传入的原列表
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
