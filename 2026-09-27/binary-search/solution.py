"""二分查找与二分答案（lower_bound / upper_bound / 精确查找 / 最小化最大值）

题意：
    1) 给定**升序**数组 a（允许重复、允许为空），回答 q 次查询：
       对每个 x 输出 lower_bound —— 第一个 >= x 的下标，
       以及 upper_bound —— 第一个 > x 的下标。
       下标为 0-based 的「插入位置」，取值范围是 [0, n]，等于 n 表示不存在。
    2) 附带演示「二分答案」：把数组按原顺序切成 k 个非空连续段，
       最小化「最大段的和」（LeetCode 410 Split Array Largest Sum）。
       这部分只在内置断言里验证，不参与 stdin IO。

思路：
    二分查找的本质是**在一个单调的判定函数上找分界点**。把候选位置想成
    [False, False, ..., True, True] 这样的一段，用左闭右开区间 [lo, hi)
    维护「还没确定」的部分，每次取中点 mid 把区间砍一半：
      - lower_bound：判定 P(i) = (a[i] >= x)，找第一个使 P 为真的 i；
      - upper_bound：判定 P(i) = (a[i] >  x)，找第一个使 P 为真的 i。
    写成 `if a[mid] < x: lo = mid + 1 else: hi = mid` 的形式，
    mid 永远落在 [lo, hi) 内，既不会死循环也不会越界。

    这两个函数只差一个比较符号（`<` 与 `<=`），是最容易写错的地方。
    有了它们之后：
      - 精确查找 = lower_bound 取出位置后判等，不等就是 -1；
      - 等于 x 的元素个数 = upper_bound(x) - lower_bound(x)。

    二分答案思路完全一样，只是判定函数换成「给定上限 limit，能否把数组
    切成不超过 k 段且每段和都 <= limit」。可行性关于 limit 单调（越大越容易），
    所以能二分出最小的可行 limit。下界取 max(a)（任何一段至少要装下最大元素），
    上界取 sum(a)（全部塞进一段一定可行）。
    判定用贪心：从左到右累加，一旦超过 limit 就在当前位置切一刀，
    段数超过 k 则判定不可行。贪心的正确性在于「能装就装」不会让后面变差。
    （要求数组元素非负，这样可行性对 limit 才是单调的。）

输入格式（stdin，所有数字按空白分隔即可）：
    n
    a1 a2 ... an        （升序；n = 0 时这一行直接省略）
    q
    x1 x2 ... xq
输出格式（stdout）：
    每个查询一行：`<lower> <upper>`
无 stdin 输入时运行内置断言测试并输出 `all tests passed`。
"""

import sys
from typing import Iterator, List


def lower_bound(a: List[int], x: int) -> int:
    """第一个 >= x 的下标；不存在则返回 len(a)。时间 O(log n)，空间 O(1)。"""
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:      # a[mid] 太小，答案在右半边
            lo = mid + 1
        else:               # a[mid] >= x，mid 本身可能是答案
            hi = mid
    return lo


def upper_bound(a: List[int], x: int) -> int:
    """第一个 > x 的下标；不存在则返回 len(a)。时间 O(log n)，空间 O(1)。"""
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] <= x:     # 与 lower_bound 唯一的差别：相等也算「太小」
            lo = mid + 1
        else:
            hi = mid
    return lo


def binary_search(a: List[int], x: int) -> int:
    """精确查找：返回第一个等于 x 的下标，不存在返回 -1。"""
    i = lower_bound(a, x)
    return i if i < len(a) and a[i] == x else -1


def count_equal(a: List[int], x: int) -> int:
    """等于 x 的元素个数：upper_bound 减 lower_bound。"""
    return upper_bound(a, x) - lower_bound(a, x)


# ---------------- 二分答案 ----------------

def can_split(a: List[int], k: int, limit: int) -> bool:
    """贪心判定：能否切成不超过 k 段且每段和 <= limit。时间 O(n)。"""
    parts, cur = 1, 0
    for v in a:
        if v > limit:
            return False      # 单个元素就超过 limit，这一段无论如何装不下
        if cur + v <= limit:
            cur += v          # 还能装下，继续往当前段里塞
        else:
            parts += 1        # 装不下了，在这里切一刀
            cur = v
            if parts > k:
                return False
    return True


def min_max_split(a: List[int], k: int) -> int:
    """把 a 按原顺序切成 k 个非空连续段，最小化最大段的和。

    要求元素非负。时间 O(n log S)，S = sum(a) - max(a)；空间 O(1)。
    """
    if not a:
        return 0
    k = min(k, len(a))        # 段数超过元素个数没有意义
    lo, hi = max(a), sum(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if can_split(a, k, mid):
            hi = mid          # mid 可行，答案 <= mid
        else:
            lo = mid + 1      # mid 不可行，答案 > mid
    return lo


def min_max_split_brute(a: List[int], k: int) -> int:
    """对照用的 O(n^2 * k) 动态规划，仅用于小规模测试验证。

    dp[i][p] = 把前 i 个元素切成 p 段时的最小「最大段和」。
    """
    n = len(a)
    if n == 0:
        return 0
    k = min(k, n)
    pre = [0] * (n + 1)
    for i, v in enumerate(a):
        pre[i + 1] = pre[i] + v
    INF = float("inf")
    dp = [[INF] * (k + 1) for _ in range(n + 1)]
    dp[0][0] = 0
    for i in range(1, n + 1):
        for p in range(1, min(k, i) + 1):
            best = INF
            for j in range(p - 1, i):          # 最后一段是 (j, i]
                if dp[j][p - 1] == INF:
                    continue
                val = max(dp[j][p - 1], pre[i] - pre[j])
                if val < best:
                    best = val
            dp[i][p] = best
    return dp[n][k]


# ---------------- 对照用的暴力实现 ----------------

def lower_bound_brute(a: List[int], x: int) -> int:
    """线性扫描版 lower_bound，用于对拍。"""
    for i, v in enumerate(a):
        if v >= x:
            return i
    return len(a)


def upper_bound_brute(a: List[int], x: int) -> int:
    """线性扫描版 upper_bound，用于对拍。"""
    for i, v in enumerate(a):
        if v > x:
            return i
    return len(a)


# ---------------- IO ----------------

def _next_int(it: Iterator[str], default: int = 0) -> int:
    """取下一个整数；输入被截断时用默认值兜底，避免直接抛异常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    it = iter(data.split())
    n = _next_int(it)
    a = [_next_int(it) for _ in range(n)]
    q = _next_int(it)
    for _ in range(q):
        x = _next_int(it)
        print(f"{lower_bound(a, x)} {upper_bound(a, x)}")


def run_tests() -> None:
    # README 示例：a = [1, 2, 2, 2, 4, 7]
    a = [1, 2, 2, 2, 4, 7]
    assert lower_bound(a, 2) == 1          # 第一个 >= 2 的位置
    assert upper_bound(a, 2) == 4          # 第一个 > 2 的位置
    assert count_equal(a, 2) == 3          # 一共 3 个 2
    assert binary_search(a, 2) == 1
    assert binary_search(a, 3) == -1       # 3 不存在
    assert count_equal(a, 3) == 0

    # 边界：比所有元素都小 / 都大
    assert lower_bound(a, 0) == 0
    assert upper_bound(a, 0) == 0
    assert lower_bound(a, 9) == 6
    assert upper_bound(a, 9) == 6
    assert count_equal(a, 9) == 0

    # 落在两个元素之间的空隙里
    assert lower_bound(a, 3) == 4
    assert upper_bound(a, 3) == 4
    assert lower_bound(a, 5) == 5
    assert upper_bound(a, 5) == 5

    # 命中唯一元素 / 命中最大元素
    assert lower_bound(a, 1) == 0
    assert upper_bound(a, 1) == 1
    assert lower_bound(a, 7) == 5
    assert upper_bound(a, 7) == 6

    # 空数组：任何查询都返回 0
    assert lower_bound([], 5) == 0
    assert upper_bound([], 5) == 0
    assert binary_search([], 5) == -1
    assert count_equal([], 5) == 0

    # 单元素数组
    assert lower_bound([5], 5) == 0
    assert upper_bound([5], 5) == 1
    assert binary_search([5], 4) == -1

    # 二分答案：LeetCode 410 的经典用例
    assert min_max_split([7, 2, 5, 10, 8], 2) == 18      # [7,2,5] | [10,8]
    assert min_max_split([1, 2, 3, 4, 5], 2) == 9        # [1,2,3] | [4,5]
    assert min_max_split([1, 4, 4], 3) == 4              # 每段一个元素
    assert min_max_split([7, 2, 5, 10, 8], 5) == 10      # 段数 >= n 时就是 max(a)
    assert min_max_split([7, 2, 5, 10, 8], 9) == 10
    assert min_max_split([], 3) == 0                     # 空数组
    assert min_max_split([5], 1) == 5

    import random

    random.seed(20260927)

    # 随机对拍一：lower_bound / upper_bound 与线性扫描逐一比对
    for _ in range(300):
        n = random.randint(0, 12)
        arr = sorted(random.randint(0, 10) for _ in range(n))   # 故意制造重复元素
        for x in range(-2, 13):
            assert lower_bound(arr, x) == lower_bound_brute(arr, x)
            assert upper_bound(arr, x) == upper_bound_brute(arr, x)
            assert 0 <= lower_bound(arr, x) <= upper_bound(arr, x) <= n
            # 精确查找与计数要和 lower/upper 自洽
            if binary_search(arr, x) == -1:
                assert count_equal(arr, x) == 0
                assert lower_bound(arr, x) == upper_bound(arr, x)
            else:
                assert arr[binary_search(arr, x)] == x
                assert binary_search(arr, x) == lower_bound(arr, x)
                assert count_equal(arr, x) == upper_bound(arr, x) - lower_bound(arr, x)

    # 随机对拍二：二分答案与 O(n^2*k) 的 DP 暴力解比对
    for _ in range(200):
        n = random.randint(0, 8)
        arr = [random.randint(0, 20) for _ in range(n)]         # 不要求有序，非负即可
        for k in range(1, n + 1):
            got = min_max_split(arr, k)
            exp = min_max_split_brute(arr, k)
            assert got == exp, (arr, k, got, exp)
            assert can_split(arr, k, got) is True                # 最优值确实可行
            if got > 0:
                assert can_split(arr, k, got - 1) is False       # 再小一点就不行了
        if n == 0:
            assert min_max_split(arr, 1) == 0

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
