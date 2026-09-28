"""二分查找與二分答案（lower_bound / upper_bound / 精確查找 / 最小化最大值）

題意：
    1) 給定**升序**數組 a（允許重複、允許爲空），回答 q 次查詢：
       對每個 x 輸出 lower_bound —— 第一個 >= x 的下標，
       以及 upper_bound —— 第一個 > x 的下標。
       下標爲 0-based 的「插入位置」，取值範圍是 [0, n]，等於 n 表示不存在。
    2) 附帶演示「二分答案」：把數組按原順序切成 k 個非空連續段，
       最小化「最大段的和」（LeetCode 410 Split Array Largest Sum）。
       這部分只在內置斷言裡驗證，不參與 stdin IO。

思路：
    二分查找的本質是**在一個單調的判定函數上找分界點**。把候選位置想成
    [False, False, ..., True, True] 這樣的一段，用左閉右開區間 [lo, hi)
    維護「還沒確定」的部分，每次取中點 mid 把區間砍一半：
      - lower_bound：判定 P(i) = (a[i] >= x)，找第一個使 P 爲真的 i；
      - upper_bound：判定 P(i) = (a[i] >  x)，找第一個使 P 爲真的 i。
    寫成 `if a[mid] < x: lo = mid + 1 else: hi = mid` 的形式，
    mid 永遠落在 [lo, hi) 內，既不會死循環也不會越界。

    這兩個函數隻差一個比較符號（`<` 與 `<=`），是最容易寫錯的地方。
    有了它們之後：
      - 精確查找 = lower_bound 取出位置後判等，不等就是 -1；
      - 等於 x 的元素個數 = upper_bound(x) - lower_bound(x)。

    二分答案思路完全一樣，只是判定函數換成「給定上限 limit，能否把數組
    切成不超過 k 段且每段和都 <= limit」。可行性關於 limit 單調（越大越容易），
    所以能二分出最小的可行 limit。下界取 max(a)（任何一段至少要裝下最大元素），
    上界取 sum(a)（全部塞進一段一定可行）。
    判定用貪心：從左到右累加，一旦超過 limit 就在當前位置切一刀，
    段數超過 k 則判定不可行。貪心的正確性在於「能裝就裝」不會讓後面變差。
    （要求數組元素非負，這樣可行性對 limit 才是單調的。）

輸入格式（stdin，所有數字按空白分隔即可）：
    n
    a1 a2 ... an        （升序；n = 0 時這一行直接省略）
    q
    x1 x2 ... xq
輸出格式（stdout）：
    每個查詢一行：`<lower> <upper>`
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import sys
from typing import Iterator, List


def lower_bound(a: List[int], x: int) -> int:
    """第一個 >= x 的下標；不存在則返回 len(a)。時間 O(log n)，空間 O(1)。"""
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:      # a[mid] 太小，答案在右半邊
            lo = mid + 1
        else:               # a[mid] >= x，mid 本身可能是答案
            hi = mid
    return lo


def upper_bound(a: List[int], x: int) -> int:
    """第一個 > x 的下標；不存在則返回 len(a)。時間 O(log n)，空間 O(1)。"""
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] <= x:     # 與 lower_bound 唯一的差別：相等也算「太小」
            lo = mid + 1
        else:
            hi = mid
    return lo


def binary_search(a: List[int], x: int) -> int:
    """精確查找：返回第一個等於 x 的下標，不存在返回 -1。"""
    i = lower_bound(a, x)
    return i if i < len(a) and a[i] == x else -1


def count_equal(a: List[int], x: int) -> int:
    """等於 x 的元素個數：upper_bound 減 lower_bound。"""
    return upper_bound(a, x) - lower_bound(a, x)


# ---------------- 二分答案 ----------------

def can_split(a: List[int], k: int, limit: int) -> bool:
    """貪心判定：能否切成不超過 k 段且每段和 <= limit。時間 O(n)。"""
    parts, cur = 1, 0
    for v in a:
        if v > limit:
            return False      # 單個元素就超過 limit，這一段無論如何裝不下
        if cur + v <= limit:
            cur += v          # 還能裝下，繼續往當前段裏塞
        else:
            parts += 1        # 裝不下了，在這裡切一刀
            cur = v
            if parts > k:
                return False
    return True


def min_max_split(a: List[int], k: int) -> int:
    """把 a 按原順序切成 k 個非空連續段，最小化最大段的和。

    要求元素非負。時間 O(n log S)，S = sum(a) - max(a)；空間 O(1)。
    """
    if not a:
        return 0
    k = min(k, len(a))        # 段數超過元素個數沒有意義
    lo, hi = max(a), sum(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if can_split(a, k, mid):
            hi = mid          # mid 可行，答案 <= mid
        else:
            lo = mid + 1      # mid 不可行，答案 > mid
    return lo


def min_max_split_brute(a: List[int], k: int) -> int:
    """對照用的 O(n^2 * k) 動態規劃，僅用於小規模測試驗證。

    dp[i][p] = 把前 i 個元素切成 p 段時的最小「最大段和」。
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
            for j in range(p - 1, i):          # 最後一段是 (j, i]
                if dp[j][p - 1] == INF:
                    continue
                val = max(dp[j][p - 1], pre[i] - pre[j])
                if val < best:
                    best = val
            dp[i][p] = best
    return dp[n][k]


# ---------------- 對照用的暴力實現 ----------------

def lower_bound_brute(a: List[int], x: int) -> int:
    """線性掃描版 lower_bound，用於對拍。"""
    for i, v in enumerate(a):
        if v >= x:
            return i
    return len(a)


def upper_bound_brute(a: List[int], x: int) -> int:
    """線性掃描版 upper_bound，用於對拍。"""
    for i, v in enumerate(a):
        if v > x:
            return i
    return len(a)


# ---------------- IO ----------------

def _next_int(it: Iterator[str], default: int = 0) -> int:
    """取下一個整數；輸入被截斷時用默認值兜底，避免直接拋異常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
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
    assert lower_bound(a, 2) == 1          # 第一個 >= 2 的位置
    assert upper_bound(a, 2) == 4          # 第一個 > 2 的位置
    assert count_equal(a, 2) == 3          # 一共 3 個 2
    assert binary_search(a, 2) == 1
    assert binary_search(a, 3) == -1       # 3 不存在
    assert count_equal(a, 3) == 0

    # 邊界：比所有元素都小 / 都大
    assert lower_bound(a, 0) == 0
    assert upper_bound(a, 0) == 0
    assert lower_bound(a, 9) == 6
    assert upper_bound(a, 9) == 6
    assert count_equal(a, 9) == 0

    # 落在兩個元素之間的空隙裏
    assert lower_bound(a, 3) == 4
    assert upper_bound(a, 3) == 4
    assert lower_bound(a, 5) == 5
    assert upper_bound(a, 5) == 5

    # 命中唯一元素 / 命中最大元素
    assert lower_bound(a, 1) == 0
    assert upper_bound(a, 1) == 1
    assert lower_bound(a, 7) == 5
    assert upper_bound(a, 7) == 6

    # 空數組：任何查詢都返回 0
    assert lower_bound([], 5) == 0
    assert upper_bound([], 5) == 0
    assert binary_search([], 5) == -1
    assert count_equal([], 5) == 0

    # 單元素數組
    assert lower_bound([5], 5) == 0
    assert upper_bound([5], 5) == 1
    assert binary_search([5], 4) == -1

    # 二分答案：LeetCode 410 的經典用例
    assert min_max_split([7, 2, 5, 10, 8], 2) == 18      # [7,2,5] | [10,8]
    assert min_max_split([1, 2, 3, 4, 5], 2) == 9        # [1,2,3] | [4,5]
    assert min_max_split([1, 4, 4], 3) == 4              # 每段一個元素
    assert min_max_split([7, 2, 5, 10, 8], 5) == 10      # 段數 >= n 時就是 max(a)
    assert min_max_split([7, 2, 5, 10, 8], 9) == 10
    assert min_max_split([], 3) == 0                     # 空數組
    assert min_max_split([5], 1) == 5

    import random

    random.seed(20260927)

    # 隨機對拍一：lower_bound / upper_bound 與線性掃描逐一比對
    for _ in range(300):
        n = random.randint(0, 12)
        arr = sorted(random.randint(0, 10) for _ in range(n))   # 故意製造重複元素
        for x in range(-2, 13):
            assert lower_bound(arr, x) == lower_bound_brute(arr, x)
            assert upper_bound(arr, x) == upper_bound_brute(arr, x)
            assert 0 <= lower_bound(arr, x) <= upper_bound(arr, x) <= n
            # 精確查找與計數要和 lower/upper 自洽
            if binary_search(arr, x) == -1:
                assert count_equal(arr, x) == 0
                assert lower_bound(arr, x) == upper_bound(arr, x)
            else:
                assert arr[binary_search(arr, x)] == x
                assert binary_search(arr, x) == lower_bound(arr, x)
                assert count_equal(arr, x) == upper_bound(arr, x) - lower_bound(arr, x)

    # 隨機對拍二：二分答案與 O(n^2*k) 的 DP 暴力解比對
    for _ in range(200):
        n = random.randint(0, 8)
        arr = [random.randint(0, 20) for _ in range(n)]         # 不要求有序，非負即可
        for k in range(1, n + 1):
            got = min_max_split(arr, k)
            exp = min_max_split_brute(arr, k)
            assert got == exp, (arr, k, got, exp)
            assert can_split(arr, k, got) is True                # 最優值確實可行
            if got > 0:
                assert can_split(arr, k, got - 1) is False       # 再小一點就不行了
        if n == 0:
            assert min_max_split(arr, 1) == 0

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
