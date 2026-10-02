"""差分陣列與前綴和：一維區間加 / 二維子矩陣加 / 子矩陣求和 / 和為 k 的子數組個數

題意：
    給定一個長度 n 的整數陣列，做 m 次「區間加」操作；再給定一個 R×C 的矩陣，
    做 M 次「子矩陣加」操作；最後回答一次子矩陣求和查詢，並統計「和恰為 k 的連續子數組」個數。
    要求所有修改都做到 O(1) 攤銷、所有查詢都做到 O(1)。

思路：
    前綴和與差分是一對互逆的操作，是「區間批量修改 / 區間批量查詢」的基礎工具：

      - **差分**：`d[i] = a[i] − a[i−1]`（一維）。對原陣列做 `[l, r] += x`，
        等價於 `d[l] += x; d[r+1] −= x` 兩個**單點**修改 —— 區間加從 O(len) 降到 O(1)。
        所有操作做完後，對 d 求一次前綴和就還原出最終陣列（O(n)）。
        二維同理，只是式子從 2 項變成 4 項：
        `d[r1][c1] += v; d[r2+1][c1] −= v; d[r1][c2+1] −= v; d[r2+1][c2+1] += v`。
      - **前綴和**：`S[i] = a[0] + … + a[i−1]`。區間 `[l, r]` 的和 = `S[r+1] − S[l]`（O(1)）。
        二維：`S[i][j]` 是左上角 (0,0) 到 (i−1,j−1) 的矩形和，
        子矩陣 `[r1..r2][c1..c2]` 的和 = `S[r2+1][c2+1] − S[r1][c2+1] − S[r2+1][c1] + S[r1][c1]`
        （容斥：減掉兩塊多的，加回被減兩次的那塊）。
      - **和為 k 的子數組個數**（LeetCode 560）：區間和 `sum(i..j) = S[j+1] − S[i]`，
        問題變成「有多少對 (i, j) 滿足 `S[j+1] − S[i] = k`」。
        從左往右掃，用哈希表記錄**已經出現過的前綴和**的次數：
        當前前綴和為 s 時，需要補上 `s − k` 出現過幾次，就有幾個以當前位置結尾的合法子數組。
        初始要放 `{0: 1}`，代表「空前綴」，否則從下標 0 開始的子數組會漏掉。
        注意這個技巧**不能**用滑動窗口（陣列含負數時窗口不具單調性）。

    一句話記憶：**差分用來「改」，前綴和用來「查」**；
    兩者都把「區間」操作拆成端點上的單點操作，因此都只要 O(1) 每次。

輸入格式（stdin，數字以空白分隔即可）：
    n m
    a1 a2 … an
    l1 r1 x1        （m 行，1-indexed 閉區間）
    …
    k
    R C M
    b(1,1) … b(1,C) （R 行，每行 C 個整數）
    …
    r1 c1 r2 c2 v   （M 行，1-indexed 閉區間）
    …
    qr1 qc1 qr2 qc2 （最後一行：子矩陣求和查詢，1-indexed 閉區間）
輸出格式（stdout）：
    第 1 行：一維區間加後的陣列（n 個數，空格分隔；n = 0 時輸出空行）
    第 2 行：該陣列的前綴和（n 個數）
    第 3 行：和恰為 k 的連續子數組個數
    第 4 … R+3 行：二維子矩陣加後的矩陣（每行 C 個數，空格分隔；R = 0 時不輸出）
    最後一行：查詢子矩陣的元素和
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import List, Tuple


# ---------------------------------------------------------------- 輸入掃描
class Scanner:
    """把 stdin 切成整數 token 序列；輸入截斷時補 0，保證不會中途崩掉。"""

    def __init__(self, data: str) -> None:
        self.tok = data.split()
        self.pos = 0

    def nxt(self) -> int:
        if self.pos < len(self.tok):
            v = int(self.tok[self.pos])
            self.pos += 1
            return v
        return 0


# ---------------------------------------------------------------- 一維差分
def build_diff1(n: int, a: List[int]) -> List[int]:
    """由原陣列構造長度 n+1 的差分陣列（第 n 位是哨兵，吸收 r = n−1 的減項）。"""
    d = [0] * (n + 1)
    for i in range(n):
        d[i] += a[i]
        d[i + 1] -= a[i]
    return d


def range_add1(d: List[int], n: int, l: int, r: int, x: int) -> None:
    """一次區間加 [l, r] += x（0-indexed 閉區間）：只動兩個端點，O(1)。"""
    if n <= 0:
        return
    l = max(0, l)
    r = min(n - 1, r)
    if l > r:
        return
    d[l] += x
    d[r + 1] -= x                 # r+1 <= n，哨兵位一定存在


def restore1(d: List[int], n: int) -> List[int]:
    """對差分陣列求前綴和，還原出最終陣列。"""
    res = [0] * n
    cur = 0
    for i in range(n):
        cur += d[i]
        res[i] = cur
    return res


def range_add_naive(n: int, a: List[int], ops: List[Tuple[int, int, int]]) -> List[int]:
    """暴力版區間加，只在測試裡當基準。"""
    res = list(a)
    for l, r, x in ops:
        for i in range(max(0, l), min(n - 1, r) + 1):
            res[i] += x
    return res


def prefix_sum1(a: List[int]) -> List[int]:
    """前綴和：out[i] = a[0] + … + a[i]。"""
    out: List[int] = []
    s = 0
    for x in a:
        s += x
        out.append(s)
    return out


def count_subarray_sum_k(a: List[int], k: int) -> int:
    """和恰為 k 的連續子數組個數：前綴和 + 哈希表，O(n)。"""
    freq = {0: 1}                 # 空前綴
    s = 0
    ans = 0
    for x in a:
        s += x
        ans += freq.get(s - k, 0)
        freq[s] = freq.get(s, 0) + 1
    return ans


def count_subarray_sum_k_naive(a: List[int], k: int) -> int:
    """暴力枚舉所有子數組，O(n^2)，只在測試裡當基準。"""
    n = len(a)
    ans = 0
    for i in range(n):
        s = 0
        for j in range(i, n):
            s += a[j]
            if s == k:
                ans += 1
    return ans


# ---------------------------------------------------------------- 二維差分
def build_diff2(R: int, C: int, b: List[List[int]]) -> List[List[int]]:
    """由原矩陣構造 (R+1)×(C+1) 的二維差分陣列（多出來的一行一列是哨兵）。"""
    d = [[0] * (C + 1) for _ in range(R + 1)]
    for i in range(R):
        for j in range(C):
            v = b[i][j]
            d[i][j] += v
            d[i + 1][j] -= v
            d[i][j + 1] -= v
            d[i + 1][j + 1] += v
    return d


def submatrix_add2(d: List[List[int]], R: int, C: int,
                   r1: int, c1: int, r2: int, c2: int, v: int) -> None:
    """一次子矩陣加 [r1..r2]×[c1..c2] += v：動四個角，O(1)。"""
    if R <= 0 or C <= 0:
        return
    r1 = max(0, r1)
    c1 = max(0, c1)
    r2 = min(R - 1, r2)
    c2 = min(C - 1, c2)
    if r1 > r2 or c1 > c2:
        return
    d[r1][c1] += v
    d[r2 + 1][c1] -= v
    d[r1][c2 + 1] -= v
    d[r2 + 1][c2 + 1] += v       # 被減兩次，加回來


def restore2(d: List[List[int]], R: int, C: int) -> List[List[int]]:
    """二維前綴和還原：g[i][j] = d[i][j] + 上 + 左 − 左上。"""
    g = [[0] * C for _ in range(R)]
    for i in range(R):
        for j in range(C):
            up = g[i - 1][j] if i else 0
            left = g[i][j - 1] if j else 0
            diag = g[i - 1][j - 1] if (i and j) else 0
            g[i][j] = d[i][j] + up + left - diag
    return g


def prefix2(g: List[List[int]], R: int, C: int) -> List[List[int]]:
    """二維前綴和表 S，S[i][j] = 左上角 (0,0) 到 (i−1,j−1) 的矩形和。"""
    S = [[0] * (C + 1) for _ in range(R + 1)]
    for i in range(R):
        for j in range(C):
            S[i + 1][j + 1] = g[i][j] + S[i][j + 1] + S[i + 1][j] - S[i][j]
    return S


def submatrix_sum(S: List[List[int]], R: int, C: int,
                  r1: int, c1: int, r2: int, c2: int) -> int:
    """子矩陣求和（容斥），O(1)。"""
    if R <= 0 or C <= 0:
        return 0
    r1 = max(0, r1)
    c1 = max(0, c1)
    r2 = min(R - 1, r2)
    c2 = min(C - 1, c2)
    if r1 > r2 or c1 > c2:
        return 0
    return S[r2 + 1][c2 + 1] - S[r1][c2 + 1] - S[r2 + 1][c1] + S[r1][c1]


def submatrix_add_naive(R: int, C: int, b: List[List[int]],
                        ops: List[Tuple[int, int, int, int, int]]) -> List[List[int]]:
    """暴力版子矩陣加，只在測試裡當基準。"""
    g = [row[:] for row in b]
    for r1, c1, r2, c2, v in ops:
        for i in range(max(0, r1), min(R - 1, r2) + 1):
            for j in range(max(0, c1), min(C - 1, c2) + 1):
                g[i][j] += v
    return g


def submatrix_sum_naive(g: List[List[int]], R: int, C: int,
                        r1: int, c1: int, r2: int, c2: int) -> int:
    """暴力版子矩陣求和，只在測試裡當基準。"""
    if R <= 0 or C <= 0:
        return 0
    r1 = max(0, r1)
    c1 = max(0, c1)
    r2 = min(R - 1, r2)
    c2 = min(C - 1, c2)
    if r1 > r2 or c1 > c2:
        return 0
    return sum(g[i][j] for i in range(r1, r2 + 1) for j in range(c1, c2 + 1))


# ---------------------------------------------------------------- IO 模式
def run_io(data: str) -> None:
    sc = Scanner(data)
    n = sc.nxt()
    m = sc.nxt()
    a = [sc.nxt() for _ in range(max(0, n))]
    ops = []
    for _ in range(max(0, m)):
        l = sc.nxt() - 1          # 輸入是 1-indexed，轉 0-indexed
        r = sc.nxt() - 1
        x = sc.nxt()
        ops.append((l, r, x))

    d1 = build_diff1(max(0, n), a)
    for l, r, x in ops:
        range_add1(d1, max(0, n), l, r, x)
    arr = restore1(d1, max(0, n))
    print(" ".join(map(str, arr)))
    print(" ".join(map(str, prefix_sum1(arr))))

    k = sc.nxt()
    print(count_subarray_sum_k(arr, k))

    R = sc.nxt()
    C = sc.nxt()
    M = sc.nxt()
    b = [[sc.nxt() for _ in range(max(0, C))] for _ in range(max(0, R))]
    ops2 = []
    for _ in range(max(0, M)):
        r1 = sc.nxt() - 1
        c1 = sc.nxt() - 1
        r2 = sc.nxt() - 1
        c2 = sc.nxt() - 1
        v = sc.nxt()
        ops2.append((r1, c1, r2, c2, v))
    qr1 = sc.nxt() - 1
    qc1 = sc.nxt() - 1
    qr2 = sc.nxt() - 1
    qc2 = sc.nxt() - 1

    if R > 0 and C > 0:
        d2 = build_diff2(R, C, b)
        for op in ops2:
            submatrix_add2(d2, R, C, *op)
        grid = restore2(d2, R, C)
        for row in grid:
            print(" ".join(map(str, row)))
        print(submatrix_sum(prefix2(grid, R, C), R, C, qr1, qc1, qr2, qc2))
    else:
        print(0)


# ---------------------------------------------------------------- 測試
def run_tests() -> None:
    # 固定用例：README 的示例
    a = [1, 2, 3, 4, 5]
    ops = [(1, 3, 10), (3, 4, 1)]         # 0-indexed：[1,3] += 10，[3,4] += 1
    d = build_diff1(5, a)
    for l, r, x in ops:
        range_add1(d, 5, l, r, x)
    got = restore1(d, 5)
    assert got == range_add_naive(5, a, ops) == [1, 12, 13, 15, 6]
    assert prefix_sum1(got) == [1, 13, 26, 41, 47]
    assert count_subarray_sum_k(got, 13) == 2              # [1, 12] 與 [13]
    assert count_subarray_sum_k(got, 13) == count_subarray_sum_k_naive(got, 13)

    # 空陣列
    assert restore1(build_diff1(0, []), 0) == []
    assert prefix_sum1([]) == []
    assert count_subarray_sum_k([], 0) == 0
    assert count_subarray_sum_k([], 7) == 0

    # 前綴和為 0 的經典陷阱：全 0 陣列 + k = 0 → 每個子數組都算
    assert count_subarray_sum_k([0, 0, 0], 0) == 6          # 3*4/2
    assert count_subarray_sum_k([3, 4, -7, 1], 0) == 1      # 只有 [3, 4, -7]
    assert count_subarray_sum_k([3, 4, -7, 1], 0) == count_subarray_sum_k_naive([3, 4, -7, 1], 0)

    # 二維固定用例：3×3 全 0，兩次子矩陣加
    b = [[0] * 3 for _ in range(3)]
    ops2 = [(0, 0, 1, 1, 5), (1, 1, 2, 2, 3)]               # 0-indexed
    d2 = build_diff2(3, 3, b)
    for op in ops2:
        submatrix_add2(d2, 3, 3, *op)
    g = restore2(d2, 3, 3)
    assert g == submatrix_add_naive(3, 3, b, ops2)
    assert g == [[5, 5, 0], [5, 8, 3], [0, 3, 3]]
    S = prefix2(g, 3, 3)
    assert submatrix_sum(S, 3, 3, 0, 0, 2, 2) == submatrix_sum_naive(g, 3, 3, 0, 0, 2, 2) == 32
    assert submatrix_sum(S, 3, 3, 1, 1, 1, 1) == 8
    assert submatrix_sum(S, 3, 3, 0, 0, 0, 2) == 10

    # 隨機對拍
    random.seed(20261002)
    for _ in range(500):
        n = random.randint(0, 8)
        m = random.randint(0, 5)
        a = [random.randint(-9, 9) for _ in range(n)]
        ops = []
        for _ in range(m):
            if n == 0:
                ops.append((0, -1, random.randint(-5, 5)))   # 空區間，應被忽略
            else:
                l = random.randint(0, n - 1)
                r = random.randint(l, n - 1)
                ops.append((l, r, random.randint(-9, 9)))
        # 偶爾混入越界區間，驗證裁剪邏輯
        if n > 0 and random.random() < 0.3:
            ops.append((-3, n + 3, 7))

        d = build_diff1(n, a)
        for l, r, x in ops:
            range_add1(d, n, l, r, x)
        arr = restore1(d, n)
        assert arr == range_add_naive(n, a, ops)

        k = random.randint(-12, 12)
        assert count_subarray_sum_k(arr, k) == count_subarray_sum_k_naive(arr, k)

        R = random.randint(0, 4)
        C = random.randint(0, 4)
        b = [[random.randint(-9, 9) for _ in range(C)] for _ in range(R)]
        ops2 = []
        for _ in range(random.randint(0, 4)):
            if R == 0 or C == 0:
                ops2.append((0, 0, -1, -1, random.randint(-5, 5)))
            else:
                r1 = random.randint(0, R - 1)
                r2 = random.randint(r1, R - 1)
                c1 = random.randint(0, C - 1)
                c2 = random.randint(c1, C - 1)
                ops2.append((r1, c1, r2, c2, random.randint(-9, 9)))
        d2 = build_diff2(R, C, b)
        for op in ops2:
            submatrix_add2(d2, R, C, *op)
        g = restore2(d2, R, C)
        if R > 0 and C > 0:
            assert g == submatrix_add_naive(R, C, b, ops2)
            S = prefix2(g, R, C)
            for _ in range(3):
                r1 = random.randint(0, R - 1)
                r2 = random.randint(r1, R - 1)
                c1 = random.randint(0, C - 1)
                c2 = random.randint(c1, C - 1)
                assert submatrix_sum(S, R, C, r1, c1, r2, c2) == \
                    submatrix_sum_naive(g, R, C, r1, c1, r2, c2)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
