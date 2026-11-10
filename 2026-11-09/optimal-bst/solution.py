# -*- coding: utf-8 -*-
"""最優二叉搜尋樹 Optimal BST 與 Knuth 四邊形不等式優化

題意：n 個已排序的鍵 1..n，第 i 個鍵被搜尋的權重為 p[i]。
     建一棵二叉搜尋樹（中序即 1..n），第 i 個鍵的深度為其節點深度（根為 1），
     總成本 = Σ p[i] × depth(key i)。求最小總成本與整棵樹的根。

樸素區間 DP（O(n³)）：
  dp[i][j] = 鍵 i..j 建成子樹的最小成本，w[i][j] = Σ p[i..j]
  dp[i][j] = min_{r=i..j} { dp[i][r-1] + dp[r+1][j] } + w[i][j]

Knuth 優化（O(n²)）：
  權重前綴和滿足四邊形不等式（QI）與區間單調性 ⟹ 最優根 s[i][j] 滿足
  s[i][j-1] ≤ s[i][j] ≤ s[i+1][j]，列舉根時只需在該窗口內找，
  每個狀態平均攤還 O(1) 次比較，總複雜度從 O(n³) 降到 O(n²)。

IO 模式（stdin 有內容時）：
  n
  n 個整數 p[1..n]（1 ≤ n ≤ 3000，0 ≤ p[i] ≤ 10^9）
  輸出兩行：第一行最小成本，第二行整樹根 s[1][n]（1-based）。
  n = 0 時輸出 "0 0"。
無 stdin 時跑內建斷言：O(n³) 樸素 DP 與 Knuth O(n²) 版對拍
（成本必須相等），並檢查 s 陣列的單調性邊界，最後輸出 all tests passed。
"""
import sys
from random import Random

INF = float("inf")


def optbst_brute(p):
    """樸素 O(n³) 區間 DP，平手取最小根。回傳 (最小成本, 整樹根)。"""
    n = len(p) - 1  # p 為 1-indexed
    if n == 0:
        return 0, 0
    w = [[0] * (n + 1) for _ in range(n + 1)]
    dp = [[0] * (n + 2) for _ in range(n + 2)]
    for i in range(1, n + 1):
        w[i][i] = p[i]
        dp[i][i] = p[i]
    for length in range(2, n + 1):
        for i in range(1, n - length + 2):
            j = i + length - 1
            w[i][j] = w[i][j - 1] + p[j]
            best = INF
            bestr = i
            for r in range(i, j + 1):
                cost = dp[i][r - 1] + dp[r + 1][j]
                if cost < best:
                    best = cost
                    bestr = r
            dp[i][j] = best + w[i][j]
    return dp[1][n], 1 if n == 1 else _brute_root(p, 1, n, w, dp)


def _brute_root(p, i, j, w, dp):
    """回溯：取 dp 區間內平手時編號最小的最優根（用 dp 表重算）。"""
    best = INF
    bestr = i
    for r in range(i, j + 1):
        cost = dp[i][r - 1] + dp[r + 1][j]
        if cost < best:
            best = cost
            bestr = r
    return bestr


def optbst_knuth(p):
    """Knuth O(n²) 版。回傳 (最小成本, s[1][n])。"""
    n = len(p) - 1
    if n == 0:
        return 0, 0
    w = [[0] * (n + 1) for _ in range(n + 1)]
    dp = [[0] * (n + 2) for _ in range(n + 2)]
    s = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        w[i][i] = p[i]
        dp[i][i] = p[i]
        s[i][i] = i
    for length in range(2, n + 1):
        for i in range(1, n - length + 2):
            j = i + length - 1
            w[i][j] = w[i][j - 1] + p[j]
            lo = s[i][j - 1]
            hi = s[i + 1][j]
            best = INF
            bestr = lo
            for r in range(lo, hi + 1):
                cost = dp[i][r - 1] + dp[r + 1][j]
                if cost < best:
                    best = cost
                    bestr = r
            dp[i][j] = best + w[i][j]
            s[i][j] = bestr
    return dp[1][n], s[1][n]


def run_tests():
    rng = Random(20261109)

    # 1) 固定用例（期望值由實跑取得）
    p = [0, 2, 3, 1, 4]  # 1-indexed
    cost_b, root_b = optbst_brute(p)
    cost_k, root_k = optbst_knuth(p)
    assert cost_b == cost_k == 18, (cost_b, cost_k)
    assert root_b == root_k == 2, (root_b, root_k)

    p = [0, 5]
    assert optbst_brute(p) == optbst_knuth(p) == (5, 1)

    p = [0, 1, 1, 1]
    cost_b, root_b = optbst_brute(p)
    cost_k, root_k = optbst_knuth(p)
    assert cost_b == cost_k == 5, (cost_b, cost_k)
    assert root_b == root_k == 2, (root_b, root_k)

    # 2) 隨機對拍：Knuth 成本必須與樸素 O(n³) 相等；Knuth 的根必為最優根
    for tc in range(300):
        n = rng.randint(1, 25)
        p = [0] + [rng.randint(0, 20) for _ in range(n)]
        cost_b, _ = optbst_brute(p)
        cost_k, root_k = optbst_knuth(p)
        assert cost_b == cost_k, "case %d: %d vs %d" % (tc, cost_b, cost_k)
        # 以 Knuth 根 r 組樹：成本 = w[1..n] + 最優 dp[i][r-1] + 最優 dp[r+1][j]
        assert _root_cost(p, root_k) == cost_b, "case %d root not optimal" % tc
    assert _root_cost([0, 2, 3, 1, 4], 2) == 18
    assert _root_cost([0, 2, 3, 1, 4], 3) == 21

    # 3) s 陣列單調性：s[i][j-1] ≤ s[i][j] ≤ s[i+1][j]
    for tc in range(50):
        n = rng.randint(2, 30)
        p = [0] + [rng.randint(0, 20) for _ in range(n)]
        _check_s_monotone(p)


def _root_cost(p, root):
    """以指定根 r 組整棵樹（左右子樹各自取最優）的總成本，用於驗證根最優。"""
    n = len(p) - 1
    w = [[0] * (n + 1) for _ in range(n + 1)]
    dp = [[0] * (n + 2) for _ in range(n + 2)]
    for i in range(1, n + 1):
        w[i][i] = p[i]
        dp[i][i] = p[i]
    for length in range(2, n + 1):
        for i in range(1, n - length + 2):
            j = i + length - 1
            w[i][j] = w[i][j - 1] + p[j]
            dp[i][j] = min(dp[i][r - 1] + dp[r + 1][j]
                           for r in range(i, j + 1)) + w[i][j]
    return w[1][n] + dp[1][root - 1] + dp[root + 1][n]


def _check_s_monotone(p):
    """重建 Knuth 的 s 表並檢查單調性（僅檢查，不回傳）。"""
    n = len(p) - 1
    w = [[0] * (n + 1) for _ in range(n + 1)]
    dp = [[0] * (n + 2) for _ in range(n + 2)]
    s = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        w[i][i] = p[i]
        dp[i][i] = p[i]
        s[i][i] = i
    for length in range(2, n + 1):
        for i in range(1, n - length + 2):
            j = i + length - 1
            w[i][j] = w[i][j - 1] + p[j]
            lo, hi = s[i][j - 1], s[i + 1][j]
            best, bestr = INF, lo
            for r in range(lo, hi + 1):
                cost = dp[i][r - 1] + dp[r + 1][j]
                if cost < best:
                    best, bestr = cost, r
            dp[i][j] = best + w[i][j]
            s[i][j] = bestr
            assert s[i][j - 1] <= s[i][j] <= s[i + 1][j], "s monotonicity broken"


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        print("all tests passed")
        return
    n = int(data[0])
    if n == 0:
        print("0 0")
        return
    p = [0] + [int(x) for x in data[1:1 + n]]
    cost, root = optbst_knuth(p)
    print(cost)
    print(root)


if __name__ == "__main__":
    main()
