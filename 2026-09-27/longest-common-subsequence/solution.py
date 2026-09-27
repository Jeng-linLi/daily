"""最长公共子序列（LCS，动态规划 + 回溯还原）

题意：
    给定两个字符串 a、b，求它们的最长公共子序列的长度，并输出一条达到该长度的
    子序列（子序列不要求连续，只要求保持相对顺序）。

思路：
    定义 dp[i][j] = a 的前 i 个字符与 b 的前 j 个字符的 LCS 长度。
    看最后一对字符 a[i-1] 与 b[j-1]，只有两种情况：
      - 相等：这个字符一定可以接在 a[:i-1] 与 b[:j-1] 的 LCS 后面，
              所以 dp[i][j] = dp[i-1][j-1] + 1。
              （可以证明最优解总能取这个字符：若某个最优解没用它，把该解的
                最后一个字符换成它依然合法且长度不变。）
      - 不等：它俩不可能同时出现在同一个匹配里，
              dp[i][j] = max(dp[i-1][j], dp[i][j-1])。
    边界 dp[0][j] = dp[i][0] = 0（一边为空，LCS 长度为 0）。
    填表顺序按 i、j 从小到大，保证用到的状态都已经算好。

    只求长度时空间可以压到两行（甚至一行 + 对角线变量），因为 dp[i][*] 只依赖
    dp[i-1][*]；但要**还原具体方案**必须保留整张表，再从 dp[m][n] 往回走：
      - a[i-1] == b[j-1]：这个字符属于 LCS，记下来，i、j 都减一；
      - 否则往 dp 值大的方向走（相等时优先走 i，即丢弃 a[i-1]）。
    回溯是倒着走的，所以收集到的字符要反转一次。

    注意 LCS 通常**不唯一**（"abcbdab" 与 "bdcaba" 有多条长度为 4 的解），
    回溯只保证给出其中一条；测试因此断言「长度等于最优值」且
    「该串确实是两边的公共子序列」，而不锁死具体是哪一条。

输入格式（stdin）：
    第一行：字符串 a
    第二行：字符串 b（可以是空行）
输出格式（stdout）：
    第一行：LCS 长度
    第二行：一条达到该长度的公共子序列（长度为 0 时输出空行）
无 stdin 输入时运行内置断言测试并输出 `all tests passed`。
"""

import sys
from functools import lru_cache
from typing import List, Tuple


def lcs_length(a: str, b: str) -> int:
    """两行滚动数组，只求长度。时间 O(m*n)，空间 O(min(m, n))。"""
    if len(a) < len(b):
        a, b = b, a          # 让 b 成为较短的那个，滚动数组更省空间
    m, n = len(a), len(b)
    prev = [0] * (n + 1)
    cur = [0] * (n + 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
            else:
                cur[j] = prev[j] if prev[j] >= cur[j - 1] else cur[j - 1]
        prev, cur = cur, prev   # 交换，下一行复用上一行的空间
        cur[0] = 0
    return prev[n]


def lcs_with_string(a: str, b: str) -> Tuple[int, str]:
    """保留完整 dp 表并回溯出一条 LCS。时间 O(m*n)，空间 O(m*n)。"""
    m, n = len(a), len(b)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    chars: List[str] = []
    i, j = m, n
    while i > 0 and j > 0:
        if a[i - 1] == b[j - 1]:
            chars.append(a[i - 1])          # 这个字符属于 LCS
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1                          # 丢弃 a[i-1]（相等时优先走 i）
        else:
            j -= 1                          # 丢弃 b[j-1]
    chars.reverse()
    return dp[m][n], "".join(chars)


def is_subsequence(sub: str, s: str) -> bool:
    """判断 sub 是否为 s 的子序列。"""
    it = iter(s)
    return all(ch in it for ch in sub)


def lcs_brute(a: str, b: str) -> int:
    """对照用的带记忆化递归（指数级搜索 + 剪枝），仅用于小规模测试验证。"""

    @lru_cache(maxsize=None)
    def go(i: int, j: int) -> int:
        if i == 0 or j == 0:
            return 0
        if a[i - 1] == b[j - 1]:
            return go(i - 1, j - 1) + 1
        return max(go(i - 1, j), go(i, j - 1))

    res = go(len(a), len(b))
    go.cache_clear()
    return res


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    lines = data.splitlines()
    a = lines[0] if len(lines) > 0 else ""
    b = lines[1] if len(lines) > 1 else ""
    length, sub = lcs_with_string(a, b)
    print(length)
    print(sub)


def run_tests() -> None:
    # README 示例：abcde 与 ace 的 LCS 是 ace，长度 3
    length, sub = lcs_with_string("abcde", "ace")
    assert length == 3
    assert sub == "ace"
    assert lcs_length("abcde", "ace") == 3
    assert lcs_brute("abcde", "ace") == 3

    # 完全相同：LCS 就是自身
    assert lcs_length("abc", "abc") == 3
    assert lcs_with_string("abc", "abc") == (3, "abc")

    # 没有公共字符：空串
    assert lcs_length("abc", "def") == 0
    assert lcs_with_string("abc", "def") == (0, "")

    # 一边为空
    assert lcs_length("", "abc") == 0
    assert lcs_with_string("", "abc") == (0, "")
    assert lcs_with_string("abc", "") == (0, "")
    assert lcs_with_string("", "") == (0, "")

    # 经典用例：长度为 4（"bdab" / "bcba" 等都算对，只断言长度与合法性）
    length, sub = lcs_with_string("abcbdab", "bdcaba")
    assert length == 4
    assert len(sub) == 4
    assert is_subsequence(sub, "abcbdab")
    assert is_subsequence(sub, "bdcaba")
    assert lcs_brute("abcbdab", "bdcaba") == 4

    # a 是 b 的子序列：LCS 就是 a
    assert lcs_length("ace", "abcde") == 3
    assert lcs_with_string("ace", "abcde")[1] == "ace"

    # 大小写敏感
    assert lcs_length("Abc", "abc") == 2

    # 重复字符
    assert lcs_length("aaaa", "aa") == 2
    assert lcs_length("aab", "aba") == 2

    import random

    random.seed(20260927)
    alphabet = "abc"

    # 随机对拍：两行滚动版 = 完整表版 = 记忆化暴力版，且还原出的串确实合法
    for _ in range(300):
        a = "".join(random.choice(alphabet) for _ in range(random.randint(0, 8)))
        b = "".join(random.choice(alphabet) for _ in range(random.randint(0, 8)))
        d1 = lcs_length(a, b)
        d2, sub = lcs_with_string(a, b)
        d3 = lcs_brute(a, b)
        assert d1 == d2 == d3                       # 三个版本长度一致
        assert len(sub) == d2                       # 还原出的串长度就是最优值
        assert is_subsequence(sub, a)               # 确实是 a 的子序列
        assert is_subsequence(sub, b)               # 确实是 b 的子序列
        assert 0 <= d2 <= min(len(a), len(b))       # 长度落在合理区间内

    # 对称性与上界性质
    for _ in range(100):
        a = "".join(random.choice(alphabet) for _ in range(random.randint(0, 6)))
        b = "".join(random.choice(alphabet) for _ in range(random.randint(0, 6)))
        assert lcs_length(a, b) == lcs_length(b, a)
        assert lcs_with_string(a, b)[0] == lcs_with_string(b, a)[0]
        assert 0 <= lcs_length(a, b) <= min(len(a), len(b))

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
