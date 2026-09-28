"""最長公共子序列（LCS，動態規劃 + 回溯還原）

題意：
    給定兩個字符串 a、b，求它們的最長公共子序列的長度，並輸出一條達到該長度的
    子序列（子序列不要求連續，只要求保持相對順序）。

思路：
    定義 dp[i][j] = a 的前 i 個字符與 b 的前 j 個字符的 LCS 長度。
    看最後一對字符 a[i-1] 與 b[j-1]，只有兩種情況：
      - 相等：這個字符一定可以接在 a[:i-1] 與 b[:j-1] 的 LCS 後面，
              所以 dp[i][j] = dp[i-1][j-1] + 1。
              （可以證明最優解總能取這個字符：若某個最優解沒用它，把該解的
                最後一個字符換成它依然合法且長度不變。）
      - 不等：它倆不可能同時出現在同一個匹配裏，
              dp[i][j] = max(dp[i-1][j], dp[i][j-1])。
    邊界 dp[0][j] = dp[i][0] = 0（一邊爲空，LCS 長度爲 0）。
    填表順序按 i、j 從小到大，保證用到的狀態都已經算好。

    只求長度時空間可以壓到兩行（甚至一行 + 對角線變量），因爲 dp[i][*] 只依賴
    dp[i-1][*]；但要**還原具體方案**必須保留整張表，再從 dp[m][n] 往回走：
      - a[i-1] == b[j-1]：這個字符屬於 LCS，記下來，i、j 都減一；
      - 否則往 dp 值大的方向走（相等時優先走 i，即丟棄 a[i-1]）。
    回溯是倒着走的，所以收集到的字符要反轉一次。

    注意 LCS 通常**不唯一**（"abcbdab" 與 "bdcaba" 有多條長度爲 4 的解），
    回溯只保證給出其中一條；測試因此斷言「長度等於最優值」且
    「該串確實是兩邊的公共子序列」，而不鎖死具體是哪一條。

輸入格式（stdin）：
    第一行：字符串 a
    第二行：字符串 b（可以是空行）
輸出格式（stdout）：
    第一行：LCS 長度
    第二行：一條達到該長度的公共子序列（長度爲 0 時輸出空行）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import sys
from functools import lru_cache
from typing import List, Tuple


def lcs_length(a: str, b: str) -> int:
    """兩行滾動數組，只求長度。時間 O(m*n)，空間 O(min(m, n))。"""
    if len(a) < len(b):
        a, b = b, a          # 讓 b 成爲較短的那個，滾動數組更省空間
    m, n = len(a), len(b)
    prev = [0] * (n + 1)
    cur = [0] * (n + 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
            else:
                cur[j] = prev[j] if prev[j] >= cur[j - 1] else cur[j - 1]
        prev, cur = cur, prev   # 交換，下一行復用上一行的空間
        cur[0] = 0
    return prev[n]


def lcs_with_string(a: str, b: str) -> Tuple[int, str]:
    """保留完整 dp 表並回溯出一條 LCS。時間 O(m*n)，空間 O(m*n)。"""
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
            chars.append(a[i - 1])          # 這個字符屬於 LCS
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1                          # 丟棄 a[i-1]（相等時優先走 i）
        else:
            j -= 1                          # 丟棄 b[j-1]
    chars.reverse()
    return dp[m][n], "".join(chars)


def is_subsequence(sub: str, s: str) -> bool:
    """判斷 sub 是否爲 s 的子序列。"""
    it = iter(s)
    return all(ch in it for ch in sub)


def lcs_brute(a: str, b: str) -> int:
    """對照用的帶記憶化遞歸（指數級搜索 + 剪枝），僅用於小規模測試驗證。"""

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
    """按統一輸入輸出格式處理 stdin 數據。"""
    lines = data.splitlines()
    a = lines[0] if len(lines) > 0 else ""
    b = lines[1] if len(lines) > 1 else ""
    length, sub = lcs_with_string(a, b)
    print(length)
    print(sub)


def run_tests() -> None:
    # README 示例：abcde 與 ace 的 LCS 是 ace，長度 3
    length, sub = lcs_with_string("abcde", "ace")
    assert length == 3
    assert sub == "ace"
    assert lcs_length("abcde", "ace") == 3
    assert lcs_brute("abcde", "ace") == 3

    # 完全相同：LCS 就是自身
    assert lcs_length("abc", "abc") == 3
    assert lcs_with_string("abc", "abc") == (3, "abc")

    # 沒有公共字符：空串
    assert lcs_length("abc", "def") == 0
    assert lcs_with_string("abc", "def") == (0, "")

    # 一邊爲空
    assert lcs_length("", "abc") == 0
    assert lcs_with_string("", "abc") == (0, "")
    assert lcs_with_string("abc", "") == (0, "")
    assert lcs_with_string("", "") == (0, "")

    # 經典用例：長度爲 4（"bdab" / "bcba" 等都算對，只斷言長度與合法性）
    length, sub = lcs_with_string("abcbdab", "bdcaba")
    assert length == 4
    assert len(sub) == 4
    assert is_subsequence(sub, "abcbdab")
    assert is_subsequence(sub, "bdcaba")
    assert lcs_brute("abcbdab", "bdcaba") == 4

    # a 是 b 的子序列：LCS 就是 a
    assert lcs_length("ace", "abcde") == 3
    assert lcs_with_string("ace", "abcde")[1] == "ace"

    # 大小寫敏感
    assert lcs_length("Abc", "abc") == 2

    # 重複字符
    assert lcs_length("aaaa", "aa") == 2
    assert lcs_length("aab", "aba") == 2

    import random

    random.seed(20260927)
    alphabet = "abc"

    # 隨機對拍：兩行滾動版 = 完整表版 = 記憶化暴力版，且還原出的串確實合法
    for _ in range(300):
        a = "".join(random.choice(alphabet) for _ in range(random.randint(0, 8)))
        b = "".join(random.choice(alphabet) for _ in range(random.randint(0, 8)))
        d1 = lcs_length(a, b)
        d2, sub = lcs_with_string(a, b)
        d3 = lcs_brute(a, b)
        assert d1 == d2 == d3                       # 三個版本長度一致
        assert len(sub) == d2                       # 還原出的串長度就是最優值
        assert is_subsequence(sub, a)               # 確實是 a 的子序列
        assert is_subsequence(sub, b)               # 確實是 b 的子序列
        assert 0 <= d2 <= min(len(a), len(b))       # 長度落在合理區間內

    # 對稱性與上界性質
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
