"""Z 函數（Z Algorithm）與擴展 KMP

題意：
    給定字串 s（文本）與 p（模式串），求：
      1. s 的 Z 陣列：z[i] = s 與 s[i:] 的最長公共前綴長度（定義 z[0] = |s|）
      2. p 的 Z 陣列
      3. p 在 s 中所有出現的起始位置（0-indexed）
      4. 擴展 KMP 的 ext 陣列：ext[i] = p 與 s[i:] 的最長公共前綴長度
      5. s 的最短週期（最小的 k 使 s 由某個長度 k 的塊重複整數次構成）
    以上全部同時用「O(n) 演算法」與「O(n²) 暴力」計算並互相驗證。

思路：
    ### Z 函數：用已經算過的區間當跳板
    維護一個「Z-box」[l, r]，它是當前已知的、最靠右的匹配區間，
    也就是 s[l..r] == s[0..r-l]。
    算 z[i] 時：
      - 若 i > r：Z-box 蓋不到 i，只能從 0 開始老實比（while 迴圈）。
      - 若 i ≤ r：因為 s[i..r] 已經等於 s[0..r-i]，
        我們可以**直接複製** z[i-l]（同一個相對位置先前算過的值），
        但最多只能複製到 r 為止 → 取 `min(r - i + 1, z[i - l])`，
        再從這個起點繼續往後比。
    關鍵是**每個字元最多被 while 迴圈比對一次**：一旦匹配超過 r，
    就把 Z-box 右端推進到新的位置，而 r 只增不減。故總複雜度 O(n)。
    這跟 KMP 的 prefix function 是對偶的兩套方法：KMP 預處理模式串、
    掃描時不回頭；Z 函數則把「自己和自己的每個後綴」一次算完。

    ### 擴展 KMP：一次 Z 搞定所有後綴
    要求 p 與 s 的每個後綴 s[i:] 的 LCP，只要在
    `p + '#' + s`（'#' 是不出現在字串裡的哨兵）上跑一次 Z 函數：
    位置 |p|+1+i 對應 s[i:]，取其 z 值再跟 |p| 取 min 即可（哨兵會擋住越界）。
    這比「對每個 i 都跑一次 KMP」快得多，是 O(|p| + |s|)。

    ### 最短週期
    若 s 有週期 k（且 k | n），則 s[k:] == s[:n-k]，也就是 z[k] == n - k。
    反過來也成立。所以枚舉 n 的所有因數 k，取第一個滿足 z[k] == n-k 的就是答案。
    這是 Z 函數在「字串週期性」問題上的標準應用（另一個常見應用是
    判斷字串是否為某個字串的冪、以及字串壓縮）。

應用場景：
    - 單模式串匹配：O(n+m)，比 KMP 更容易寫對，且一次拿到所有後綴的匹配長度
    - 字串週期 / 冪次判定、字串壓縮（求最短週期）、重複子串偵測
    - 擴展 KMP 用於：最長公共前綴的大量查詢、DNA 序列比對、
      抄襲檢測中的「從每個位置開始能匹配多長」
    - 搭配後綴結構可做「不同子串計數」「最長重複子串」等進階題

複雜度：
    Z 陣列                    O(n) 時間 / O(n) 空間
    擴展 KMP（所有後綴 LCP）  O(|p| + |s|) 時間 / O(|p| + |s|) 空間
    模式匹配（所有出現位置）  O(|p| + |s| + 出現次數)
    最短週期                  O(n)（枚舉因數 + Z 陣列）
"""

import random
import re
import sys
from typing import List

WORD_RE = re.compile(r"^[a-z]+$")

MAXL = 2000          # 單串最大長度（IO 模式下的截斷上限）
SENTINEL = "#"       # 哨兵字元：只可能出現在拼接串裡，不會在輸入中出現


# ---------------------------------------------------------------- 核心演算法


def z_algorithm(s: str) -> List[int]:
    """Z 函數：z[i] = LCP(s, s[i:])，約定 z[0] = len(s)。O(n)。"""
    n = len(s)
    if n == 0:
        return []
    z = [0] * n
    z[0] = n
    l, r = 0, 0                       # 當前最右的 Z-box [l, r]（閉區間）
    for i in range(1, n):
        if i <= r:
            # s[i..r] 已知等於 s[0..r-i]，直接複製舊答案，但不要超出 r
            z[i] = min(r - i + 1, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] - 1 > r:          # 匹配突破了舊的右端，推進 Z-box
            l, r = i, i + z[i] - 1
    return z


def extended_kmp(text: str, pat: str) -> List[int]:
    """擴展 KMP：回傳 ext[i] = LCP(pat, text[i:])，長度為 len(text)。"""
    n = len(text)
    if n == 0:
        return []
    m = len(pat)
    if m == 0:
        return [0] * n
    z = z_algorithm(pat + SENTINEL + text)
    return [min(z[m + 1 + i], m) for i in range(n)]


def z_match(text: str, pat: str) -> List[int]:
    """pat 在 text 中所有出現的起始位置（0-indexed，遞增）。"""
    n, m = len(text), len(pat)
    if m == 0:
        return list(range(n + 1))
    ext = extended_kmp(text, pat)
    return [i for i in range(n - m + 1) if ext[i] == m]


def smallest_period(s: str) -> int:
    """s 的最短週期：最小的 k ≥ 1 使 s 恰好由長度 k 的塊重複整數次構成。"""
    n = len(s)
    if n == 0:
        return 0
    z = z_algorithm(s)
    for k in range(1, n):             # k = n 必然是週期，留作兜底（z[n] 不存在）
        if n % k == 0 and z[k] == n - k:
            return k
    return n                          # 恆成立：k = n 一定是一個週期


# ---------------------------------------------------------------- 暴力對拍


def z_brute(s: str) -> List[int]:
    n = len(s)
    res = [0] * n
    for i in range(n):
        k = 0
        while i + k < n and s[k] == s[i + k]:
            k += 1
        res[i] = k
    if n:
        res[0] = n
    return res


def ext_brute(text: str, pat: str) -> List[int]:
    n, m = len(text), len(pat)
    res = []
    for i in range(n):
        k = 0
        while k < m and i + k < n and pat[k] == text[i + k]:
            k += 1
        res.append(k)
    return res


def match_brute(text: str, pat: str) -> List[int]:
    n, m = len(text), len(pat)
    if m == 0:
        return list(range(n + 1))
    return [i for i in range(n - m + 1) if text[i:i + m] == pat]


def period_brute(s: str) -> int:
    n = len(s)
    for k in range(1, n + 1):
        if n % k:
            continue
        if all(s[i] == s[i % k] for i in range(n)):
            return k
    return n


# ---------------------------------------------------------------- 輸出


def join(nums: List[int]) -> str:
    return " ".join(str(v) for v in nums)


def run(text: str, pat: str) -> str:
    """回傳 6 行輸出；第 6 行是「演算法 vs 暴力」的一致性旗標。"""
    zs = z_algorithm(text)
    zp = z_algorithm(pat)
    occ = z_match(text, pat)
    ext = extended_kmp(text, pat)
    per = smallest_period(text)

    ok = 1
    if zs != z_brute(text):
        ok = 0
    if zp != z_brute(pat):
        ok = 0
    if occ != match_brute(text, pat):
        ok = 0
    if ext != ext_brute(text, pat):
        ok = 0
    if per != period_brute(text):
        ok = 0
    # 出現位置必須與 ext 陣列自洽：ext[i] == |p| ⟺ i 是合法起點且完全匹配
    m = len(pat)
    if m > 0 and [i for i in range(len(text) - m + 1) if ext[i] == m] != occ:
        ok = 0

    return "\n".join([join(zs), join(zp), join(occ), join(ext), str(per), str(ok)])


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例 ----
    assert z_algorithm("aaaaa") == [5, 4, 3, 2, 1]
    assert z_algorithm("ababab") == [6, 0, 4, 0, 2, 0]
    assert z_algorithm("") == []
    assert z_algorithm("a") == [1]
    assert z_algorithm("abcababc") == [8, 0, 0, 2, 0, 3, 0, 0]

    assert z_match("abababab", "abab") == [0, 2, 4]
    assert z_match("aaaaa", "aa") == [0, 1, 2, 3]
    assert z_match("abc", "xyz") == []
    assert z_match("abc", "") == [0, 1, 2, 3]
    assert z_match("", "a") == []
    assert z_match("", "") == [0]

    assert extended_kmp("aaaaa", "aa") == [2, 2, 2, 2, 1]
    assert extended_kmp("abcabc", "abc") == [3, 0, 0, 3, 0, 0]
    assert extended_kmp("abc", "") == [0, 0, 0]

    assert smallest_period("ababab") == 2
    assert smallest_period("aaaa") == 1
    assert smallest_period("abcab") == 5
    assert smallest_period("abcabcabc") == 3
    assert smallest_period("") == 0

    # ---- 隨機對拍 ----
    rng = random.Random(20261015)
    for _ in range(600):
        alpha = rng.choice(["ab", "abc", "a", "abcdefgh"])
        text = "".join(rng.choice(alpha) for _ in range(rng.randint(0, 60)))
        pat = "".join(rng.choice(alpha) for _ in range(rng.randint(0, 12)))
        assert z_algorithm(text) == z_brute(text)
        assert z_algorithm(pat) == z_brute(pat)
        assert extended_kmp(text, pat) == ext_brute(text, pat)
        assert z_match(text, pat) == match_brute(text, pat)
        assert smallest_period(text) == period_brute(text)
        assert run(text, pat).endswith("\n1")

    # ---- 長串只跑線性演算法（不跑 O(n²) 暴力）----
    long_text = ("ab" * 800)
    assert smallest_period(long_text) == 2
    assert z_match(long_text, "abab") == list(range(0, len(long_text) - 4 + 1, 2))
    assert z_algorithm("a" * 1500)[1] == 1499


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出（latin-1），避免 Windows 把 \\n 換成 \\r\\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def clean(tok: str) -> str:
    """只接受純小寫字母的 token，否則視為空字串；長度截斷到 MAXL。"""
    if not WORD_RE.match(tok):
        return ""
    return tok[:MAXL]


def run_io(raw: str) -> None:
    toks = raw.split()
    text = clean(toks[0]) if len(toks) >= 1 else ""
    pat = clean(toks[1]) if len(toks) >= 2 else ""
    emit(run(text, pat) + "\n")


def main() -> None:
    sys.setrecursionlimit(1000000)
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
