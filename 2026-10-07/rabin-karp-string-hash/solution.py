"""Rabin-Karp 與滾動哈希（子串匹配 / 最長公共子串 / 最長回文子串 / 不同子串計數 / k-gram 相似度）

題意：
    給定兩個字符串 text 與 pattern（ASCII），用「多項式滾動哈希（Rabin-Karp fingerprint）」解決五個問題：
      1. pattern 在 text 中的所有出現位置（與 KMP、暴力法三方對拍）；
      2. 兩串的最長公共子串（長度 + 最靠左的那一個）；
      3. text 的最長回文子串（長度 + 最靠左的那一個）；
      4. 兩串各自的不同子串個數；
      5. 兩串的 3-gram Jaccard 相似度（**應用**：文件指紋 / 抄襲檢測），以最簡分數 `p/q` 輸出。

思路：
    ### 多項式滾動哈希
    把字符串看成 base 進制的多項式：
        H(s) = ( ord(s[0])*base^(n-1) + ord(s[1])*base^(n-2) + ... + ord(s[n-1]) ) mod M
    等價的遞推寫法是 `h[i+1] = h[i]*base + ord(s[i])`，預處理 O(n)。
    則任意子串 `s[l:r]` 的指紋為
        h[r] - h[l] * base^(r-l)      (mod M)
    也就是「把前綴左移對齊後相減」，O(1) 得到任意子串指紋。
    這就是 Rabin-Karp 的核心：窗口滑動時指紋可以增量更新，不必重新掃一遍子串。

    ### 為什麼用雙哈希 + 真實比對
    單模哈希有生日碰撞風險：n 個子串時碰撞概率約 n^2 / (2M)。
    這裡同時取兩個質數模數 1e9+7 與 1e9+9（碰撞概率降到 ~1e-18），
    並且在每個哈希命中的位置**再真實比對一次字符**（verify），
    因此輸出結果與暴力法逐字節一致，不存在哈希碰撞導致的誤判。

    ### 五個子問題
      - **子串匹配**：滑動長度 |pat| 的窗口比指紋 → 平均 O(n+m)，最壞與暴力同階 O(n·m) 但常數極小。
      - **最長公共子串**：「存在長度 L 的公共子串」對 L 單調 ⟹ 二分 L，
        判定用哈希集合 O(n+m)，總 O((n+m) log min(n,m))，比 O(n·m) 的二維 DP 省空間。
      - **最長回文子串**：枚舉 2n-1 個回文中心，**對每個中心二分回文半徑**
        （半徑單調：半徑 r 是回文 ⟹ 半徑 r-1 也是），把「正串窗口」與「反串對應窗口」的指紋比對
        → O(n log n)。注意不能直接對「長度」二分：存在長度 4 的回文不代表存在長度 3 的回文
        （例如 `baab`），「存在長度恰為 L 的回文」對 L **不單調**。
      - **不同子串個數**：枚舉長度 L，把 n-L+1 個窗口指紋塞進集合 → O(n^2) 時間。
      - **k-gram Jaccard（應用）**：把文本切成長度 k 的滑動片段集合（shingle），
        相似度 = |A∩B| / |A∪B|。這是抄襲檢測、網頁去重（Broder's shingling）的經典做法，
        也是「把字符串問題轉成集合問題」的典型應用。

應用場景：
    抄襲 / 重複內容檢測（shingling + Jaccard / MinHash）、生物序列比對（DNA 片段 fingerprint）、
    大文件差分同步（rsync 的弱滾動校驗和）、編譯器與 IDE 的增量字符串搜索。

複雜度：
    記 n = len(text)，m = len(pattern)。
      預處理            O(n + m) 時間、O(n + m) 空間
      子串匹配          平均 O(n + m)，最壞 O(n·m)
      最長公共子串      O((n + m) log min(n, m))
      最長回文子串      O(n log n)
      不同子串個數      O(n^2)（小規模驗證用）
      k-gram 相似度     O(n + m)

輸入格式（stdin）：
    第 1 行：text
    第 2 行：pattern
    （兩行都可以是空行，代表空串；不足兩行時缺的部分視為空串；多餘的行忽略）
輸出格式（stdout）：
    第 1 行：pattern 在 text 中的出現次數
    第 2 行：所有起始位置（0-indexed，空格分隔；無則輸出空行）
    第 3 行：最長公共子串長度
    第 4 行：最長公共子串（最靠左的那個；長度為 0 時輸出空行）
    第 5 行：text 的最長回文子串長度
    第 6 行：最長回文子串（最靠左的那個；長度為 0 時輸出空行）
    第 7 行：text 的不同子串個數
    第 8 行：pattern 的不同子串個數
    第 9 行：3-gram Jaccard 相似度，最簡分數 `p/q`
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from math import gcd
from typing import List, Set, Tuple

MOD1 = 1_000_000_007
MOD2 = 1_000_000_009
BASE = 911382323          # 對兩個模數都是奇數且小於模數，保證 base^(r-l) 可逆意義下的良好分佈
K_SHINGLE = 3             # 應用：k-gram 的窗口長度


# ---------------------------------------------------------------- 滾動哈希


class RollingHash:
    """多項式滾動哈希：h[i] 為前 i 個字符的指紋，支援 O(1) 取任意子串指紋。"""

    __slots__ = ("n", "h1", "h2", "p1", "p2")

    def __init__(self, s: str) -> None:
        n = len(s)
        self.n = n
        self.h1 = [0] * (n + 1)
        self.h2 = [0] * (n + 1)
        self.p1 = [1] * (n + 1)
        self.p2 = [1] * (n + 1)
        for i, ch in enumerate(s):
            c = ord(ch)
            self.h1[i + 1] = (self.h1[i] * BASE + c) % MOD1
            self.h2[i + 1] = (self.h2[i] * BASE + c) % MOD2
            self.p1[i + 1] = self.p1[i] * BASE % MOD1
            self.p2[i + 1] = self.p2[i] * BASE % MOD2

    def get(self, l: int, r: int) -> Tuple[int, int]:
        """子串 s[l:r] 的雙指紋。"""
        return ((self.h1[r] - self.h1[l] * self.p1[r - l]) % MOD1,
                (self.h2[r] - self.h2[l] * self.p2[r - l]) % MOD2)


def polynomial_hash(s: str) -> Tuple[int, int]:
    """整串的指紋（用於集合類比對）。"""
    a1 = a2 = 0
    for ch in s:
        a1 = (a1 * BASE + ord(ch)) % MOD1
        a2 = (a2 * BASE + ord(ch)) % MOD2
    return (a1, a2)


# ---------------------------------------------------------------- 子串匹配


def rabin_karp(text: str, pat: str) -> List[int]:
    """Rabin-Karp：返回 pat 在 text 中所有出現的起始位置（升序）。"""
    n, m = len(text), len(pat)
    if m == 0 or m > n:
        return []
    rh_t = RollingHash(text)
    rh_p = RollingHash(pat)
    target = rh_p.get(0, m)
    hits: List[int] = []
    for i in range(n - m + 1):
        if rh_t.get(i, i + m) == target and text[i:i + m] == pat:
            hits.append(i)
    return hits


def kmp_search(text: str, pat: str) -> List[int]:
    """KMP：作為 Rabin-Karp 的獨立對拍實現。"""
    n, m = len(text), len(pat)
    if m == 0 or m > n:
        return []
    pi = [0] * m
    for i in range(1, m):
        j = pi[i - 1]
        while j > 0 and pat[i] != pat[j]:
            j = pi[j - 1]
        if pat[i] == pat[j]:
            j += 1
        pi[i] = j
    hits: List[int] = []
    j = 0
    for i in range(n):
        while j > 0 and text[i] != pat[j]:
            j = pi[j - 1]
        if text[i] == pat[j]:
            j += 1
        if j == m:
            hits.append(i - m + 1)
            j = pi[j - 1]
    return hits


def naive_search(text: str, pat: str) -> List[int]:
    """暴力法：三方對拍的基準。"""
    n, m = len(text), len(pat)
    if m == 0 or m > n:
        return []
    return [i for i in range(n - m + 1) if text[i:i + m] == pat]


# ---------------------------------------------------------------- 最長公共子串


def _has_common_len(a: str, b: str, rha: RollingHash, rhb: RollingHash, length: int) -> bool:
    """是否存在長度為 length 的公共子串。"""
    if length <= 0:
        return True
    if length > len(a) or length > len(b):
        return False
    seen: Set[Tuple[int, int]] = set()
    for j in range(len(b) - length + 1):
        seen.add(rhb.get(j, j + length))
    for i in range(len(a) - length + 1):
        if rha.get(i, i + length) in seen:
            if b.find(a[i:i + length]) >= 0:      # 真實比對，杜絕哈希碰撞誤判
                return True
    return False


def longest_common_substring(a: str, b: str) -> Tuple[int, str]:
    """最長公共子串：二分長度 + 哈希集合，返回 (長度, 最靠左的實例)。"""
    if not a or not b:
        return 0, ""
    rha, rhb = RollingHash(a), RollingHash(b)
    lo, hi = 0, min(len(a), len(b))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if _has_common_len(a, b, rha, rhb, mid):
            lo = mid
        else:
            hi = mid - 1
    if lo == 0:
        return 0, ""
    length = lo
    for i in range(len(a) - length + 1):          # 取最靠左的 i
        if b.find(a[i:i + length]) >= 0:
            return length, a[i:i + length]
    return 0, ""


def lcs_bruteforce(a: str, b: str) -> Tuple[int, str]:
    """最長公共子串暴力版（三方對拍基準）。"""
    best, bi = 0, 0
    for i in range(len(a)):
        for j in range(len(b)):
            k = 0
            while i + k < len(a) and j + k < len(b) and a[i + k] == b[j + k]:
                k += 1
            if k > best:
                best, bi = k, i
    return best, (a[bi:bi + best] if best else "")


# ---------------------------------------------------------------- 最長回文子串


def _is_pal_window(s: str, rhs: RollingHash, rhr: RollingHash, l: int, r: int) -> bool:
    """判斷 s[l:r] 是否為回文：正反哈希相等後再真實比對一次。"""
    n = len(s)
    if l >= r:                                    # 空窗口（偶中心半徑 0）視為回文
        return True
    fwd = rhs.get(l, r)
    rev = rhr.get(n - r, n - l)                   # 反串中對應的同一段（已反轉）
    if fwd != rev:
        return False
    return s[l:r] == s[l:r][::-1]                 # 真實比對，杜絕哈希碰撞誤判


def longest_palindrome_hash(s: str) -> Tuple[int, str]:
    """最長回文子串：枚舉中心 + 二分回文半徑，返回 (長度, 最靠左的實例)。"""
    n = len(s)
    if n == 0:
        return 0, ""
    rhs, rhr = RollingHash(s), RollingHash(s[::-1])
    best_len, best_i = 1, 0

    # 奇數長度：中心是一個字符 c，半徑 r，窗口 [c-r, c+r+1)
    for c in range(n):
        lo, hi = 0, min(c, n - 1 - c)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if _is_pal_window(s, rhs, rhr, c - mid, c + mid + 1):
                lo = mid
            else:
                hi = mid - 1
        length = 2 * lo + 1
        start = c - lo
        if length > best_len or (length == best_len and start < best_i):
            best_len, best_i = length, start

    # 偶數長度：中心在 c-1 與 c 之間，半徑 r >= 1，窗口 [c-r, c+r)
    for c in range(1, n):
        lo, hi = 0, min(c, n - c)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if _is_pal_window(s, rhs, rhr, c - mid, c + mid):
                lo = mid
            else:
                hi = mid - 1
        if lo >= 1:
            length = 2 * lo
            start = c - lo
            if length > best_len or (length == best_len and start < best_i):
                best_len, best_i = length, start

    return best_len, s[best_i:best_i + best_len]


def palindrome_bruteforce(s: str) -> Tuple[int, str]:
    """最長回文子串暴力版（中心擴展，對拍基準）。"""
    n = len(s)
    best, bi = 0, 0
    for c in range(n):
        for a, b in ((c, c), (c, c + 1)):
            i, j = a, b
            while i >= 0 and j < n and s[i] == s[j]:
                i -= 1
                j += 1
            if j - i - 1 > best or (j - i - 1 == best and i + 1 < bi):
                best, bi = j - i - 1, i + 1
    return best, (s[bi:bi + best] if best else "")


# ---------------------------------------------------------------- 不同子串 / 應用


def count_distinct_substrings(s: str) -> int:
    """不同子串個數：枚舉長度 + 指紋集合，O(n^2)。"""
    n = len(s)
    seen: Set[Tuple[int, int]] = set()
    rh = RollingHash(s)
    for length in range(1, n + 1):
        for i in range(n - length + 1):
            seen.add(rh.get(i, i + length))
    return len(seen)


def distinct_substrings_bruteforce(s: str) -> int:
    """不同子串個數暴力版（直接收集真實子串）。"""
    return len({s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)})


def shingles(s: str, k: int = K_SHINGLE) -> Set[str]:
    """k-gram 集合（shingle）；長度不足 k 時整串自成一格，空串為空集。"""
    if not s:
        return set()
    if len(s) <= k:
        return {s}
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def shingle_jaccard(a: str, b: str, k: int = K_SHINGLE) -> str:
    """應用：k-gram Jaccard 相似度，以最簡分數 p/q 返回（避免浮點誤差導致兩語言不一致）。"""
    sa, sb = shingles(a, k), shingles(b, k)
    inter = len(sa & sb)
    union = len(sa | sb)
    if union == 0:
        return "1/1"
    g = gcd(inter, union)
    return f"{inter // g}/{union // g}"


# ---------------------------------------------------------------- IO 與測試


def run_io(raw: str) -> None:
    lines = raw.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    lines = [ln.rstrip("\r") for ln in lines]
    text = lines[0] if len(lines) > 0 else ""
    pat = lines[1] if len(lines) > 1 else ""

    pos = rabin_karp(text, pat)
    clen, csub = longest_common_substring(text, pat)
    plen, psub = longest_palindrome_hash(text)
    out = [
        str(len(pos)),
        " ".join(str(x) for x in pos),
        str(clen),
        csub,
        str(plen),
        psub,
        str(count_distinct_substrings(text)),
        str(count_distinct_substrings(pat)),
        shingle_jaccard(text, pat),
    ]
    sys.stdout.write("\n".join(out) + "\n")


def run_tests() -> None:
    # ---- 固定用例 ----
    assert rabin_karp("abracadabra", "abra") == [0, 7]
    assert rabin_karp("aaaaa", "aa") == [0, 1, 2, 3]
    assert rabin_karp("abc", "") == []
    assert rabin_karp("", "a") == []
    assert rabin_karp("abc", "abcd") == []
    assert rabin_karp("a", "a") == [0]
    assert kmp_search("abracadabra", "abra") == [0, 7]
    assert naive_search("abracadabra", "abra") == [0, 7]

    assert longest_common_substring("abcdef", "zcdemf") == (3, "cde")
    assert longest_common_substring("", "abc") == (0, "")
    assert longest_common_substring("abc", "def") == (0, "")

    assert longest_palindrome_hash("babad") == (3, "bab")
    assert longest_palindrome_hash("cbbd") == (2, "bb")
    assert longest_palindrome_hash("") == (0, "")
    assert longest_palindrome_hash("a") == (1, "a")

    # 整串指紋與前綴哈希口徑一致
    rh = RollingHash("abracadabra")
    assert polynomial_hash("abracadabra") == rh.get(0, 11)
    assert polynomial_hash("") == rh.get(0, 0)

    assert count_distinct_substrings("aaa") == 3
    assert count_distinct_substrings("") == 0
    assert shingle_jaccard("abcde", "abcde") == "1/1"
    assert shingle_jaccard("abcde", "fghij") == "0/1"

    # ---- 隨機對拍：三種匹配實現必須一致 ----
    random.seed(20261007)
    alphabet = "ab"
    for _ in range(400):
        n = random.randint(0, 30)
        m = random.randint(0, 5)
        text = "".join(random.choice(alphabet) for _ in range(n))
        pat = "".join(random.choice(alphabet) for _ in range(m))
        a = rabin_karp(text, pat)
        b = kmp_search(text, pat)
        c = naive_search(text, pat)
        assert a == b == c, (text, pat, a, b, c)

    # ---- 隨機對拍：最長公共子串 / 最長回文子串 / 不同子串 ----
    for _ in range(300):
        n = random.randint(0, 14)
        m = random.randint(0, 14)
        a = "".join(random.choice(alphabet) for _ in range(n))
        b = "".join(random.choice(alphabet) for _ in range(m))
        got = longest_common_substring(a, b)
        want = lcs_bruteforce(a, b)
        assert got[0] == want[0], (a, b, got, want)
        if want[0] > 0:
            assert got[1] == want[1] or (len(got[1]) == want[0] and got[1] in b)

        gp = longest_palindrome_hash(a)
        wp = palindrome_bruteforce(a)
        assert gp[0] == wp[0], (a, gp, wp)
        assert gp[1] == gp[1][::-1] and gp[1] in a
        assert len(gp[1]) == gp[0]

        assert count_distinct_substrings(a) == distinct_substrings_bruteforce(a)
        assert count_distinct_substrings(b) == distinct_substrings_bruteforce(b)

    # ---- 應用：k-gram Jaccard 與直接集合計算對拍 ----
    for _ in range(200):
        n = random.randint(0, 20)
        m = random.randint(0, 20)
        a = "".join(random.choice("abc") for _ in range(n))
        b = "".join(random.choice("abc") for _ in range(m))
        sa, sb = shingles(a), shingles(b)
        inter = len(sa & sb)
        union = len(sa | sb)
        if union == 0:
            assert shingle_jaccard(a, b) == "1/1"
        else:
            g = gcd(inter, union)
            assert shingle_jaccard(a, b) == f"{inter // g}/{union // g}"

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
