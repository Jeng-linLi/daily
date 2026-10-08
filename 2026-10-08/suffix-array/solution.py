"""後綴陣列（Suffix Array）：倍增構造 + Kasai LCP + 最長重複子串 / 不同子串數 / BWT / 子串搜尋

題意：
    給定字串 s（長度 n）與一個模式串 pat，求：
      1. **後綴陣列 sa**：把 s 的 n 個後綴按字典序排序後，記錄每個後綴的起始下標；
      2. **名次陣列 rank**：rank[i] = 後綴 s[i:] 在 sa 中的排名（sa 的逆陣列）；
      3. **LCP 陣列**：lcp[i] = 排名相鄰的兩個後綴 sa[i] 與 sa[i-1] 的最長公共前綴長度（lcp[0] = 0）；
      4. **最長重複子串**：在 s 中出現至少兩次的最長子串（長度、內容、起始下標）；
      5. **不同子串數**：s 中本質不同的子串共有幾個；
      6. **BWT 變換**：Burrows-Wheeler Transform，bzip2 的核心；
      7. **子串搜尋**：pat 在 s 中出現幾次、分別在哪些位置（升序）。

思路：
    ### 為什麼要後綴陣列
    後綴樹的空間與常數都很大，後綴陣列用「一個整數陣列」達到同樣的效果：
    所有後綴**排序**之後，具有相同前綴的後綴會聚成連續區間，
    於是大量子串問題都能用二分搜尋 + LCP 在 O(|pat| log n) 內解決。

    ### 倍增構造 O(n log n)
    直接對 n 個後綴做字串排序是 O(n² log n)。倍增法每輪把「已排好序的長度 k 前綴」
    組合成「長度 2k 前綴」的排序鍵：
        key(i) = ( rank[i], rank[i + k] if i + k < n else -1 )
    兩個長度 k 的排名拼起來就是長度 2k 的排名；越界的補 -1（空串字典序最小）。
    每輪把 k 加倍，最多 log n 輪；當所有排名兩兩不同（rank[sa[n-1]] == n-1）即可提前結束。
    用比較排序是 O(n log² n)，用計數排序（基數排序）可降到 O(n log n)；
    兩種寫法最終順序**完全一致**——因為最後一輪的排序鍵兩兩不同，比較排序的結果唯一。

    ### Kasai 求 LCP，O(n)
    設 h[i] = 後綴 s[i:] 與「它在 sa 中的前一名」的 LCP 長度。
    關鍵性質：h[i+1] ≥ h[i] - 1（把兩個後綴都砍掉首字元，公共前綴最多少 1）。
    因此只需從 h[i] - 1 開始繼續比對，指標總位移是 O(n)，整體線性。

    ### 最長重複子串 = max(lcp)
    任何重複出現的子串都是某兩個後綴的公共前綴，而相鄰排名的 LCP 是全域最大的候選
    （因為 LCP(sa[i], sa[j]) = min(lcp[i+1..j])），所以答案就是 max(lcp)。

    ### 不同子串數 = n(n+1)/2 - Σ lcp
    全部子串有 n(n+1)/2 個。按排名看，後綴 sa[i] 會貢獻「長度 1..n-sa[i]」這些子串，
    其中前 lcp[i] 個已經被前一名貢獻過了，扣掉 Σ lcp 即為本質不同的個數。

    ### BWT
    bwt[i] = s[(sa[i] - 1 + n) % n]，也就是每個後綴往前挪一格的那一撇字元。
    它把重複子串集中成連續的相同字元，配合 MTF + 霍夫曼編碼就是 bzip2。

    ### 子串搜尋 = 兩次二分
    所有以 pat 開頭的後綴在 sa 中必定是一段**連續區間**。
    用兩次二分找下界（第一個後綴 ≥ pat）與上界（第一個後綴 > pat），
    區間長度就是出現次數 → **O(|pat| log n)**，比 KMP 的 O(n + |pat|) 慢一點，
    但同一份 SA 可以服務成千上萬次查詢（這是後綴陣列真正的優勢）。

    ### 確定性（保證 Python 與 C++ 輸出逐字節一致）
      - 排序鍵相同 → 排序結果唯一；
      - 最長重複子串平手時取 sa 序最前者（第一次出現較早者）；
      - 出現位置升序輸出；BWT 按 sa 順序輸出。

應用場景：
    基因序列比對（BLAST 的核心索引）、全文檢索與搜尋引擎的倒排索引、
    最長重複片段檢測（抄襲 / 版權查重）、bzip2 壓縮（BWT）、
    字串壓縮的 LZ77 匹配尋找、生物資訊的 reads 比對。

複雜度：
    倍增構造       O(n log n)（本實作用比較排序，為 O(n log² n)）、O(n) 空間
    Kasai LCP      O(n) 時間、O(n) 空間
    最長重複子串   O(n)
    不同子串數     O(n)
    BWT            O(n)
    子串搜尋       O(|pat| log n)
    暴力對拍       O(n² log n)（僅測試用，n ≤ 12）

輸入格式（stdin，全部以空白分隔）：
    s       第一個 token 為待處理字串（缺則視為空串）
    pat     第二個 token 為模式串（缺則視為空串）
輸出格式（stdout）：
    第 1 行：後綴陣列 sa（空格分隔）
    第 2 行：名次陣列 rank（空格分隔）
    第 3 行：LCP 陣列（空格分隔）
    第 4 行：最長重複子串長度
    第 5 行：最長重複子串（不存在時為空行）
    第 6 行：最長重複子串的起始下標（-1 表示不存在）
    第 7 行：不同子串數
    第 8 行：BWT 變換結果（空串時為空行）
    第 9 行：pat 的出現次數
    第 10 行：pat 的出現位置（升序，空格分隔；無則為空行）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import List, Tuple

NEG = -1                    # 越界的排名鍵（-1 比任何合法排名都小，對應「空串最小」）


# ---------------------------------------------------------------- 構造與 LCP


def build_sa(s: str) -> List[int]:
    """倍增法構造後綴陣列，返回 sa（把 n 個後綴按字典序排序後的起始下標）。"""
    n = len(s)
    if n == 0:
        return []
    sa = list(range(n))
    rank = [ord(ch) for ch in s]
    k = 1
    while k < n:
        def key(i: int, rnk=rank, kk=k) -> Tuple[int, int]:
            r2 = rnk[i + kk] if i + kk < n else NEG
            return (rnk[i], r2)

        sa.sort(key=key)
        new_rank = [0] * n
        new_rank[sa[0]] = 0
        for i in range(1, n):
            new_rank[sa[i]] = new_rank[sa[i - 1]] + (0 if key(sa[i]) == key(sa[i - 1]) else 1)
        rank = new_rank
        if rank[sa[-1]] == n - 1:        # 排名已兩兩不同，排序完成
            break
        k <<= 1
    return sa


def sa_bruteforce(s: str) -> List[int]:
    """暴力構造：直接對所有後綴做字串排序，O(n² log n)，僅供對拍。"""
    return sorted(range(len(s)), key=lambda i: s[i:])


def build_rank(sa: List[int]) -> List[int]:
    """由 sa 得到名次陣列 rank（sa 的逆置換）。"""
    rank = [0] * len(sa)
    for i, p in enumerate(sa):
        rank[p] = i
    return rank


def lcp_naive(s: str, i: int, j: int) -> int:
    """直接求 s[i:] 與 s[j:] 的最長公共前綴長度，O(n)，僅供對拍。"""
    n = len(s)
    h = 0
    while i + h < n and j + h < n and s[i + h] == s[j + h]:
        h += 1
    return h


def kasai_lcp(s: str, sa: List[int]) -> List[int]:
    """Kasai 演算法，O(n)：lcp[i] = LCP(sa[i], sa[i-1])，lcp[0] = 0。"""
    n = len(s)
    lcp = [0] * n
    if n == 0:
        return lcp
    rank = build_rank(sa)
    h = 0
    for i in range(n):
        r = rank[i]
        if r == 0:
            h = 0
            continue
        j = sa[r - 1]
        while i + h < n and j + h < n and s[i + h] == s[j + h]:
            h += 1
        lcp[r] = h
        if h > 0:
            h -= 1                       # 下一個後綴最多只會少 1，這是線性的關鍵
    return lcp


# ---------------------------------------------------------------- 應用


def longest_repeated_substring(s: str, sa: List[int], lcp: List[int]) -> Tuple[int, int, str]:
    """最長重複子串，返回 (長度, 起始下標, 子串)；不存在則 (0, -1, "")。

    平手時取 sa 序最前者（即第一次出現位置在 sa 中較靠前）。
    """
    best = 0
    start = -1
    for i in range(len(lcp)):
        if lcp[i] > best:                # 嚴格更優才更新 → 平手取最前，保證確定性
            best = lcp[i]
            start = sa[i]
    return best, start, s[start:start + best] if best > 0 else ""


def longest_repeated_bruteforce(s: str) -> Tuple[int, int]:
    """暴力求最長重複子串（n ≤ 12）：枚舉所有 (i, j) 求 LCP，平手取起始下標較小者。"""
    n = len(s)
    best = 0
    start = -1
    for i in range(n):
        for j in range(i + 1, n):
            h = lcp_naive(s, i, j)
            if h == 0:
                continue
            si = min(i, j)
            if h > best or (h == best and (start == -1 or si < start)):
                best = h
                start = si
    return best, start


def count_distinct_substrings(s: str, lcp: List[int]) -> int:
    """本質不同的子串數 = n(n+1)/2 - Σ lcp。"""
    n = len(s)
    return n * (n + 1) // 2 - sum(lcp)


def count_distinct_bruteforce(s: str) -> int:
    """暴力統計不同子串數：直接建集合，僅供對拍。"""
    seen = set()
    for i in range(len(s)):
        for j in range(i + 1, len(s) + 1):
            seen.add(s[i:j])
    return len(seen)


def burrows_wheeler(s: str, sa: List[int]) -> str:
    """BWT：bwt[i] = s[(sa[i] - 1 + n) % n]，把重複子串聚成連續相同字元。"""
    n = len(s)
    if n == 0:
        return ""
    return "".join(s[(p - 1) % n] for p in sa)


def sa_search(s: str, sa: List[int], pat: str) -> Tuple[int, List[int]]:
    """在 sa 上二分搜尋 pat，返回 (出現次數, 出現位置升序列表)。

    所有以 pat 開頭的後綴在 sa 中是一段連續區間，兩次二分即可定位。
    """
    n = len(s)
    if not pat:
        return 0, []
    L = len(pat)

    # 下界：第一個「前綴 ≥ pat」的後綴
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        p = sa[mid]
        if s[p:p + L] < pat:
            lo = mid + 1
        else:
            hi = mid
    left = lo

    # 上界：第一個「前綴 > pat」的後綴（s 的尾端截斷視為較小，故 == pat 時往右走）
    lo, hi = left, n
    while lo < hi:
        mid = (lo + hi) // 2
        p = sa[mid]
        if s[p:p + L] > pat:
            hi = mid
        else:
            lo = mid + 1
    right = lo

    occ = sa[left:right]
    return right - left, sorted(occ)


def search_bruteforce(s: str, pat: str) -> Tuple[int, List[int]]:
    """暴力搜尋所有出現位置，僅供對拍。"""
    if not pat:
        return 0, []
    pos = []
    start = 0
    while True:
        idx = s.find(pat, start)
        if idx < 0:
            break
        pos.append(idx)
        start = idx + 1
    return len(pos), pos


# ---------------------------------------------------------------- IO 與測試


def run_io(raw: str) -> None:
    toks = raw.split()
    s = toks[0] if len(toks) > 0 else ""
    pat = toks[1] if len(toks) > 1 else ""

    sa = build_sa(s)
    rank = build_rank(sa)
    lcp = kasai_lcp(s, sa)
    blen, bstart, bsub = longest_repeated_substring(s, sa, lcp)
    distinct = count_distinct_substrings(s, lcp)
    bwt = burrows_wheeler(s, sa)
    cnt, positions = sa_search(s, sa, pat)

    out = [
        " ".join(str(x) for x in sa),
        " ".join(str(x) for x in rank),
        " ".join(str(x) for x in lcp),
        str(blen),
        bsub,
        str(bstart),
        str(distinct),
        bwt,
        str(cnt),
        " ".join(str(x) for x in positions),
    ]
    sys.stdout.write("\n".join(out) + "\n")


def run_tests() -> None:
    # ---- 空串 ----
    assert build_sa("") == []
    assert kasai_lcp("", []) == []
    assert longest_repeated_substring("", [], []) == (0, -1, "")
    assert count_distinct_substrings("", []) == 0
    assert burrows_wheeler("", []) == ""
    assert sa_search("", [], "a") == (0, [])

    # ---- 單字元 ----
    assert build_sa("a") == [0]
    assert kasai_lcp("a", [0]) == [0]
    assert longest_repeated_substring("a", [0], [0]) == (0, -1, "")
    assert count_distinct_substrings("a", [0]) == 1
    assert burrows_wheeler("a", [0]) == "a"
    assert sa_search("a", [0], "a") == (1, [0])
    assert sa_search("a", [0], "b") == (0, [])

    # ---- 經典：banana ----
    b = "banana"
    sa_b = build_sa(b)
    assert sa_b == [5, 3, 1, 0, 4, 2]
    lcp_b = kasai_lcp(b, sa_b)
    assert lcp_b == [0, 1, 3, 0, 0, 2]
    assert build_rank(sa_b) == [3, 2, 5, 1, 4, 0]
    assert longest_repeated_substring(b, sa_b, lcp_b) == (3, 1, "ana")
    assert count_distinct_substrings(b, lcp_b) == count_distinct_bruteforce(b)
    assert sa_search(b, sa_b, "ana") == (2, [1, 3])
    assert sa_search(b, sa_b, "na") == (2, [2, 4])
    assert sa_search(b, sa_b, "banana") == (1, [0])
    assert sa_search(b, sa_b, "bananas") == (0, [])

    # ---- 經典：mississippi ----
    m = "mississippi"
    sa_m = build_sa(m)
    lcp_m = kasai_lcp(m, sa_m)
    assert longest_repeated_substring(m, sa_m, lcp_m) == (4, 1, "issi")
    assert count_distinct_substrings(m, lcp_m) == 53
    assert count_distinct_substrings(m, lcp_m) == count_distinct_bruteforce(m)
    assert sa_search(m, sa_m, "ssi") == (2, [2, 5])
    assert sa_search(m, sa_m, "issi") == (2, [1, 4])

    # ---- 全同字元：aaaa ----
    a4 = "aaaa"
    sa_a = build_sa(a4)
    lcp_a = kasai_lcp(a4, sa_a)
    assert lcp_a == [0, 1, 2, 3]
    assert longest_repeated_substring(a4, sa_a, lcp_a) == (3, 0, "aaa")
    assert count_distinct_substrings(a4, lcp_a) == 4
    assert burrows_wheeler(a4, sa_a) == "aaaa"
    assert sa_search(a4, sa_a, "aa") == (3, [0, 1, 2])

    # ---- BWT 經典：abracadabra ----
    ab = "abracadabra"
    sa_ab = build_sa(ab)
    assert sa_ab == [10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2]
    lcp_ab = kasai_lcp(ab, sa_ab)
    assert lcp_ab == [0, 1, 4, 1, 1, 0, 3, 0, 0, 0, 2]
    assert burrows_wheeler(ab, sa_ab) == "rdarcaaaabb"
    assert longest_repeated_substring(ab, sa_ab, lcp_ab) == (4, 0, "abra")
    assert count_distinct_substrings(ab, lcp_ab) == 54

    # ---- 隨機對拍 ----
    rnd = random.Random(20261008)
    alphabet = "ab"
    for _ in range(200):
        n = rnd.randint(0, 12)
        s = "".join(rnd.choice(alphabet) for _ in range(n))
        sa = build_sa(s)
        assert sa == sa_bruteforce(s)
        assert len(sa) == n and sorted(sa) == list(range(n))
        rank = build_rank(sa)
        for i in range(n):
            assert rank[sa[i]] == i
        lcp = kasai_lcp(s, sa)
        for i in range(1, n):
            assert lcp[i] == lcp_naive(s, sa[i], sa[i - 1])
        bl, bs, bsub = longest_repeated_substring(s, sa, lcp)
        bl2, bs2 = longest_repeated_bruteforce(s)
        assert bl == bl2, (s, bl, bl2)
        if bl > 0:
            assert s[bs:bs + bl] == bsub
            # 至少還有一個「不同起點」的出現（注意 str.count 只算不重疊的，不能用）
            occ = [j for j in range(n - bl + 1) if s[j:j + bl] == bsub]
            assert bs in occ and len(occ) >= 2
        assert count_distinct_substrings(s, lcp) == count_distinct_bruteforce(s)
        bwt = burrows_wheeler(s, sa)
        assert len(bwt) == n
        assert sorted(bwt) == sorted(s)                     # BWT 是 s 的一個排列
        for plen in range(1, 4):
            pat = "".join(rnd.choice(alphabet) for _ in range(plen))
            c1, p1 = sa_search(s, sa, pat)
            c2, p2 = search_bruteforce(s, pat)
            assert (c1, p1) == (c2, p2), (s, pat, c1, p1, c2, p2)

    # ---- 三元字母表，提高多樣性 ----
    for _ in range(100):
        n = rnd.randint(0, 14)
        s = "".join(rnd.choice("abc") for _ in range(n))
        sa = build_sa(s)
        assert sa == sa_bruteforce(s)
        lcp = kasai_lcp(s, sa)
        for i in range(1, n):
            assert lcp[i] == lcp_naive(s, sa[i], sa[i - 1])
        assert count_distinct_substrings(s, lcp) == count_distinct_bruteforce(s)
        bl, bs, bsub = longest_repeated_substring(s, sa, lcp)
        bl2, _ = longest_repeated_bruteforce(s)
        assert bl == bl2


if __name__ == "__main__":
    data = sys.stdin.read()
    if data.strip():
        run_io(data)
    else:
        run_tests()
        print("all tests passed")
