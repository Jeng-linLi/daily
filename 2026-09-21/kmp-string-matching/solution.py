"""KMP 字符串匹配（Knuth-Morris-Pratt）

在 text 中查找 pattern 首次出現的位置；若不存在返回 -1。
核心：先對 pattern 求前綴函數（部分匹配表），失配時利用已匹配信息跳過不可能的位置，
從而把樸素算法的 O(n*m) 降到 O(n+m)。
"""

from typing import List


def build_prefix(pattern: str) -> List[int]:
    """構造前綴函數 pi：pi[i] = pattern[0..i] 的最長真前綴同時也是後綴的長度。"""
    pi = [0] * len(pattern)
    j = 0  # 當前已匹配的前綴長度
    for i in range(1, len(pattern)):
        # 失配時回退到更短的前綴，直到能接上或回到 0
        while j > 0 and pattern[i] != pattern[j]:
            j = pi[j - 1]
        if pattern[i] == pattern[j]:
            j += 1
        pi[i] = j
    return pi


def kmp_search(text: str, pattern: str) -> int:
    """返回 pattern 在 text 中首次出現的下標，不存在返回 -1。空 pattern 視爲出現在 0。"""
    if not pattern:
        return 0
    pi = build_prefix(pattern)
    j = 0  # 當前在 pattern 上匹配到的長度
    for i, ch in enumerate(text):
        while j > 0 and ch != pattern[j]:
            j = pi[j - 1]  # 關鍵：不回退 text 指針，只回退 pattern 指針
        if ch == pattern[j]:
            j += 1
        if j == len(pattern):
            return i - len(pattern) + 1
    return -1


def kmp_search_all(text: str, pattern: str) -> List[int]:
    """返回 pattern 在 text 中所有出現位置（允許重疊）。"""
    if not pattern:
        return list(range(len(text) + 1))
    pi = build_prefix(pattern)
    hits, j = [], 0
    for i, ch in enumerate(text):
        while j > 0 and ch != pattern[j]:
            j = pi[j - 1]
        if ch == pattern[j]:
            j += 1
        if j == len(pattern):
            hits.append(i - len(pattern) + 1)
            j = pi[j - 1]  # 繼續尋找下一次匹配
    return hits


if __name__ == "__main__":
    assert kmp_search("ababcabcabababd", "ababd") == 10
    assert kmp_search("hello", "ll") == 2
    assert kmp_search("aaaaa", "bba") == -1
    assert kmp_search("abc", "") == 0
    assert kmp_search_all("ababab", "aba") == [0, 2]
    assert build_prefix("ababaca") == [0, 0, 1, 2, 3, 0, 1]
    print("all tests passed")
