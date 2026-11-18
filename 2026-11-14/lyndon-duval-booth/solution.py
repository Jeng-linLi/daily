#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lyndon 分解（Duval 演算法）與循環串最小／最大表示法（Booth 演算法）。

命令（IO 模式，從 stdin 讀入；S 為非空 token）：
  lyndon    S   -> Lyndon 分解，各因子以單個空格分隔
  minrot    S   -> 最小表示的起始下標（0-based）與該循環串
  maxrot    S   -> 最大表示的起始下標（0-based）與該循環串

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * Lyndon 詞：嚴格小於它所有**真後綴**的非空串。任何串都能唯一分解為
    c1 c2 ... ck，每個 ci 都是 Lyndon 詞，且 c1 >= c2 >= ... >= ck（字典序）。
  * Duval：從 i 開始延伸「週期性主鏈」——比較 s[k] 與 s[j]，相等則 k++（週期不變），
    s[k] < s[j] 則 k 回到 i（更有希望的新週期），s[k] > s[j] 則主鏈斷開。
    主鏈長度 j-i、週期 j-k，於是從 i 起連續切出 (j-k) 長度的因子直到超過 k。
    全程 i、j 單調前進，O(n)。
  * Booth：把 s 複製成 s+s，維護兩個候選起點 i、j 與已匹配長度 k。
    s[i+k] > s[j+k] 時 i 這一段全部不可能（i 到 i+k 都輸給對應的 j 段），
    直接跳到 i+k+1；反之跳 j。k 從不倒退，O(n)。
  * 最大表示只是把 Booth 的比較方向反過來（s[i+k] < s[j+k] 時淘汰 i）。
"""

import sys


def duval(s):
    """Lyndon 分解：回傳因子列表（非空串）。"""
    n = len(s)
    res = []
    i = 0
    while i < n:
        j, k = i + 1, i          # k 是主鏈的週期對應位置
        while j < n and s[k] <= s[j]:
            if s[k] < s[j]:
                k = i            # 出現更大的字元：週期重設
            else:
                k += 1           # 相等：週期延續
            j += 1
        # 主鏈 break 掉，長度為 j-i 的週期節，可切出若干個長 (j-k) 的因子
        while i <= k:
            res.append(s[i:i + j - k])
            i += j - k
    return res


def booth(s, want_max=False):
    """回傳循環串最小（或最大）表示的起始下標。空串回傳 0。"""
    n = len(s)
    if n == 0:
        return 0
    ss = s + s
    i, j, k = 0, 1, 0
    while i < n and j < n and k < n:
        a, b = ss[i + k], ss[j + k]
        if a == b:
            k += 1
            continue
        # want_max 時「a 較小」才淘汰 i，其餘與最小版對稱
        if (a < b) if want_max else (a > b):
            i += k + 1
            if i == j:
                i += 1
        else:
            j += k + 1
            if i == j:
                j += 1
        k = 0
    return min(i, j)


def rotation(s, start):
    """取循環串從 start 開始、長度 len(s) 的那一段。"""
    n = len(s)
    if n == 0:
        return ''
    return s[start:] + s[:start] if start < n else s


# ----------------------------- 暴力參考實作 -----------------------------

def is_lyndon(t):
    """t 嚴格小於所有真後綴。"""
    return all(t < t[p:] for p in range(1, len(t)))


def brute_decompose(s):
    """貪心：每次取最長的 Lyndon 前綴（O(n^3)，僅供對拍）。"""
    res = []
    i = 0
    n = len(s)
    while i < n:
        best = -1
        for j in range(n, i, -1):
            if is_lyndon(s[i:j]):
                best = j
                break
        if best < 0:      # 單個字元一定是 Lyndon 詞，正常不會走到
            best = i + 1
        res.append(s[i:best])
        i = best
    return res


def brute_rotation(s, want_max=False):
    """枚舉所有循環移位取字典序最小（最大）者，平手取最小下標。"""
    n = len(s)
    if n == 0:
        return 0, ''
    best_i, best = 0, None
    for p in range(n):
        cur = s[p:] + s[:p]
        if best is None or (cur > best if want_max else cur < best):
            best_i, best = p, cur
    return best_i, best


# ------------------------------- IO 模式 -------------------------------

def run_io(data):
    toks = data.split()
    pos = 0
    out = []

    def nxt():
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    while True:
        cmd = nxt()
        if cmd is None:
            break
        if cmd not in ('lyndon', 'minrot', 'maxrot'):
            break
        s = nxt()
        if s is None:
            break
        if cmd == 'lyndon':
            out.append(' '.join(duval(s)))
        elif cmd == 'minrot':
            p = booth(s, False)
            out.append('%d %s' % (p, rotation(s, p)))
        else:
            p = booth(s, True)
            out.append('%d %s' % (p, rotation(s, p)))
    return ''.join(line + '\n' for line in out)


# ------------------------------ 內建測試 ------------------------------

def run_tests():
    # 1) 固定用例
    assert duval('abcde') == ['abcde']
    assert duval('banana') == ['b', 'an', 'an', 'a']
    assert duval('aaaa') == ['a', 'a', 'a', 'a']
    assert duval('ababb') == ['ababb']
    assert booth('banana') == 5                 # "abanan"
    assert rotation('banana', 5) == 'abanan'
    assert booth('cba') == 2                    # "acb"
    assert booth('a') == 0
    assert booth('banana', True) == 2           # "nanaba"（n 之後 b < n 輸給下一輪）
    assert booth('abab') == 0                   # 平手取最小下標
    assert duval('cabca') == ['c', 'abc', 'a']

    # 2) 隨機對拍：Duval vs 貪心暴力；Booth vs 全枚舉
    state = 20261114

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    for _ in range(400):
        n = rnd(24) + 1
        alpha = rnd(3) + 1
        s = ''.join(chr(ord('a') + rnd(alpha)) for _ in range(n))
        got = duval(s)
        want = brute_decompose(s)
        assert got == want, (s, got, want)
        # Lyndon 分解的性質：每項是 Lyndon 詞且非遞增
        assert all(is_lyndon(t) for t in got), (s, got)
        assert all(got[t] >= got[t + 1] for t in range(len(got) - 1)), (s, got)
        assert ''.join(got) == s
        pi, ps = brute_rotation(s, False)
        assert booth(s, False) == pi, (s, booth(s, False), pi)
        assert rotation(s, booth(s, False)) == ps
        qi, qs = brute_rotation(s, True)
        assert booth(s, True) == qi, (s, booth(s, True), qi)
        assert rotation(s, booth(s, True)) == qs

    # 3) 週期串與長串（考驗 k 不倒退的前提）
    for rep in range(1, 40):
        s = 'ab' * rep
        assert booth(s) == 0
        assert duval(s) == ['ab'] * rep
    big = 'z' * 5000 + 'a' + 'z' * 4999
    assert booth(big) == 5000          # 最小表示從那個 'a' 開始
    assert duval('a' * 3000) == ['a'] * 3000


def main():
    data = sys.stdin.read()
    if data.strip():
        # 二進位輸出：避免 Windows 文字模式把 \n 轉成 \r\n（與 C++ 版不一致）
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
