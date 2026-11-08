# -*- coding: utf-8 -*-
"""回文自動機 Palindromic Tree（Eertree）

題意：給定字串 s，統計其所有回文子字串：
  1) 相異回文子字串個數 k（經典結論：≤ |s|）；
  2) 每個相異回文子字串的出現次數（按出現位置計，可重疊）；
  3) 每個相異回文子字串首次出現的起始下標。

回文自動機性質：
  - 每個節點代表一個相異回文子字串；另有兩個根：
    節點 0 是「虛根」（len = -1），節點 1 是「空串根」（len = 0）。
  - 每次加入一個字元，最多新增一個相異回文（以該位置結尾的最長回文），
    其餘以該位置結尾的回文都是它的真後綴、由 fail（後綴連結）鏈給出。
  - 出現次數：先在建立節點時計 +1，再沿 fail 樹自葉向根累加，
    節點 v 的 cnt 即該回文在整串中的出現次數。

IO 模式（stdin 有內容時）：
  一行字串 s（僅小寫字母）
  輸出：
    第一行：k（相異回文子字串個數）
    第二行：total（回文子字串總出現次數，可重疊）
    接下來 k 行，按 (長度, 首次起始下標) 升冪：
      len start count
無 stdin 時跑內建斷言（與 O(n²) 暴力枚舉全部子字串對拍）並輸出 all tests passed。
"""
import sys


def eertree(s):
    """建立回文自動機，回傳 [(len, first_start, count), ...] 按 (len, start) 升冪。"""
    n = len(s)
    maxn = n + 3
    length = [0] * maxn          # 節點代表的回文長度
    link = [0] * maxn            # fail：最長真回文後綴
    cnt = [0] * maxn             # 出現次數（建好後沿 fail 樹累加）
    nxt = [dict() for _ in range(maxn)]  # 轉移邊
    endpos = [0] * maxn          # 節點建立時的結束下標（= 首次出現的結尾）

    length[0], link[0] = -1, 0   # 虛根
    length[1], link[1] = 0, 0    # 空串根
    last, size = 1, 2

    def get_link(v, i):
        # 沿 fail 鏈找「兩側都能擴展出 s[i]」的節點
        while True:
            l = length[v]
            if i - l - 1 >= 0 and s[i - l - 1] == s[i]:
                return v
            v = link[v]

    for i, ch in enumerate(s):
        cur = get_link(last, i)
        if ch in nxt[cur]:
            # 該回文已存在
            last = nxt[cur][ch]
            cnt[last] += 1
        else:
            now = size
            size += 1
            length[now] = length[cur] + 2
            endpos[now] = i
            cnt[now] = 1
            nxt[cur][ch] = now
            if length[now] == 1:
                # 長度 1 的回文，fail 指向空串根
                # （若走 get_link(link[cur]) 會指回自己，必須特判）
                link[now] = 1
            else:
                link[now] = nxt[get_link(link[cur], i)][ch]
            last = now

    # 沿 fail 樹自底向上累加出現次數
    for v in range(size - 1, 1, -1):
        cnt[link[v]] += cnt[v]

    nodes = [(length[v], endpos[v] - length[v] + 1, cnt[v]) for v in range(2, size)]
    nodes.sort()  # (len, first_start) 唯一確定一個相異回文
    return nodes


def brute(s):
    """O(n²·n) 暴力：枚舉全部子字串判回文，統計 (len, first_start, count)。"""
    n = len(s)
    info = {}  # 子字串 -> [count, first_start]
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub = s[i:j]
            if sub == sub[::-1]:
                if sub not in info:
                    info[sub] = [0, i]
                info[sub][0] += 1
    return sorted((len(p), st, c) for p, (c, st) in info.items())


def run_tests():
    fixed = ["", "a", "aa", "ab", "aba", "aaaa", "abacaba", "abcba", "aabaa",
             "banana", "abcabc", "abab", "racecar", "aaaaaa", "abaaba"]
    for s in fixed:
        assert eertree(s) == brute(s), "fixed case failed: %r" % s

    # 自製 LCG 隨機對拍（兩語言一致；此處只與暴力比對）
    seed = 20261108
    mod = 1 << 31

    def rnd(n):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) % mod
        return seed % n

    for _ in range(300):
        n = rnd(24)
        s = "".join(chr(ord('a') + rnd(3)) for _ in range(n))
        assert eertree(s) == brute(s), "random case failed: %r" % s

    # 性質檢查：相異回文數 ≤ |s|；總出現次數 = Σcnt
    for s in fixed:
        nodes = eertree(s)
        assert len(nodes) <= max(len(s), 1)
        assert sum(c for _, _, c in nodes) == sum(1 for i in range(len(s))
                                                  for j in range(i + 1, len(s) + 1)
                                                  if s[i:j] == s[i:j][::-1])


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        print("all tests passed")
        return
    s = data[0].decode("ascii")
    nodes = eertree(s)
    out = []
    out.append(str(len(nodes)))
    out.append(str(sum(c for _, _, c in nodes)))
    for ln, st, c in nodes:
        out.append("%d %d %d" % (ln, st, c))
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
