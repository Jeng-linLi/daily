# -*- coding: utf-8 -*-
"""同餘最短路（Congruence Shortest Path）

題意：給定正整數 b_1..b_n 與多個查詢 x，求最小的 y >= x 使得
      y 可以表示成 b_i 的非負整係數和（y = a_1*b_1 + ... + a_n*b_n, a_i >= 0）。
      若不存在（該同餘類永遠湊不出）輸出 -1。

思路：
  記 m = min(b_i)。把每個非負整數按 mod m 分類，建圖：
    節點為剩餘類 0..m-1，邊 r -> (r + b_i) mod m，邊權 +b_i。
  從 0 跑 Dijkstra，dist[r] 即「湊得出且 ≡ r (mod m) 的最小值」。
  注意：m 本身是某個 b_i，故同餘類內的可湊值恰為 {dist[r] + t*m}；
  但查詢 y >= x 不要求 y 與 x 同餘，所以答案要對全部 m 個剩餘類取 min：
    對每個 dist[r] < INF 的類，該類中首個 >= x 的可湊值為
      dist[r]                       （dist[r] >= x）
      dist[r] + ceil((x-dist[r])/m)*m（否則）
  複雜度 O(n*m*log m + q*m)。

IO 模式（stdin 有內容時）：
  n
  b_1 ... b_n            （正整數，n <= 20，m = min <= 2000）
  q
  x_1 ... x_q            （q <= 200，x <= 1e9）
  輸出 q 行，每行一個答案。
無 stdin 時跑內建斷言並輸出 all tests passed。
"""
import sys
import heapq

INF = 1 << 62


def dijkstra_residues(bs, m):
    """回傳 dist[r]：可湊出且 ≡ r (mod m) 的最小值。"""
    dist = [INF] * m
    dist[0] = 0
    pq = [(0, 0)]
    while pq:
        d, r = heapq.heappop(pq)
        if d > dist[r]:
            continue
        for b in bs:
            nr = (r + b) % m
            nd = d + b
            if nd < dist[nr]:
                dist[nr] = nd
                heapq.heappush(pq, (nd, nr))
    return dist


def solve_query(x, dist, m):
    """最小的可湊 y >= x：對每個剩餘類取該類中首個 >= x 的可湊值，再取 min。"""
    best = INF
    for r in range(m):
        f = dist[r]
        if f >= INF:
            continue
        if f >= x:
            c = f
        else:
            t = -((f - x) // m)  # ceil((x - f) / m) >= 1
            c = f + t * m
        if c < best:
            best = c
    return -1 if best >= INF else best


def brute(bs, m, lim):
    """暴力：完全背包式可達性 DP，lim 內逐一標記。"""
    reach = [False] * (lim + 1)
    reach[0] = True
    for y in range(1, lim + 1):
        for b in bs:
            if y >= b and reach[y - b]:
                reach[y] = True
                break
    bf = [INF] * m
    for y in range(lim + 1):
        if reach[y]:
            r = y % m
            if y < bf[r]:
                bf[r] = y
    return bf


def run_tests():
    seed = 20261106

    def rnd(mv):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF
        return seed % mv

    # 固定用例
    fixed = [[3, 5], [4, 6], [7], [6, 10, 15], [5, 7, 11]]
    for bs in fixed:
        bs = sorted(bs)
        m = bs[0]
        dist = dijkstra_residues(bs, m)
        lim = m * max(bs)
        bf = brute(bs, m, lim)
        assert dist == bf, (bs, dist, bf)
        for x in range(0, lim + 1):
            # 暴力查詢：從 x 往上找到第一個可湊值
            want = -1
            for y in range(x, lim + 1):
                r = y % m
                if bf[r] != INF and y >= bf[r]:
                    want = y
                    break
            got = solve_query(x, dist, m)
            assert got == want, (bs, x, want, got)

    # 隨機用例（m <= 40，b <= 60，lim 可控）
    for _ in range(200):
        n = rnd(5) + 1
        bs = sorted(rnd(60) + 1 for _ in range(n))
        m = bs[0]
        dist = dijkstra_residues(bs, m)
        lim = m * max(bs)
        bf = brute(bs, m, lim)
        assert dist == bf, (bs, dist, bf)
        for _ in range(10):
            x = rnd(lim + 1)
            want = -1
            for y in range(x, lim + 1):
                r = y % m
                if bf[r] != INF and y >= bf[r]:
                    want = y
                    break
            got = solve_query(x, dist, m)
            assert got == want, (bs, x, want, got)

    print("all tests passed")


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        return
    it = iter(data)
    n = int(next(it))
    bs = [int(next(it)) for _ in range(n)]
    q = int(next(it))
    xs = [int(next(it)) for _ in range(q)]
    m = min(bs)
    dist = dijkstra_residues(bs, m)
    out = [str(solve_query(x, dist, m)) for x in xs]
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
