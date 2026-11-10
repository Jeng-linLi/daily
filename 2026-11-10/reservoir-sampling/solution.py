#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""蓄水池抽樣 Reservoir Sampling（Algorithm R）—— 流式均勻採樣。

命令（IO 模式，從 stdin 讀入）：
  sample n k seed
      -> 對 1..n 的流做蓄水池抽樣（隨機決策由自製 LCG 產生），
         輸出被選中的 k 個編號（升序，空格分隔）。

應用場景：A/B 測試分組、日誌流抽樣、大數據集單遍均勻子採樣——
只需 O(k) 內存、單遍掃描、事先不知道流長 n 也能保證每個元素
被選中概率恰為 k/n。

算法（Algorithm R）：
  * 前 k 個元素直接入池；
  * 第 i 個元素（i > k）以 k/i 的概率替換池中隨機一個位置。
  * 歸納可證：任意時刻池中 k 個元素是已見前綴的均勻樣本。

stdin 無內容時跑內建斷言測試，輸出 all tests passed。
"""

import sys


class LCG:
    """自製線性同餘產生器（與倉庫其他題目一致，兩語言行為完全相同）。"""

    def __init__(self, seed):
        self.state = seed & 0x7FFFFFFF

    def next_rand(self, m):
        self.state = (self.state * 1103515245 + 12345) & 0x7FFFFFFF
        return self.state % m


def reservoir_sample(n, k, seed):
    """Algorithm R：返回選中的編號列表（順序為池中順序）。"""
    rng = LCG(seed)
    res = []
    for i in range(1, n + 1):
        if i <= k:
            res.append(i)
        elif rng.next_rand(i) < k:
            res[rng.next_rand(k)] = i
    return res


def solve(n, k, seed):
    """邊界處理 + 抽樣 + 升序輸出。"""
    if n <= 0 or k <= 0:
        return []
    if k >= n:
        return list(range(1, n + 1))
    return sorted(reservoir_sample(n, k, seed))


def run_io(data):
    tokens = data.split()
    pos = 0
    out = []

    def nxt():
        nonlocal pos
        if pos >= len(tokens):
            return None
        t = tokens[pos]
        pos += 1
        return t

    def nxt_int():
        t = nxt()
        if t is None:
            return None
        try:
            return int(t)
        except ValueError:
            return None

    while True:
        cmd = nxt()
        if cmd is None:
            break
        if cmd == 'sample':
            a = nxt_int()
            b = nxt_int()
            c = nxt_int()
            if a is None or b is None or c is None:
                break
            out.append(' '.join(str(x) for x in solve(a, b, c)))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    # 1) 結構性質：恰好 k 個、無重複、升序
    state = 20261110

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    for _ in range(200):
        n = rnd(50) + 1
        k = rnd(n + 1)  # k 可以為 0 或 n
        seed = rnd(1000000)
        r = solve(n, k, seed)
        assert len(r) == min(k, n), (n, k, r)
        assert r == sorted(set(r)), (n, k, r)
        assert all(1 <= x <= n for x in r), (n, k, r)

    # 2) 邊界：k=0 空、k>=n 全選、n=1
    assert solve(10, 0, 1) == []
    assert solve(10, 10, 42) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert solve(10, 99, 42) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert solve(1, 1, 7) == [1]
    assert solve(0, 0, 7) == []

    # 3) 確定性：同一 seed 輸出一致
    assert solve(20, 5, 12345) == solve(20, 5, 12345)

    # 4) 均勻性 Monte Carlo：n=10, k=3，每個元素期望 30% 入選
    #    20000 次實驗 -> 每個編號期望 6000 次，二項 std 約 65，容差 ±10% = 9σ
    trials = 20000
    counts = [0] * 11  # 1..10
    for t in range(trials):
        for x in reservoir_sample(10, 3, t + 1):
            counts[x] += 1
    for x in range(1, 11):
        lo, hi = 5400, 6600
        assert lo <= counts[x] <= hi, (x, counts[x])

    # 5) 組合多樣性：120 種 3-子集應全部出現過（期望每種 166.7 次）
    seen = set()
    for t in range(trials):
        seen.add(tuple(sorted(reservoir_sample(10, 3, t + 1000000))))
    assert len(seen) == 120, len(seen)


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.write(run_io(data))
    else:
        run_tests()
        sys.stdout.write('all tests passed\n')


if __name__ == '__main__':
    main()
