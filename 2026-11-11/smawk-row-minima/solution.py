#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMAWK 演算法：完全單調矩陣（totally monotone matrix）逐列最小值搜尋。

命令（IO 模式，從 stdin 讀入）：
  minima n m x0 ... x_{n-1} y0 ... y_{m-1}
      -> 定義矩陣 A[i][j] = (x[i] - y[j])^2，其中 x、y 都是遞增序列。
         此時 A 滿足反 Monge 條件（逆 Monge）：
           A[i][j] + A[i+1][j+1] >= A[i][j+1] + A[i+1][j]
         進而是完全單調矩陣：若 A[i][j] >= A[i][j']（j < j'），
         則對所有 i' > i 也有 A[i'][j] >= A[i'][j']。
      -> 對每一列 i 輸出一行「最小值 argmin」：
         即對每個 x[i] 找出使距離平方最小的 y[j]（最近鄰配對）。

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * 逐列掃全部欄是 O(n·m)；SMAWK 利用完全單調性做到 O(n + m) 次取值。
  * Reduce 步驟：只用奇數編號的列，維護一個候選欄 stack——
    當前奇數列在 stack 頂欄不優於新欄（>= 語義，平手丟舊欄）時彈棧；
    stack 滿了就丟棄新欄。被丟棄的欄對「所有剩下的列」都被某存活欄支配，
    因此每列的最小值一定留在候選欄裡。
  * 遞迴解出奇數列的 argmin；偶數列的 argmin 被相鄰兩個奇數列的
    argmin 上下夾住（argmin 隨列號單調不減），在該視窗內二分：
    比較相鄰兩欄，較小方向保留（<= 時往左收，> 時往右收）。
  * 應用：Monge/凹凸代價函數的 1D/1D 動態規劃加速、
    斷點搜尋（breakpoint searching）、最近鄰配對等。
"""

import sys


def smawk(rows, cols, A):
    """在完全單調矩陣 A 上，回傳 dict：row -> argmin（限於 cols 中）。

    A(i, j) 是取值函數（延遲求值，SMAWK 的複雜度以取值次數計）。
    """
    if len(rows) == 1:
        # 單列：直接線性掃描（O(m)，總複雜度仍是 O(n + m)）
        r = rows[0]
        best = cols[0]
        for c in cols:
            if A(r, c) < A(r, best):
                best = c
        return {r: best}

    # --- Reduce：用奇數列篩選候選欄 ---
    odd = rows[1::2]
    stack = []
    for c in cols:
        while stack and A(odd[len(stack) - 1], stack[-1]) >= A(odd[len(stack) - 1], c):
            stack.pop()
        if len(stack) < len(odd):
            stack.append(c)
    odd_ans = smawk(odd, stack, A)

    # --- Interpolate：偶數列在相鄰奇數列答案之間二分 ---
    col_pos = {c: t for t, c in enumerate(cols)}
    ans = dict(odd_ans)
    for idx, r in enumerate(rows):
        if idx % 2 == 1:
            continue
        lo = col_pos[odd_ans[rows[idx - 1]]] if idx > 0 else 0
        hi = col_pos[odd_ans[rows[idx + 1]]] if idx + 1 < len(rows) else len(cols) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if A(r, cols[mid]) <= A(r, cols[mid + 1]):
                hi = mid
            else:
                lo = mid + 1
        ans[r] = cols[lo]
    return ans


def nearest_squares(x, y):
    """對每個 x[i] 找最接近的 y[j]，回傳 [(最小值平方, argmin), ...]（SMAWK 版）。

    注意：先把 y 去重成嚴格遞增序列再跑 SMAWK——y 有重複時，(x-y)^2 依
    「索引」看在下降段會出現平台（如 16,16,16,9,9,4,4,1），不是嚴格單峰，
    相鄰比較二分會停在平台邊緣而漏掉真正的最小值。去重後每列嚴格單峰
    （至多一組對稱平手），二分才成立；結果再映射回原 y 的首次出現索引。
    """
    n = len(x)
    if n == 0 or len(y) == 0:
        return []
    yd = []
    first_idx = {}
    for j, v in enumerate(y):
        if v not in first_idx:
            first_idx[v] = j
            yd.append(v)

    def A(i, j):
        d = x[i] - yd[j]
        return d * d

    argmins = smawk(list(range(n)), list(range(len(yd))), A)
    out = []
    for i in range(n):
        j = argmins[i]
        d = x[i] - yd[j]
        out.append((d * d, first_idx[yd[j]]))
    return out


def nearest_squares_brute(x, y):
    """O(n·m) 暴力版（取最左 argmin），僅供內建測試對拍。"""
    out = []
    for xi in x:
        best_j = 0
        best_v = (xi - y[0]) ** 2
        for j in range(1, len(y)):
            v = (xi - y[j]) ** 2
            if v < best_v:
                best_v = v
                best_j = j
        out.append((best_v, best_j))
    return out


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
        if cmd == 'minima':
            n = nxt_int()
            m = nxt_int()
            if n is None or m is None:
                break
            x = []
            for _ in range(n):
                v = nxt_int()
                if v is None:
                    break
                x.append(v)
            y = []
            for _ in range(m):
                v = nxt_int()
                if v is None:
                    break
                y.append(v)
            if len(x) < n or len(y) < m:
                break
            for v, j in nearest_squares(x, y):
                out.append('%d %d' % (v, j))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    # 1) 固定小例（期望值以實跑結果為準：x=4 最近 y=2，x=9 最近 y=7）
    res = nearest_squares([1, 4, 9], [2, 7])
    assert res == [(1, 0), (4, 0), (4, 1)], res
    assert nearest_squares([], [1, 2]) == []
    assert nearest_squares([1], []) == []

    # 2) 與 O(nm) 暴力對拍：隨機遞增序列，小值域製造大量平手
    state = 20261111

    def rnd(mm):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mm

    for _ in range(300):
        n = rnd(30) + 1
        m = rnd(30) + 1
        base = rnd(10)
        x = sorted(base + rnd(15) for _ in range(n))
        y = sorted(base + rnd(15) for _ in range(m))
        got = nearest_squares(x, y)
        want = nearest_squares_brute(x, y)
        assert [g[0] for g in got] == [w[0] for w in want], (x, y, got, want)
        # argmin 位置上也必須真的是最小值（位置可能與暴力不同，值必須相同）
        for i, (v, j) in enumerate(got):
            assert v == (x[i] - y[j]) ** 2
            assert v == want[i][0]

    # 3) 較大規模：值必須正確、位置必須自洽（平手時位置可能與暴力不同）
    for _ in range(20):
        n = rnd(300) + 1
        m = rnd(300) + 1
        x = sorted(rnd(10 ** 6) for _ in range(n))
        y = sorted(rnd(10 ** 6) for _ in range(m))
        got = nearest_squares(x, y)
        want = nearest_squares_brute(x, y)
        assert [g[0] for g in got] == [w[0] for w in want]
        for i, (v, j) in enumerate(got):
            assert (x[i] - y[j]) ** 2 == v


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
