# -*- coding: utf-8 -*-
"""博弈論與 Sprague-Grundy 函數（mex DP / 組合遊戲和）

題意：多筆遊戲，每筆有 k「堆」，每堆屬於兩種類型之一：
  type 0 取石遊戲：一次可從該堆取 1..3 顆（取完即消失）；
  type 1 拆分遊戲：一次可把大小 x 的堆拆成 i 與 x-i（1 ≤ i < x）兩堆。

每回合選一堆走一步，不能走者輸。求先手是否必勝，並輸出每堆的 SG 值。

理論要點：
  - SG(x) = mex{ SG(後繼局面) }；拆分遊戲的後繼是兩堆之和，
    其 SG 值為 g(i) xor g(x-i)（組合遊戲和的定理）。
  - 取石遊戲（1..3）有閉式解：g(x) = x mod 4。
  - 組合遊戲和（一堆堆獨立走）整體 SG = 各堆 SG 的 XOR；
    XOR ≠ 0 ⟺ 先手必勝。

IO 模式（stdin 有內容時）：
  t
  每筆：k，接著 k 行「type x」（type ∈ {0,1}，0 ≤ x ≤ 3000）
  每筆輸出兩行：
    第一行：k 個 SG 值（空格分隔，與輸入堆順序一致）
    第二行：first（先手必勝）或 second（後手必勝）
無 stdin 時跑內建斷言（與記憶化搜索的極小極大暴力對拍：
「SG ≠ 0 ⟺ 先手必勝」在雙類型混合的隨機小局面上逐一驗證）並輸出 all tests passed。
"""
import sys

LIMIT = 3000


def mex(s):
    m = 0
    while m in s:
        m += 1
    return m


def grundy_tables(limit):
    """通用 mex DP 同時計算兩種遊戲的 SG 表。

    take: g[x] = mex{g[x-1], g[x-2], g[x-3]}
    split: g[x] = mex{g[i] xor g[x-i] | 1 <= i < x}
    """
    gt = [0] * (limit + 1)
    gs = [0] * (limit + 1)
    for x in range(1, limit + 1):
        s1 = set()
        for r in (1, 2, 3):
            if 0 <= x - r:
                s1.add(gt[x - r])
        gt[x] = mex(s1)
        s2 = set()
        for i in range(1, x):
            s2.add(gs[i] ^ gs[x - i])
        gs[x] = mex(s2)
    return gt, gs


GT, GS = grundy_tables(LIMIT)


def sg(t, x):
    return GT[x] if t == 0 else GS[x]


def state_key(state):
    """局面正規化：排序後的 (type, x) 序列，x=0 的堆直接丟棄。"""
    return tuple(sorted((t, x) for t, x in state if x > 0))


def has_moves(state):
    for t, x in state:
        if t == 0 and x > 0:
            return True
        if t == 1 and x >= 2:
            return True
    return False


WIN = {}


def win(state):
    """極小極大暴力：當前玩家有一步能把對手逼進必敗局則必勝。"""
    key = state_key(state)
    if not key:
        return False  # 無子可走：當前玩家輸
    if key in WIN:
        return WIN[key]
    st = list(key)
    res = False
    for idx in range(len(st)):
        t, x = st[idx]
        if t == 0:
            for r in (1, 2, 3):
                if x - r < 0:
                    continue
                nxt = st[:idx] + st[idx + 1:]
                if x - r > 0:
                    nxt.append((0, x - r))
                if not win(nxt):
                    res = True
                    break
        else:
            if x >= 2:
                for i in range(1, x):
                    nxt = st[:idx] + st[idx + 1:] + [(1, i), (1, x - i)]
                    if not win(nxt):
                        res = True
                        break
        if res:
            break
    WIN[key] = res
    return res


def run_tests():
    # 1) 取石閉式解驗證：g(x) = x mod 4
    for x in range(LIMIT + 1):
        assert GT[x] == x % 4, "take-away formula failed at %d" % x

    # 2) 混合局面：XOR != 0 ⟺ 先手必勝（與極小極大暴力對拍）
    seed = 20261108
    mod = 1 << 31

    def rnd(n):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) % mod
        return seed % n

    for _ in range(400):
        k = 1 + rnd(4)
        state = []
        for _ in range(k):
            state.append((rnd(2), rnd(11)))  # x ∈ [0, 10]，總石數小
        xor = 0
        for t, x in state:
            xor ^= sg(t, x)
        assert (xor != 0) == win(state), "mixed position failed: %r" % (state,)

    # 3) 固定手驗用例（實跑取值：split 的 SG 是 0/1 交替——
    #    偶數堆任何拆法都是 奇+偶 或 偶+偶，xor 恆為 0 → mex{0}=1；
    #    奇數堆任何拆法都是 奇+偶，xor 恆為 1 → mex{1}=0）
    assert [sg(0, 5), sg(0, 6), sg(0, 7)] == [1, 2, 3]  # 取石: x mod 4
    assert sg(1, 4) == 1 and sg(1, 6) == 1 and sg(1, 7) == 0


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        print("all tests passed")
        return
    pos = 0

    def read_int():
        nonlocal pos
        v = int(data[pos])
        pos += 1
        return v

    t = read_int()
    out = []
    for _ in range(t):
        k = read_int()
        piles = []
        vals = []
        xor = 0
        for _ in range(k):
            tp = read_int()
            x = read_int()
            if tp not in (0, 1) or x < 0 or x > LIMIT:
                tp, x = 0, 0  # 非法輸入按空堆處理
            g = sg(tp, x)
            piles.append((tp, x))
            vals.append(g)
            xor ^= g
        out.append(" ".join(map(str, vals)))
        out.append("first" if xor != 0 else "second")
    if out:
        sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
