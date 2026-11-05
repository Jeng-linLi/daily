"""KM 演算法：二分圖最大權完美匹配（Kuhn–Munkres / 匈牙利算法帶權版）

題意：
    給定 n × n 的權值矩陣 w[i][j]（可為負），要在左部與右部之間
    找出一個「完美匹配」（每個左部點恰好配一個右部點），
    使匹配邊權值之和最大。

    輸入：第一行 n；接下來 n 行、每行 n 個整數（權值矩陣）。
    輸出：一個整數 = 最大權完美匹配的總權值。

思路：
    ### 1. 頂標（potentials）與相等子圖
    給每個左部點 i 維護頂標 u[i]，每個右部點 j 維護頂標 v[j]，
    始終保持不等式 u[i] + v[j] >= w[i][j]（可行性 / 對偶上界）。
    只考慮「緊繃邊」w[i][j] == u[i] + v[j] 構成的相等子圖；
    由對偶理論，相等子圖上的完美匹配權值 = Σu + Σv = 全局最優。
    初始頂標取 u[i] = max_j w[i][j]、v[j] = 0（顯然可行）。

    ### 2. 逐步增廣 + 頂標調整
    為每個左部點 i 嘗試在相等子圖裡找增廣路：
      - 找得到：直接增廣（匈牙利算法的部分）。
      - 找不到：記錄每個未訪問右部點 j 的最小「鬆弛量」
        minv[j] = min{ u[i0] + v[j] - w[i0][j] }，取全域最小的 delta；
        把所有已訪問左部點的 u 減 delta、已訪問右部點的 v 加 delta。
        這一步保證：
          (a) 原有緊繃邊仍然緊繃（兩端同時動，和不變）；
          (b) 可行性不破壞：樹上行 × 非樹列的鬆弛量減 delta 後仍 >= 0
              （delta 是這些鬆弛量的最小值）；非樹行 × 樹列鬆弛量反而增大；
          (c) 至少一條「探索過的失敗交錯邊」鬆弛量歸零（變緊繃），
              相等子圖擴大。
        調整後繼續找增廣路，直到成功為止。

    ### 3. 複雜度
    每個左部點增廣一次，每次增廣最多調整 O(n) 輪頂標，
    每輪掃描右部點 O(n)，總計 O(n^3)。

    ### 4. 對拍
    內建測試用「枚舉所有排列」的 O(n·n!) 暴力（n ≤ 6），
    覆蓋負權、全等值、極端不對稱矩陣與隨機矩陣，逐一比對。

複雜度：
    時間 O(n^3)，空間 O(n^2)（存矩陣）。

應用場景（Applications）：
    - **任務指派**：n 個人 / n 個任務的成本或收益矩陣，找總收益最大的指派
    - **供需配對**：倉庫與門市的一對一調撥、廣告位與投放計劃的匹配
    - 是線上判題系統「Assignment Problem」的標準解法；
      匈牙利算法（無權版）見本倉庫 2026-10-08 Dinic 的二分圖匹配部分
"""

import sys
from itertools import permutations
from typing import List, Optional, Sequence, Tuple

NMAX = 300
W_LIM = 10 ** 9
VAL_LIM = 10 ** 18
INF = 1 << 62


def km_max_weight_matching(w: Sequence[Sequence[int]], n: int) -> int:
    """KM 演算法（e-maxx 潛能版），O(n^3)，回傳最大權完美匹配總權值。

    w 是 1-indexed 使用：w[i][j] 為左部 i 配右部 j 的權值，
    w[0][*] 與 w[*][0] 不使用（佔位 0）。
    """
    # u: 左部頂標；v: 右部頂標；p[j]: 與右部 j 匹配的左部點（0 = 無）
    # 可行性不變式：u[i] + v[j] >= w[i][j]；緊繃邊 <=> 鬆弛量 u+v-w == 0
    u = [0] * (n + 1)
    v = [0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)
    # 初始頂標：u[i] = max_j w[i][j]、v[j] = 0（顯然可行）
    for i in range(1, n + 1):
        u[i] = max(w[i][j] for j in range(1, n + 1)) if n else 0
    for i in range(1, n + 1):
        p[0] = i                    # 虛擬右部 0 掛上當前左部點
        j0 = 0                      # 當前正在擴展的右部點
        minv = [INF] * (n + 1)      # 每個未訪問右部點的最小鬆弛量
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = INF
            j1 = 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = u[i0] + v[j] - w[i0][j]   # 鬆弛量 >= 0
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0        # 記錄走到 j 的來源右部點
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            # 頂標調整：訪問側 u -= delta、v += delta；未訪問側鬆弛量同步縮小
            for j in range(0, n + 1):
                if used[j]:
                    u[p[j]] -= delta
                    v[j] += delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:          # j0 尚未匹配 → 增廣路打通
                break
        # 沿 way 回溯，翻轉增廣路上的匹配邊
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    ans = 0
    for j in range(1, n + 1):
        ans += w[p[j]][j]
    return ans


def brute_max_weight_matching(w: Sequence[Sequence[int]], n: int) -> int:
    """暴力：枚舉右部點的全排列，O(n · n!)。"""
    if n == 0:
        return 0
    best = -VAL_LIM * 2
    cols = list(range(1, n + 1))
    for perm in permutations(cols):
        s = 0
        for i in range(1, n + 1):
            s += w[i][perm[i - 1]]
        if s > best:
            best = s
    return best


# --------------------------------------------------------------- 內建測試

def run_tests() -> None:
    # 固定用例 1：經典指派（LeetCode 題面值，先跑暴力再寫死）
    w1 = [[0, 0, 0, 0, 0],
          [0, 4, 1, 3, 2],
          [0, 2, 0, 5, 3],
          [0, 3, 2, 2, 3],
          [0, 4, 0, 2, 1]]
    got = km_max_weight_matching(w1, 4)
    want = brute_max_weight_matching(w1, 4)
    assert got == want, f"case1 got={got} want={want}"

    # 固定用例 2：含負權（n=3，最優是避開最負的組合）
    w2 = [[0, 0, 0, 0],
          [0, -5, -2, 3],
          [0, 1, -7, 2],
          [0, 0, 4, -1]]
    got = km_max_weight_matching(w2, 3)
    want = brute_max_weight_matching(w2, 3)
    assert got == want, f"case2 got={got} want={want}"

    # 固定用例 3：全零矩陣
    w3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    assert km_max_weight_matching(w3, 2) == 0

    # 固定用例 4：單點
    assert km_max_weight_matching([[0, 0], [0, 7]], 1) == 7

    # 隨機對拍（n ≤ 6）
    rng_state = 12345

    def rnd(mod: int) -> int:
        nonlocal rng_state
        rng_state = (rng_state * 1103515245 + 12345) % (1 << 31)
        return rng_state % mod

    for it in range(200):
        n = rnd(6) + 1
        w: List[List[int]] = [[0] * (n + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                w[i][j] = rnd(21) - 10        # 權值範圍 [-10, 10]
        got = km_max_weight_matching(w, n)
        want = brute_max_weight_matching(w, n)
        assert got == want, f"random#{it} n={n} got={got} want={want} w={w}"


# --------------------------------------------------------------- IO 模式

def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；非純整數或超出 ±10^18 回 None。"""
    if not tok:
        return None
    i = 0
    if tok[0] in "+-":
        i = 1
    if i == len(tok):
        return None
    for ch in tok[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(tok)
    if val > VAL_LIM:
        return None
    return val


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    n = _parse_int(nxt() or "")
    if n is None or not (0 <= n <= NMAX):
        bad()
        return
    w: List[List[int]] = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            x = _parse_int(nxt() or "")
            if x is None or abs(x) > W_LIM:
                bad()
                return
            w[i][j] = x
    emit(str(km_max_weight_matching(w, n)) + "\n")


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
