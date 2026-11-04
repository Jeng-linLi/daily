"""可持久化 01-Trie（Persistent Binary Trie）：區間最大異或

題意：
    給定長度為 n 的數列 a[1..n]（0 <= a_i < 2^20），以及 q 筆查詢：
        l r x  —— 求 max{ a_i XOR x : l <= i <= r }
    輸入：第一行 n q，第二行 n 個整數，接下來 q 行每行 l r x（1 起算）。
    輸出：每筆查詢一行答案。

思路：
    ### 1. 為什麼要「可持久化」
    靜態全局最大異或用普通 01-Trie 就能做：把 x 按位從高到低，
    每一步貪心走「與 x 該位相反」的孩子。但本題限制了**區間** [l, r]，
    普通單棵 Trie 無法知道「哪些數落在區間裡」。

    ### 2. 前綴版本 + 可持久化
    令版本 i 的 Trie = 依次插入 a_1, ..., a_i 之後的 Trie（版本 0 是空樹）。
    插入時只「複製被修改路徑上的節點」，其餘孩子直接共享舊版本的指針，
    因此每個版本只多 O(B) 個節點（B = 位數 = 20），總節點數 O(nB)。
    每個節點維護 cnt：**以它為根的子樹裡被插入過的數的個數**。

    ### 3. 版本差 = 區間
    「值落在 a[l..r] 中的數」= 版本 r 減去版本 l-1：
    對 Trie 上同一個前綴節點位置，cnt_r - cnt_{l-1} 就是區間內
    帶該前綴的數的個數。查詢時同時在兩個版本上向下走：
        想走的方向 d（與 x 該位相反）若 cntR[d] - cntOld[d] > 0，
        說明區間內存在帶該前綴的數 → 走 d 並把答案該位置 1；
        否則只能走 x 該位的方向。
    每筆查詢 O(B)。

    ### 4. 邊界
    l = 1 時舊版本是版本 0（空樹，所有 cnt = 0），邏輯自然成立；
    區間內至少一個數（輸入校驗保證 l <= r），貪心走一定有落點。

複雜度：
    建樹   O(n * B)  節點（B = 20）
    查詢   每筆 O(B)
    空間   O(n * B)

應用場景（Applications）：
    - **區間異或查詢**：遊戲數值對比、訊號 / 校驗碼比對中「差異最大」的樣本
    - **離線日誌分析**：在時間視窗 [l, r] 內找與給定指紋差異最大的記錄
    - **安全**：在時間段內找與已知簽名漢明距離最遠（異或最大）的特徵
    - **資料庫**：位圖索引上的範圍 + 異或組合查詢
"""

import sys
from typing import List, Optional, Sequence, Tuple

BITS = 20
VAL_MAX = 1 << BITS          # 1048576
NMAX = 50000
QMAX = 50000


class PersistentTrie:
    """可持久化 01-Trie：節點 0 是空節點（ch = [0,0], cnt = 0）。"""

    def __init__(self) -> None:
        self.ch: List[List[int]] = [[0, 0]]
        self.cnt: List[int] = [0]
        self.roots: List[int] = [0]   # roots[i] = 插入前 i 個數後的根

    def _new_node(self) -> int:
        self.ch.append([0, 0])
        self.cnt.append(0)
        return len(self.ch) - 1

    def insert(self, prev_root: int, val: int) -> int:
        """在版本 prev_root 基礎上插入 val，回傳新版本根節點。"""
        cur = self._new_node()
        root = cur
        p = prev_root
        self.cnt[cur] = self.cnt[p] + 1
        for b in range(BITS - 1, -1, -1):
            d = (val >> b) & 1
            o = d ^ 1
            # 另一側孩子直接共享舊版本
            self.ch[cur][o] = self.ch[p][o]
            q = self.ch[p][d]
            self.ch[cur][d] = self._new_node()
            cur = self.ch[cur][d]
            p = q
            self.cnt[cur] = self.cnt[p] + 1
        return root

    def build(self, values: Sequence[int]) -> None:
        for v in values:
            self.roots.append(self.insert(self.roots[-1], v))

    def query_max_xor(self, l: int, r: int, x: int) -> int:
        """max{ a_i XOR x : l <= i <= r }（l / r 為 1 起算，l <= r）。"""
        u = self.roots[r]          # 新版本
        v = self.roots[l - 1]      # 舊版本（要扣掉的部分）
        ch, cnt = self.ch, self.cnt
        ans = 0
        for b in range(BITS - 1, -1, -1):
            d = (x >> b) & 1
            want = d ^ 1
            cr = ch[u][want]
            co = ch[v][want]
            if cnt[cr] - cnt[co] > 0:
                ans |= 1 << b
                u, v = cr, co
            else:
                u, v = ch[u][d], ch[v][d]
        return ans


def brute_force(a: Sequence[int], l: int, r: int, x: int) -> int:
    """參考實作：O(r-l+1) 直接枚舉。"""
    best = -1
    for i in range(l - 1, r):
        cur = a[i] ^ x
        if cur > best:
            best = cur
    return best


def run_tests() -> None:
    # ---- 固定用例 ----
    t = PersistentTrie()
    t.build([1, 2, 3])
    # a = [1,2,3]，與 0 的最大異或：1^2 = 3
    assert t.query_max_xor(1, 3, 0) == 3
    # 與 1 的最大異或：2^1 = 3
    assert t.query_max_xor(1, 3, 1) == 3
    # 區間只有一個數：答案就是 a_i ^ x
    assert t.query_max_xor(2, 2, 5) == 2 ^ 5
    # 區間 [3,3]：3 ^ 7 = 4
    assert t.query_max_xor(3, 3, 7) == 4

    t2 = PersistentTrie()
    t2.build([0, 0, 0])
    assert t2.query_max_xor(1, 3, 0) == 0
    assert t2.query_max_xor(1, 3, 7) == 7

    t3 = PersistentTrie()
    t3.build([10, 10, 10, 3])
    # [1,3] 全是 10 → 10 ^ x
    assert t3.query_max_xor(1, 3, 0) == 10
    # [1,4] 有 3 → max(10^0, 3^0) = 10；x=9 時 max(10^9=3, 3^9=10) = 10
    assert t3.query_max_xor(1, 4, 9) == 10

    # ---- 與暴力對拍（自製 LCG，兩語言同一序列）----
    state = 20261103

    def nxt(mod: int) -> int:
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mod

    for _ in range(200):
        n = 1 + nxt(12)                       # 1..12
        vals = [nxt(64) for _ in range(n)]    # < 64
        trie = PersistentTrie()
        trie.build(vals)
        for _ in range(10):
            l = 1 + nxt(n)
            r = l + nxt(n - l + 1)
            x = nxt(64)
            got = trie.query_max_xor(l, r, x)
            want = brute_force(vals, l, r, x)
            assert got == want, (vals, l, r, x, got, want)

    # 大值域測試（< 2^20）
    for _ in range(60):
        n = 1 + nxt(60)
        vals = [nxt(VAL_MAX) for _ in range(n)]
        trie = PersistentTrie()
        trie.build(vals)
        for _ in range(8):
            l = 1 + nxt(n)
            r = l + nxt(n - l + 1)
            x = nxt(VAL_MAX)
            got = trie.query_max_xor(l, r, x)
            want = brute_force(vals, l, r, x)
            assert got == want, (vals, l, r, x, got, want)


# ---------------------------------------------------------------- IO 模式

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
    if val > 10 ** 18:
        return None
    return val


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0
    out: List[str] = []

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    n_opt = _parse_int(nxt() or "")
    if n_opt is None or not (1 <= n_opt <= NMAX):
        bad()
        return
    q_opt = _parse_int(nxt() or "")
    if q_opt is None or not (1 <= q_opt <= QMAX):
        bad()
        return
    n, q = n_opt, q_opt

    a: List[int] = []
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None or not (0 <= v < VAL_MAX):
            bad()
            return
        a.append(v)

    trie = PersistentTrie()
    trie.build(a)

    for _ in range(q):
        l = _parse_int(nxt() or "")
        r = _parse_int(nxt() or "")
        x = _parse_int(nxt() or "")
        if l is None or r is None or x is None:
            bad()
            return
        if not (1 <= l <= r <= n) or not (0 <= x < VAL_MAX):
            bad()
            return
        out.append(str(trie.query_max_xor(l, r, x)))

    for line in out:
        emit(line + "\n")


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
