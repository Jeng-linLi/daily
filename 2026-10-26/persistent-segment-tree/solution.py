"""可持久化線段樹（主席樹）：靜態區間第 k 小 / 前驅 / 後繼

題意：
    給定一個靜態陣列 a[1..n]，以及 m 個詢問：

        K l r k    a[l..r] 中第 k 小的值（1 <= k <= r-l+1）
        L l r x    a[l..r] 中小於 x 的最大值（不存在輸出 none）
        R l r x    a[l..r] 中大於 x 的最小值（不存在輸出 none）

    限制：1 <= n, m <= 30000，|a_i|, |x| <= 10^9。

思路：
    ### 1. 值域離散化 + 前綴版本（主席樹）
    「區間第 k 小」是經典的靜態問題。先把所有值離散化到值域 1..M，
    對每個前綴 a[1..i] 維護一棵**權值線段樹** T_i：
    T_i 在離散值 p 上的計數 = 前 i 個元素中等於 p 的個數。
    T_i 與 T_{i-1} 只差一次單點 +1，用**函數式（持久化）**寫法
    每次只新建 O(log M) 個節點、其餘與舊版本共用 —— 這就是「主席樹」。

    ### 2. 區間查詢 = 兩個前綴版本相減
    區間 [l, r] 的權值樹就是 T_r 「減去」 T_{l-1}。查第 k 小時
    在兩個版本上同步下降：左子樹的計數差
        cnt = S[left[T_r]] - S[left[T_{l-1}]]
    即區間裡落在左半值域的元素個數。cnt >= k 往左走；
    否則往右走並把 k 換成 k - cnt。每次詢問 O(log M)。

    ### 3. 前驅 / 後繼
    小於 x 的個數 cnt_lt 可用同樣的「版本相減」在值域
    [1, bisect_left(x)) 上求和得到；cnt_lt > 0 時前驅就是第 cnt_lt 小。
    後驅同理：數出 <= x 的個數 cnt_le，若 cnt_le < r-l+1，
    後驅就是第 cnt_le+1 小。三種詢問最後都歸結到同一個 _kth。

    ### 4. 空節點 0 的技巧
    節點池下標 0 當作空節點：S[0] = L[0] = R[0] = 0。
    版本 0（空前綴）直接就是空節點，update 從空節點長出第一層，
    不需要顯式建一棵滿零樹，節點數從 O(M) 降回 O(修改次數·log M)。

複雜度（M = 離散後的值域大小）：
    建樹          時間 O(n log M)    空間 O(n log M) 個節點
    單次詢問      O(log M)

應用場景：
    - OLAP / 資料倉儲：靜態事實表上反覆出現的「區間分位數」
      「區間第 k 大」報表，預處理一次，線上 O(log n) 回答
    - 遊戲排行榜：靜態戰力榜上查某戰力區間的第 k 名、
      緊鄰的上/下一名（前驅 / 後繼）
    - 價格系統：靜態歷史報價序列上查「某時段內第 k 便宜/貴的成交價」
    - 可持久化陣列：把同樣的「路徑複製」技巧用在**下標域**，
      即得支援歷史版本單點修改 / 查詢的可持久化陣列（思路直接遷移）
    - 二維推廣：值域再套一維結構（樹套樹 / 平衡樹）即可回答
      子矩陣第 k 小
"""

import random
import sys
from bisect import bisect_left, bisect_right
from typing import List, Optional, Tuple

NMAX = 30000
MMAX = 30000
VMAX = 10 ** 9


class PersistentSegTree:
    """值域可持久化權值線段樹（主席樹）。

    節點池：L[u] / R[u] 是左右子下標，S[u] 是計數；
    下標 0 是空節點（S[0] = L[0] = R[0] = 0）。
    roots[i] 是前綴 a[1..i] 對應的根（roots[0] = 0）。
    """

    def __init__(self, arr: List[int]) -> None:
        vals = sorted(set(arr))
        self.values = [0] + vals          # 離散值域，1-indexed
        self.m = len(vals)
        self.L: List[int] = [0]
        self.R: List[int] = [0]
        self.S: List[int] = [0]
        self.roots: List[int] = [0]
        for v in arr:
            p = bisect_left(self.values, v, 1)
            self.roots.append(self._update(self.roots[-1], 1, self.m, p))

    def _new_node(self) -> int:
        self.L.append(0)
        self.R.append(0)
        self.S.append(0)
        return len(self.S) - 1

    def _update(self, u: int, lo: int, hi: int, pos: int) -> int:
        """在版本 u 的基礎上把離散值 pos 的計數 +1，返回新根。"""
        v = self._new_node()
        self.L[v] = self.L[u]
        self.R[v] = self.R[u]
        self.S[v] = self.S[u] + 1
        if lo != hi:
            mid = (lo + hi) // 2
            if pos <= mid:
                self.L[v] = self._update(self.L[u], lo, mid, pos)
            else:
                self.R[v] = self._update(self.R[u], mid + 1, hi, pos)
        return v

    def _kth(self, u: int, v: int, lo: int, hi: int, k: int) -> int:
        """兩個版本相減後的第 k 小離散值（保證存在）。"""
        if lo == hi:
            return lo
        mid = (lo + hi) // 2
        cnt = self.S[self.L[v]] - self.S[self.L[u]]
        if cnt >= k:
            return self._kth(self.L[u], self.L[v], lo, mid, k)
        return self._kth(self.R[u], self.R[v], mid + 1, hi, k - cnt)

    def _count_prefix(self, u: int, v: int, lo: int, hi: int, qr: int) -> int:
        """兩個版本相減後，離散值域 [lo, min(hi, qr)] 內的個數。"""
        if qr < lo:
            return 0
        if hi <= qr:
            return self.S[v] - self.S[u]
        mid = (lo + hi) // 2
        return (self._count_prefix(self.L[u], self.L[v], lo, mid, qr)
                + self._count_prefix(self.R[u], self.R[v], mid + 1, hi, qr))

    def kth(self, l: int, r: int, k: int) -> int:
        p = self._kth(self.roots[l - 1], self.roots[r], 1, self.m, k)
        return self.values[p]

    def predecessor(self, l: int, r: int, x: int) -> Optional[int]:
        """[l, r] 中小於 x 的最大值，不存在回 None。"""
        p = bisect_left(self.values, x, 1) - 1     # 離散值域 [1, p] 都 < x
        if p <= 0:
            return None
        cnt = self._count_prefix(self.roots[l - 1], self.roots[r], 1, self.m, p)
        if cnt == 0:
            return None
        return self.values[self._kth(self.roots[l - 1], self.roots[r], 1, self.m, cnt)]

    def successor(self, l: int, r: int, x: int) -> Optional[int]:
        """[l, r] 中大於 x 的最小值，不存在回 None。"""
        p = bisect_right(self.values, x, 1) - 1    # 離散值域 [1, p] 都 <= x
        total = r - l + 1
        cnt = self._count_prefix(self.roots[l - 1], self.roots[r], 1, self.m, p)
        if cnt >= total:
            return None
        return self.values[self._kth(self.roots[l - 1], self.roots[r], 1, self.m, cnt + 1)]


def solve(n: int, arr: List[int], ops: List[Tuple[str, int, int, int]]) -> List[str]:
    """ops 每項是 (類型, l, r, k或x)；回傳每行輸出的字串。"""
    pst = PersistentSegTree(arr)
    out: List[str] = []
    for op in ops:
        t = op[0]
        if t == "K":
            _, l, r, k = op
            if not (1 <= l <= r <= n and 1 <= k <= r - l + 1):
                out.append("error")
                continue
            out.append(str(pst.kth(l, r, k)))
        elif t in ("L", "R"):
            _, l, r, x = op
            if not (1 <= l <= r <= n and abs(x) <= VMAX):
                out.append("error")
                continue
            if t == "L":
                v = pst.predecessor(l, r, x)
            else:
                v = pst.successor(l, r, x)
            out.append("none" if v is None else str(v))
        else:
            out.append("error")
    return out


def brute(arr: List[int], ops: List[Tuple[str, int, int, int]]) -> List[str]:
    """樸素參考實作：每次詢問排序後線性掃。"""
    out: List[str] = []
    for op in ops:
        t, l, r, x = op
        if not (1 <= l <= r <= len(arr)):
            out.append("error")
            continue
        seg = sorted(arr[l - 1:r])
        if t == "K":
            if not (1 <= x <= len(seg)):
                out.append("error")
                continue
            out.append(str(seg[x - 1]))
        elif t == "L":
            if abs(x) > VMAX:
                out.append("error")
                continue
            cand = [v for v in seg if v < x]
            out.append(str(cand[-1]) if cand else "none")
        elif t == "R":
            if abs(x) > VMAX:
                out.append("error")
                continue
            cand = [v for v in seg if v > x]
            out.append(str(cand[0]) if cand else "none")
        else:
            out.append("error")
    return out


def run_tests() -> None:
    # ---- 固定用例（全部手工驗算過） ----
    assert solve(1, [5], [("K", 1, 1, 1), ("L", 1, 1, 4), ("L", 1, 1, 6),
                          ("R", 1, 1, 5), ("R", 1, 1, 4)]) == \
        ["5", "none", "5", "none", "5"]
    # a = [3,1,3,2,1]，排序後 [1,1,2,3,3]
    assert solve(5, [3, 1, 3, 2, 1], [("K", 1, 5, 3), ("L", 1, 5, 3), ("R", 1, 5, 3),
                                       ("K", 2, 4, 2), ("L", 2, 4, 3), ("R", 2, 4, 1)]) == \
        ["2", "2", "none", "2", "2", "2"]
    # 負值：[-5,0,5,-10] 排序 [-10,-5,0,5]
    assert solve(4, [-5, 0, 5, -10], [("K", 1, 4, 1), ("K", 1, 4, 4), ("K", 1, 4, 2),
                                       ("L", 1, 4, -5), ("R", 1, 4, -5), ("L", 1, 4, -10)]) == \
        ["-10", "5", "-5", "-10", "0", "none"]
    # 全部相同
    assert solve(4, [7, 7, 7, 7], [("K", 2, 3, 1), ("L", 2, 3, 8), ("R", 2, 3, 6),
                                    ("L", 2, 3, 7), ("R", 2, 3, 7)]) == \
        ["7", "7", "7", "none", "none"]
    # 非法操作
    assert solve(3, [1, 2, 3], [("K", 2, 1, 1), ("K", 1, 3, 4), ("X", 1, 3, 1)]) == \
        ["error", "error", "error"]

    # ---- 隨機對拍：主席樹 vs 樸素排序 ----
    rng = random.Random(20261026)
    for _ in range(120):
        n = rng.randint(1, 14)
        arr = [rng.randint(-25, 25) for _ in range(n)]
        ops: List[Tuple[str, int, int, int]] = []
        for _ in range(rng.randint(1, 30)):
            t = rng.choice("KLR")
            l = rng.randint(1, n)
            r = rng.randint(l, n)
            if t == "K":
                ops.append(("K", l, r, rng.randint(1, r - l + 1)))
            else:
                ops.append((t, l, r, rng.randint(-30, 30)))
        assert solve(n, arr, ops) == brute(arr, ops), (n, arr, ops)

    # ---- 稍大規模：確認持久化節點沒有互相污染 ----
    rng = random.Random(987654321)
    n = 800
    arr = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(n)]
    ops = []
    for _ in range(800):
        t = rng.choice("KLR")
        l = rng.randint(1, n)
        r = rng.randint(l, n)
        if t == "K":
            ops.append(("K", l, r, rng.randint(1, r - l + 1)))
        else:
            ops.append((t, l, r, rng.randint(-10 ** 6 - 5, 10 ** 6 + 5)))
    assert solve(n, arr, ops) == brute(arr, ops)


# ---------------------------------------------------------------- IO 模式

OP_SET = {"K", "L", "R"}


def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；越界（超出 ±10^18）回 None。"""
    s = tok
    neg = False
    i = 0
    if i < len(s) and s[i] in "+-":
        neg = s[i] == "-"
        i += 1
    if i == len(s):
        return None
    for ch in s[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(s)
    if val > 10 ** 18:
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

    out_lines: List[str] = []
    n = _parse_int(nxt() or "")
    if n is None or not (1 <= n <= NMAX):
        emit("error\n")
        return
    m = _parse_int(nxt() or "")
    if m is None or not (1 <= m <= MMAX):
        emit("error\n")
        return
    arr: List[int] = []
    ok = True
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None or abs(v) > VMAX:
            ok = False
            break
        arr.append(v)
    if not ok:
        emit("error\n")
        return
    ops: List[Tuple[str, int, int, int]] = []
    for _ in range(m):
        # 每個操作一律消費 4 個 token（類型 + 三個整數），
        # 任一環節不合法都在對應行輸出 error，但 token 仍要前進，
        # 這樣後續操作的解析才不會錯位（兩個語言同一條規則）。
        t = nxt()
        a = _parse_int(nxt() or "")
        b = _parse_int(nxt() or "")
        c = _parse_int(nxt() or "")
        if t is None or t not in OP_SET or a is None or b is None or c is None:
            out_lines.append("error")
            continue
        ops.append((t, a, b, c))
    # solve 內部會逐條校驗 l/r/k/x 的合法性
    out_lines.extend(solve(n, arr, ops))
    emit("".join(s + "\n" for s in out_lines))


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
