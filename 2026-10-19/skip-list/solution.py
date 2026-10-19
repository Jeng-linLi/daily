"""跳表（Skip List）：以隨機層高達成期望 O(log n) 的有序多重集合

題意：
    維護一個**有序多重集合**（允許重複元素），支援以下操作，並與「有序陣列 + 二分」
    的參考實作逐項對拍，確保跳表的行為完全正確。

    操作碼（每行 `op x`）：
        1  insert x    插入 x（允許重複）
        2  erase  x    刪除一個 x（不存在則忽略）
        3  count  x    x 的出現次數
        4  rank   x    嚴格小於 x 的元素個數
        5  kth    k    第 k 小的元素（k 從 1 起算），越界回傳 -1
        6  pred   x    嚴格小於 x 的最大元素，不存在回傳 -1
        7  succ   x    嚴格大於 x 的最小元素，不存在回傳 -1

思路：
    ### 1. 為什麼需要跳表
    有序陣列的插入/刪除是 O(n)；平衡樹（紅黑樹、AVL）能做到 O(log n) 但旋轉規則複雜。
    跳表用「**機率平衡**」換來簡潔：每個節點隨機一個層高 lv（P(lv ≥ k) = 2^{-(k-1)}），
    第 i 層的指標串起一條「稀疏高速公路」，越高層跳過越多節點。
    搜尋時從最高層往右走，走不過去就下降一層，期望步數 O(log n)。

    ### 2. span（跨度）欄位：讓「排名」也能 O(log n)
    純指標只能回答「在不在」，要回答 rank / kth 必須知道每條邊**跨過幾個節點**。
    每個節點在第 i 層額外存 `span[i]`：從自己沿第 i 層指標走一步，等於在第 0 層往前跳幾個節點。
    插入時若節點只出現在 0..lv-1 層，那麼更高層的每條邊都被「多跨過一個節點」：
        - i <  lv：新節點切斷舊邊，兩段跨度是 (rank[0]-rank[i]) 與 舊跨度 - (rank[0]-rank[i])
        - i ≥  lv：舊邊仍然存在，跨度直接 +1
    刪除則是逆向操作：`update[i].span[i] += x.span[i] - 1`。

    ### 3. 隨機層高必須跨語言一致
    Python 的 `random` 與 C++ 的 `mt19937` 給出的序列不同，兩邊的樹形會長得不一樣
    （雖然答案應該相同，但對拍時就無法逐位元組比對中間過程）。
    這裡兩個版本都用**同一組 LCG 常數**自製亂數，取 64 位元的第 33 位當硬幣：
        state = state * 6364136223846793005 + 1442695040888963407  (mod 2^64)
    擲到 1 就停止，否則層高 +1（上限 32 層）。

    ### 4. 不變量檢查
    對每一層的每一條邊，都用第 0 層指標實際走一趟，確認 `span[i]` 等於真實距離。
    這個檢查是 O(n²)，但能把「跨度算錯但答案碰巧對」的 bug 抓出來。

應用場景：
    - Redis 的 `zset`（有序集合）底層就是跳表（zskiplist），`ZRANK` / `ZRANGE` 靠 span 實作
    - LevelDB / RocksDB 的 MemTable（併發友善的記憶體索引）
    - 需要區間鎖細粒度、或想要無鎖（lock-free）實作的場合，跳表比平衡樹容易做
    - 任何「需要有序 + 頻繁插入刪除 + 排名查詢」的中介層（排行榜、時間軸、價格檔位）

複雜度：
    搜尋 / 插入 / 刪除     期望 O(log n) 時間
    rank / kth / pred / succ  期望 O(log n) 時間（靠 span 欄位）
    count                  期望 O(log n) 時間（count_leq - count_less）
    空間                   期望 O(n)（每個節點期望 2 條邊），最壞 O(n log n)
"""

import random
import re
import sys
from bisect import bisect_left, bisect_right, insort
from typing import List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

MASK64 = (1 << 64) - 1
LCG_MUL = 6364136223846793005
LCG_ADD = 1442695040888963407
MAX_LEVEL = 32
INIT_SEED = 0x9E3779B97F4A7C15 & MASK64

QMAX = 300          # IO 模式下的操作數上限
VMAX = 10 ** 6      # 元素值上限（夾到 [0, VMAX]，讓 -1 可以當哨兵）


# ---------------------------------------------------------------- 跳表節點


class Node:
    """一個跳表節點：key + 每層的後繼指標 + 每層的跨度。"""

    __slots__ = ("key", "forward", "span")

    def __init__(self, key: Optional[int], level: int) -> None:
        self.key = key
        self.forward: List[Optional["Node"]] = [None] * level
        self.span: List[int] = [0] * level


class SkipList:
    def __init__(self, seed: int = INIT_SEED) -> None:
        self.state = seed & MASK64
        self.head = Node(None, MAX_LEVEL)   # 哨兵頭，固定 32 層
        self.level = 1                      # 目前「有用」的最高層
        self.length = 0

    # -- 自製 LCG（與 C++ 版逐位元一致） --

    def _rand64(self) -> int:
        self.state = (self.state * LCG_MUL + LCG_ADD) & MASK64
        return self.state

    def _random_level(self) -> int:
        lv = 1
        while lv < MAX_LEVEL:
            if (self._rand64() >> 33) & 1:
                break
            lv += 1
        return lv

    # -- 內部：尋找每一層的前驅 --

    def _find_preds(self, key: int, strict: bool) -> Tuple[List[Optional[Node]], List[int]]:
        """回傳 (update, rank)：update[i] 是第 i 層的最後一個「不越界」節點，
        rank[i] 是 update[i] 前面有幾個元素。strict=True 用 `<`，False 用 `<=`。"""
        update: List[Optional[Node]] = [None] * MAX_LEVEL
        rank: List[int] = [0] * MAX_LEVEL
        x: Optional[Node] = self.head
        for i in range(self.level - 1, -1, -1):
            if i + 1 < self.level:
                rank[i] = rank[i + 1]
            while True:
                nxt = x.forward[i]
                if nxt is None:
                    break
                if (nxt.key < key) if strict else (nxt.key <= key):
                    rank[i] += x.span[i]
                    x = nxt
                else:
                    break
            update[i] = x
        return update, rank

    # -- 基本操作 --

    def insert(self, key: int) -> None:
        update, rank = self._find_preds(key, True)
        lv = self._random_level()
        if lv > self.level:
            for i in range(self.level, lv):
                update[i] = self.head
                update[i].span[i] = self.length
                rank[i] = 0
            self.level = lv
        node = Node(key, lv)
        for i in range(lv):
            node.forward[i] = update[i].forward[i]
            update[i].forward[i] = node
            node.span[i] = update[i].span[i] - (rank[0] - rank[i])
            update[i].span[i] = (rank[0] - rank[i]) + 1
        for i in range(lv, self.level):
            update[i].span[i] += 1
        self.length += 1

    def erase(self, key: int) -> bool:
        update, _rank = self._find_preds(key, True)
        x = update[0].forward[0]
        if x is None or x.key != key:
            return False
        for i in range(self.level):
            if update[i].forward[i] is x:
                update[i].span[i] += x.span[i] - 1
                update[i].forward[i] = x.forward[i]
            else:
                update[i].span[i] -= 1
        while self.level > 1 and self.head.forward[self.level - 1] is None:
            self.level -= 1
        self.length -= 1
        return True

    def count_less(self, key: int) -> int:
        cnt = 0
        x: Optional[Node] = self.head
        for i in range(self.level - 1, -1, -1):
            while True:
                nxt = x.forward[i]
                if nxt is None or nxt.key >= key:
                    break
                cnt += x.span[i]
                x = nxt
        return cnt

    def count_leq(self, key: int) -> int:
        cnt = 0
        x: Optional[Node] = self.head
        for i in range(self.level - 1, -1, -1):
            while True:
                nxt = x.forward[i]
                if nxt is None or nxt.key > key:
                    break
                cnt += x.span[i]
                x = nxt
        return cnt

    def count(self, key: int) -> int:
        return self.count_leq(key) - self.count_less(key)

    def kth(self, k: int) -> int:
        if k < 1 or k > self.length:
            return -1
        traversed = 0
        x: Optional[Node] = self.head
        for i in range(self.level - 1, -1, -1):
            while True:
                nxt = x.forward[i]
                if nxt is None or traversed + x.span[i] > k:
                    break
                traversed += x.span[i]
                x = nxt
        if x is self.head:
            return -1
        return x.key

    def predecessor(self, key: int) -> int:
        r = self.count_less(key)
        return self.kth(r) if r >= 1 else -1

    def successor(self, key: int) -> int:
        return self.kth(self.count_leq(key) + 1)

    def to_list(self) -> List[int]:
        out: List[int] = []
        x = self.head.forward[0]
        while x is not None:
            out.append(x.key)
            x = x.forward[0]
        return out

    def levels(self) -> List[int]:
        """第 0 層每個節點的層高（用來對拍隨機序列是否一致）。"""
        out: List[int] = []
        x = self.head.forward[0]
        while x is not None:
            out.append(len(x.forward))
            x = x.forward[0]
        return out

    def validate(self) -> bool:
        """對每條邊用第 0 層實走一趟，確認 span 欄位正確。"""
        for i in range(self.level):
            x: Optional[Node] = self.head
            while x.forward[i] is not None:
                target = x.forward[i]
                d = 0
                y = x
                while y is not target:
                    y = y.forward[0]
                    if y is None:
                        return False
                    d += 1
                    if d > self.length + 1:
                        return False
                if d != x.span[i]:
                    return False
                x = target
        return True


# ---------------------------------------------------------------- 參考實作（有序陣列）


class RefSet:
    def __init__(self) -> None:
        self.a: List[int] = []

    def insert(self, x: int) -> None:
        insort(self.a, x)

    def erase(self, x: int) -> bool:
        i = bisect_left(self.a, x)
        if i < len(self.a) and self.a[i] == x:
            del self.a[i]
            return True
        return False

    def count(self, x: int) -> int:
        return bisect_right(self.a, x) - bisect_left(self.a, x)

    def rank(self, x: int) -> int:
        return bisect_left(self.a, x)

    def kth(self, k: int) -> int:
        return self.a[k - 1] if 1 <= k <= len(self.a) else -1

    def predecessor(self, x: int) -> int:
        i = bisect_left(self.a, x)
        return self.a[i - 1] if i > 0 else -1

    def successor(self, x: int) -> int:
        i = bisect_right(self.a, x)
        return self.a[i] if i < len(self.a) else -1


# ---------------------------------------------------------------- 主流程


def run(ops: List[Tuple[int, int]]) -> str:
    sl = SkipList()
    ref = RefSet()
    res: List[int] = []
    ok = True

    for op, x in ops:
        if op == 1:
            sl.insert(x)
            ref.insert(x)
        elif op == 2:
            a = sl.erase(x)
            b = ref.erase(x)
            if a != b:
                ok = False
        elif op == 3:
            sv, rv = sl.count(x), ref.count(x)
            res.append(sv)
            if sv != rv:
                ok = False
        elif op == 4:
            sv, rv = sl.count_less(x), ref.rank(x)
            res.append(sv)
            if sv != rv:
                ok = False
        elif op == 5:
            sv, rv = sl.kth(x), ref.kth(x)
            res.append(sv)
            if sv != rv:
                ok = False
        elif op == 6:
            sv, rv = sl.predecessor(x), ref.predecessor(x)
            res.append(sv)
            if sv != rv:
                ok = False
        elif op == 7:
            sv, rv = sl.successor(x), ref.successor(x)
            res.append(sv)
            if sv != rv:
                ok = False

    if sl.length != len(ref.a):
        ok = False
    if sl.to_list() != ref.a:
        ok = False
    if not sl.validate():
        ok = False

    lines = []
    lines.append(" ".join(str(v) for v in res))
    tail = [str(sl.length)]
    tail.extend(str(v) for v in sl.to_list())
    lines.append(" ".join(tail))
    lines.append("1" if ok else "0")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例 1：基本插入 ----
    out = run([(1, 5), (1, 3), (1, 7), (3, 3), (4, 5), (5, 2)])
    assert out == "1 1 5\n3 3 5 7\n1\n", out

    # ---- 固定用例 2：重複元素 ----
    out = run([(1, 2), (1, 2), (1, 2), (3, 2), (4, 2), (4, 3)])
    assert out == "3 0 3\n3 2 2 2\n1\n", out

    # ---- 固定用例 3：前驅後繼 ----
    out = run([(1, 10), (1, 20), (1, 30), (6, 20), (7, 20), (6, 10), (7, 30)])
    assert out == "10 30 -1 -1\n3 10 20 30\n1\n", out

    # ---- 固定用例 4：刪除不存在 / 刪空 ----
    out = run([(1, 1), (2, 9), (2, 1), (3, 1), (5, 1)])
    assert out == "0 -1\n0\n1\n", out

    # ---- 空操作 ----
    assert run([]) == "\n0\n1\n", repr(run([]))

    # ---- 層高序列跨語言一致性（回歸用） ----
    sl = SkipList()
    want = [sl._random_level() for _ in range(8)]
    assert want == [3, 1, 2, 4, 1, 3, 1, 1], want

    # ---- 隨機對拍 ----
    rng = random.Random(20261019)
    for _ in range(700):
        q = rng.randint(0, 40)
        ops = []
        for _ in range(q):
            op = rng.randint(1, 7)
            x = rng.randint(0, 12)
            ops.append((op, x))
        out = run(ops)
        assert out.split("\n")[-2] == "1", (ops, out)

    # ---- 大量插入後隨機查詢 ----
    for trial in range(60):
        rng2 = random.Random(9000 + trial)
        ops = [(1, rng2.randint(0, 200)) for _ in range(rng2.randint(0, 120))]
        for _ in range(30):
            ops.append((rng2.randint(2, 7), rng2.randint(0, 210)))
        out = run(ops)
        assert out.split("\n")[-2] == "1", (trial, out)


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出，避免 Windows 把 \\n 換成 \\r\\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    toks = raw.split()
    idx = 0

    def next_int(default: int) -> int:
        nonlocal idx
        if idx >= len(toks):
            return default
        t = toks[idx]
        idx += 1
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-(10 ** 18), min(10 ** 18, v))

    q = next_int(0)
    if q < 0:
        q = 0
    if q > QMAX:
        q = QMAX
    ops: List[Tuple[int, int]] = []
    for _ in range(q):
        op = next_int(0)
        x = next_int(0)
        if op < 1 or op > 7:
            op = 0
        if x < 0:
            x = 0
        if x > VMAX:
            x = VMAX
        ops.append((op, x))
    emit(run(ops))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
