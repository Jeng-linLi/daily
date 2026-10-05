"""開放尋址哈希表（線性 / 二次 / 雙重探測）與 LRU Cache

題意：
    實現一個 **不使用鏈地址法** 的哈希表（string -> int），支持：
      1. `put(key, value)` 插入 / 更新；
      2. `get(key)` 查詢，返回 `(是否命中, 值)`；
      3. `erase(key)` 刪除（必須用**墓碑**標記，不能直接清空槽位）；
      4. 三種探測策略：線性探測 / 二次探測 / 雙重散列；
      5. 負載因子（含墓碑）超過 0.5 時自動擴容再散列；
      6. 統計探測次數，用來直觀比較三種策略的聚集程度。
    另外實現一個 **LRU Cache**（哈希表 + 雙向鏈表），`get` / `put` 均 O(1)。

思路：
    ### 為什麼刪除要用墓碑（tombstone）
    開放尋址裡，查找是沿著探測序列一路走到 **第一個 EMPTY 槽** 才宣告失敗。
    若刪除時直接把槽位設回 EMPTY，會把探測鏈「截斷」，導致該鏈後方的元素永遠查不到。
    所以刪除改成標記 DELETED：
      - 查找遇到 DELETED **繼續往前走**（不能停），但要記住第一個 DELETED 的位置，
        插入時可以複用；
      - 插入遇到 DELETED 時先記下來，直到碰到 EMPTY 才回頭用第一個墓碑。

    ### 三種探測策略
      - 線性探測 `h + i`：緩存友好（連續記憶體），但會出現 **一次聚集**（primary
        clustering），連續的佔用區塊越長越慢。
      - 二次探測 `h + i*i`：消除一次聚集，但仍有 **二次聚集**（同哈希值的 key 走同一條序列）。
        當表長為質數且負載因子 <= 0.5 時，`i^2 mod M` 恰好取到 `(M+1)/2` 個不同值，
        與「至少 M/2 個空槽」必有交集，因此**保證**能找到空位。
      - 雙重散列 `h + i * (1 + h' mod (M-1))`：步長與 M 互質（M 為質數），
        探測序列可以走遍整張表，聚集最少，代價是多了第二次散列、緩存局部性較差。

    ### 負載因子為什麼卡在 0.5
    一方面保證二次探測的正確性（見上），另一方面讓平均探測次數維持在 O(1)。

    ### LRU Cache
    哈希表負責 O(1) 定位，雙向鏈表負責 O(1) 維護「最近使用順序」：
      - 頭部哨兵之後是**最久未使用**的節點，尾端哨兵之前是**最近使用**的；
      - `get` 命中 → 把節點摘下掛到尾端；
      - `put` 已存在 → 更新值並掛到尾端；不存在且已滿 → 淘汰頭部後的第一個節點。

輸入格式（stdin，全部以空白分隔）：
    mode                      linear | quadratic | double
    cap                       初始容量（會被提升到 >= 2 的質數）
    m                         哈希表操作數
    m 行：put key value | get key | del key
    lru_cap                   LRU 容量
    q                         LRU 操作數
    q 行：put key value | get key
輸出格式（stdout）：
    第 1 行：哈希表最終容量
    第 2 行：存活元素個數
    第 3 行：累計探測次數
    第 4 行：所有元素，按 key 字典序，`k:v` 以空白分隔（空表輸出空行）
    接著每個 get 一行：`1 <value>`（命中）或 `0 -1`（未命中）
    再一行：LRU 目前大小
    再一行：LRU 內容，由最久未使用到最近使用，`k:v` 空白分隔（空則空行）
    再接著每個 LRU get 一行：`<value>` 或 `-1`
輸入被截斷時，缺的部分按「0 個操作 / 容量 0」處理。
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import List, Optional, Tuple

# 槽位的三種狀態
EMPTY, OCCUPIED, DELETED = 0, 1, 2
MASK32 = 0xFFFFFFFF

# 只接受 [+-]?digits，與 C++ 版本的 tryLL 完全一致；
# 其它 token 一律視為「輸入到此為止」，避免兩語言一個崩一個不崩。
_INT_RE = re.compile(r"^[+-]?[0-9]+$")


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；不是十進制整數則返回 None。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


# ---------------------------------------------------------------- 雜湊函數
def poly_hash(s: str) -> int:
    """多項式雜湊（基數 131，截斷到 32 位無號）。

    刻意不放 Python 的內建 hash()：它對字串是隨機加鹽的，每次進程結果都不同。
    這裡用確定性雜湊，才能和 C++ 版本逐字節對拍（僅對 ASCII 保證一致）。
    """
    h = 0
    for ch in s:
        h = (h * 131 + ord(ch)) & MASK32
    return h


def is_prime(x: int) -> bool:
    """試除法判質數。"""
    if x < 2:
        return False
    if x % 2 == 0:
        return x == 2
    d = 3
    while d * d <= x:
        if x % d == 0:
            return False
        d += 2
    return True


def next_prime(x: int) -> int:
    """返回 >= x 的最小質數。"""
    if x <= 2:
        return 2
    if x % 2 == 0:
        x += 1
    while not is_prime(x):
        x += 2
    return x


# ---------------------------------------------------------------- 開放尋址哈希表
class OpenAddressingHashTable:
    """開放尋址哈希表（string -> int），支援三種探測策略 + 墓碑 + 自動擴容。"""

    def __init__(self, capacity: int = 8, mode: str = "linear") -> None:
        self.cap: int = max(next_prime(capacity), 2)
        self.mode: str = mode if mode in ("linear", "quadratic", "double") else "linear"
        self.keys: List[str] = [""] * self.cap
        self.vals: List[int] = [0] * self.cap
        self.state: List[int] = [EMPTY] * self.cap
        self.size: int = 0          # 存活元素數（不含墓碑）
        self.used: int = 0          # 被佔用槽位數（含墓碑）
        self.probes: int = 0        # 累計探測次數，用來比較不同策略

    # 第 i 次探測的槽位下標
    def _probe(self, key: str, i: int) -> int:
        base = poly_hash(key) % self.cap
        if self.mode == "linear":
            return (base + i) % self.cap
        if self.mode == "quadratic":
            return (base + i * i) % self.cap
        step = 1 + (poly_hash(key) % (self.cap - 1))
        return (base + i * step) % self.cap

    def _find(self, key: str) -> Tuple[int, bool]:
        """返回 (槽位下標, 是否已存在)。

        查找必須掃到 EMPTY 才能判定「不存在」；中途的 DELETED 要記錄下來供插入複用。
        """
        first_tomb = -1
        for i in range(self.cap):
            idx = self._probe(key, i)
            self.probes += 1
            st = self.state[idx]
            if st == EMPTY:
                return (first_tomb if first_tomb >= 0 else idx), False
            if st == DELETED:
                if first_tomb < 0:
                    first_tomb = idx
            elif self.keys[idx] == key:
                return idx, True
        return (first_tomb if first_tomb >= 0 else -1), False

    def _rehash(self, new_cap: int) -> None:
        """擴容並把存活元素全部重新散列（墓碑會被清掉）。"""
        pairs = [(self.keys[i], self.vals[i])
                 for i in range(self.cap) if self.state[i] == OCCUPIED]
        self.cap = max(next_prime(new_cap), 2)
        self.keys = [""] * self.cap
        self.vals = [0] * self.cap
        self.state = [EMPTY] * self.cap
        self.size = 0
        self.used = 0
        for k, v in pairs:
            self.put(k, v)

    def put(self, key: str, value: int) -> None:
        # 負載因子含墓碑超過 0.5 就先擴容，避免探測鏈無限變長
        if (self.used + 1) * 2 > self.cap:
            self._rehash(self.cap * 2)
        idx, found = self._find(key)
        if found:
            self.vals[idx] = value
            return
        while idx < 0:                      # 理論上不會發生，留作保險
            self._rehash(self.cap * 2)
            idx, found = self._find(key)
        self.keys[idx] = key
        self.vals[idx] = value
        self.state[idx] = OCCUPIED
        self.size += 1
        self.used += 1

    def get(self, key: str) -> Tuple[bool, int]:
        idx, found = self._find(key)
        return (True, self.vals[idx]) if found else (False, -1)

    def erase(self, key: str) -> bool:
        idx, found = self._find(key)
        if not found:
            return False
        self.state[idx] = DELETED           # 墓碑：不能設回 EMPTY，否則截斷探測鏈
        self.keys[idx] = ""
        self.vals[idx] = 0
        self.size -= 1
        return True

    def items(self) -> List[Tuple[str, int]]:
        """所有存活元素，按 key 字典序（保證輸出確定）。"""
        out = [(self.keys[i], self.vals[i])
               for i in range(self.cap) if self.state[i] == OCCUPIED]
        out.sort()
        return out


# ---------------------------------------------------------------- LRU Cache
class _Node:
    """雙向鏈表節點；head / tail 兩個哨兵的 key 為 None。"""
    __slots__ = ("key", "val", "prev", "next")

    def __init__(self, key: Optional[str] = None, val: int = 0) -> None:
        self.key = key
        self.val = val
        self.prev: Optional["_Node"] = None
        self.next: Optional["_Node"] = None


class LRUCache:
    """LRU Cache：哈希表 O(1) 定位 + 雙向鏈表 O(1) 維護順序。"""

    def __init__(self, capacity: int = 2) -> None:
        self.capacity = max(capacity, 0)
        self.table = {}                     # key -> _Node
        self.head = _Node()                 # 哨兵：其後是最久未使用
        self.tail = _Node()                 # 哨兵：其前是最近使用
        self.head.next = self.tail
        self.tail.prev = self.head

    def _detach(self, node: _Node) -> None:
        assert node.prev is not None and node.next is not None
        node.prev.next = node.next
        node.next.prev = node.prev

    def _attach_tail(self, node: _Node) -> None:
        last = self.tail.prev
        assert last is not None
        last.next = node
        node.prev = last
        node.next = self.tail
        self.tail.prev = node

    def get(self, key: str) -> int:
        node = self.table.get(key)
        if node is None:
            return -1
        self._detach(node)
        self._attach_tail(node)             # 命中後變「最近使用」
        return node.val

    def put(self, key: str, value: int) -> None:
        if self.capacity == 0:
            return                          # 容量為 0 的 cache 什麼都不存
        node = self.table.get(key)
        if node is not None:
            node.val = value
            self._detach(node)
            self._attach_tail(node)
            return
        if len(self.table) >= self.capacity:
            victim = self.head.next
            assert victim is not None and victim is not self.tail
            self._detach(victim)
            assert victim.key is not None
            del self.table[victim.key]
        fresh = _Node(key, value)
        self.table[key] = fresh
        self._attach_tail(fresh)

    def size(self) -> int:
        return len(self.table)

    def snapshot(self) -> List[Tuple[str, int]]:
        """由最久未使用到最近使用。"""
        res: List[Tuple[str, int]] = []
        cur = self.head.next
        while cur is not None and cur is not self.tail:
            assert cur.key is not None
            res.append((cur.key, cur.val))
            cur = cur.next
        return res


# ---------------------------------------------------------------- IO 模式
def run_io(data: str) -> None:
    toks = data.split()
    pos = 0

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos < len(toks):
            v = toks[pos]
            pos += 1
            return v
        return None

    def next_int(default: int) -> int:
        v = parse_int(nxt())
        return default if v is None else v

    mode = nxt() or "linear"
    if mode not in ("linear", "quadratic", "double"):
        mode = "linear"
    cap = next_int(8)
    m = next_int(0)

    ht = OpenAddressingHashTable(cap, mode)
    get_lines: List[str] = []
    stopped = False
    for _ in range(m):
        op = nxt()
        if op is None:
            stopped = True
            break
        if op == "put":
            k = nxt()
            val = parse_int(nxt())
            if k is None or val is None:
                stopped = True
                break
            ht.put(k, val)
        elif op == "get":
            k = nxt()
            if k is None:
                stopped = True
                break
            ok, val = ht.get(k)
            get_lines.append("1 %d" % val if ok else "0 -1")
        elif op == "del":
            k = nxt()
            if k is None:
                stopped = True
                break
            ht.erase(k)

    lru_cap = 0 if stopped else next_int(0)
    q = 0 if stopped else next_int(0)

    cache = LRUCache(lru_cap)
    lru_lines: List[str] = []
    for _ in range(q):
        op = nxt()
        if op is None:
            break
        if op == "put":
            k = nxt()
            val = parse_int(nxt())
            if k is None or val is None:
                break
            cache.put(k, val)
        elif op == "get":
            k = nxt()
            if k is None:
                break
            lru_lines.append(str(cache.get(k)))

    print(ht.cap)
    print(ht.size)
    print(ht.probes)
    print(" ".join("%s:%d" % kv for kv in ht.items()))
    for line in get_lines:
        print(line)
    print(cache.size())
    print(" ".join("%s:%d" % kv for kv in cache.snapshot()))
    for line in lru_lines:
        print(line)


# ---------------------------------------------------------------- 測試
def run_tests() -> None:
    # 質數工具
    assert [next_prime(x) for x in (0, 1, 2, 3, 4, 10, 11, 12)] == [2, 2, 2, 3, 5, 11, 11, 13]

    # 墓碑不能截斷探測鏈：先塞滿造成碰撞，刪掉中間的，再查後面的
    for mode in ("linear", "quadratic", "double"):
        ht = OpenAddressingHashTable(5, mode)
        for i, key in enumerate(["aa", "bb", "cc", "dd"]):
            ht.put(key, i + 1)
        assert ht.size == 4
        assert ht.get("cc") == (True, 3)
        assert ht.erase("bb") is True
        assert ht.size == 3
        assert ht.get("bb") == (False, -1)
        assert ht.get("cc") == (True, 3)        # 關鍵：刪掉 bb 後 cc 仍要查得到
        assert ht.get("dd") == (True, 4)
        assert ht.erase("bb") is False          # 重複刪除
        ht.put("bb", 22)                        # 重新插入應複用墓碑
        assert ht.get("bb") == (True, 22)
        assert ht.size == 4

    # 更新既有 key 不會增加 size / used
    ht = OpenAddressingHashTable(8, "linear")
    ht.put("x", 1)
    ht.put("x", 2)
    assert ht.size == 1
    assert ht.get("x") == (True, 2)

    # 三種策略都要能正確處理大量插入 + 刪除 + 再插入
    keys = ["k%03d" % i for i in range(60)]
    for mode in ("linear", "quadratic", "double"):
        ht = OpenAddressingHashTable(4, mode)
        for i, k in enumerate(keys):
            ht.put(k, i)
        assert ht.size == 60
        # 負載因子永遠 <= 0.5（含墓碑）
        assert ht.used * 2 <= ht.cap
        for i, k in enumerate(keys):
            assert ht.get(k) == (True, i)
        for k in keys[:30]:
            assert ht.erase(k) is True
        assert ht.size == 30
        for i, k in enumerate(keys):
            if i < 30:
                assert ht.get(k) == (False, -1)
            else:
                assert ht.get(k) == (True, i)

    # 與 Python dict 隨機對拍
    random.seed(20261005)
    for mode in ("linear", "quadratic", "double"):
        for _ in range(30):
            cap = random.randint(1, 12)
            ht = OpenAddressingHashTable(cap, mode)
            ref: dict = {}
            for _ in range(300):
                key = "s%d" % random.randint(0, 40)
                op = random.randint(0, 2)
                if op == 0:
                    val = random.randint(-50, 50)
                    ht.put(key, val)
                    ref[key] = val
                elif op == 1:
                    ok, val = ht.get(key)
                    if key in ref:
                        assert ok and val == ref[key]
                    else:
                        assert (ok, val) == (False, -1)
                else:
                    gone = ht.erase(key)
                    assert gone == (key in ref)
                    ref.pop(key, None)
                assert ht.size == len(ref)
                assert ht.used * 2 <= ht.cap          # 負載因子不超 0.5
                assert ht.items() == sorted(ref.items())

    # LRU：基本行為
    cache = LRUCache(2)
    assert cache.get("a") == -1
    cache.put("a", 1)
    cache.put("b", 2)
    assert cache.get("a") == 1
    cache.put("c", 3)                       # 淘汰最久未使用的 b
    assert cache.get("b") == -1
    assert cache.get("a") == 1
    assert cache.get("c") == 3
    assert cache.size() == 2
    assert cache.snapshot() == [("a", 1), ("c", 3)]

    # LRU：get 會刷新順序
    cache = LRUCache(3)
    for k, v in (("a", 1), ("b", 2), ("c", 3)):
        cache.put(k, v)
    assert cache.get("a") == 1              # a 變成最近使用
    cache.put("d", 4)                       # 淘汰 b
    assert cache.get("b") == -1
    assert cache.snapshot() == [("c", 3), ("a", 1), ("d", 4)]

    # LRU：put 既有 key 只更新值、不淘汰
    cache = LRUCache(2)
    cache.put("a", 1)
    cache.put("b", 2)
    cache.put("a", 10)
    assert cache.size() == 2
    assert cache.get("a") == 10
    assert cache.snapshot() == [("b", 2), ("a", 10)]

    # LRU：容量 0 什麼都不存
    cache = LRUCache(0)
    cache.put("a", 1)
    assert cache.size() == 0
    assert cache.get("a") == -1
    assert cache.snapshot() == []

    # LRU：容量 1
    cache = LRUCache(1)
    cache.put("a", 1)
    cache.put("b", 2)
    assert cache.size() == 1
    assert cache.get("a") == -1
    assert cache.get("b") == 2

    # LRU 隨機對拍：用 list 當樸素基準（O(n) 但絕對不會錯）
    for _ in range(30):
        cap = random.randint(1, 8)
        cache = LRUCache(cap)
        order: List[str] = []               # 由最久未使用到最近使用
        vals: dict = {}
        for _ in range(200):
            key = "x%d" % random.randint(0, 12)
            if random.randint(0, 1) == 0:
                val = random.randint(0, 99)
                cache.put(key, val)
                if cap == 0:
                    continue
                if key in vals:
                    order.remove(key)
                    vals[key] = val
                else:
                    if len(order) >= cap:
                        vals.pop(order.pop(0), None)
                    vals[key] = val
                order.append(key)
            else:
                got = cache.get(key)
                if key in vals:
                    assert got == vals[key]
                    order.remove(key)
                    order.append(key)
                else:
                    assert got == -1
            assert cache.size() == len(vals)
            assert cache.snapshot() == [(k, vals[k]) for k in order]

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
