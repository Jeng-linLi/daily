"""後綴自動機（Suffix Automaton, SAM）

題意：
    給定一個字串 s（長度 n），建立它的後綴自動機，再回答五種經典查詢：
      1. s 有幾個**不同**的子串；
      2. 另一個字串 t 是否為 s 的子串；
      3. t 在 s 中出現幾次（重疊也算）；
      4. s 與 t 的**最長公共子串**；
      5. s 的所有不同子串中，字典序第 k 小的是哪一個。

思路：
    ### 後綴自動機是什麼
    SAM 是一個**接受 s 所有子串的最小確定有限狀態自動機**（DFA）。
    它把 O(n²) 個子串壓進 O(n) 個狀態：狀態數 ≤ 2n-1、轉移數 ≤ 3n-4。

    每個狀態 v 代表一個 **endpos 等價類**：一組子串，它們在 s 中的「結束位置集合」完全相同。
    狀態上記兩個量：
      - `length[v]`：該等價類中**最長**串的長度；
      - `link[v]` （後綴鏈）：指向「最長串的最長真後綴」所在的狀態。
    於是該狀態代表的所有子串長度恰好是區間 `(length[link[v]], length[v]]`，
    所以 **不同子串總數 = Σ_v (length[v] - length[link[v]])**（v ≠ 根）。

    ### 在線構造（extend）
    加入字元 c 時新建狀態 cur（length = length[last] + 1），
    沿著 last 的後綴鏈往上走，把「還沒有 c 轉移」的狀態補上 c → cur。
    若撞上已有 c 轉移的狀態 p（轉到 q），分兩種情況：
      - `length[p] + 1 == length[q]`：直接 `link[cur] = q`；
      - 否則 q 同時代表長短不一的串，必須**克隆**出一個只代表短串的新狀態 clone，
        並把沿途指向 q 的 c 轉移改指到 clone。
    攤還 O(1)，整體 O(n)。克隆是 SAM 唯一的難點，也是它為什麼能做到線性。

    ### 出現次數：後綴鏈上的「計數上傳」
    構造時每建一個新狀態 cur 就 `cnt[cur] = 1`（代表一個結束位置），克隆的 cnt 為 0。
    之後按 `length` **由大到小**遍歷所有狀態，把 `cnt[v]` 加到 `cnt[link[v]]` 上。
    因為後綴鏈一定指向長度更小的狀態，這樣掃一次就完成拓撲累積，
    `cnt[v]` 即為該等價類每個子串的出現次數（根節點的 cnt 恰為 |s|）。

    ### 第 k 小子串：DAG 上數路徑
    從根出發的每一條路徑都唯一對應一個子串（自動機是確定的，且每條路徑標籤都是子串），
    所以「不同子串數」=「從根出發的不同路徑數」。定義
        `dp[v] = Σ_c (1 + dp[nxt[v][c]])`
    （選字元 c 得到一個長度 1 的串，再加上從下一個狀態繼續走的所有可能），
    按 length 由大到小算完即可。
    要找第 k 小（1-indexed）就從根往下：對每個狀態按**字元排序**逐個試，
    若 `k ≤ 1 + dp[u]` 就往那條轉移走（`k -= 1`），否則 `k -= 1 + dp[u]` 換下一個字元。
    這本質上是前綴樹的前序走訪，順序正是字典序。

    ### 為什麼要限制輸入長度
    不同子串數是 O(n²)，n = 3000 時約 4.5×10⁶，兩個語言都用 64 位元整數才安全，
    因此本題把 s、t 都截斷到 3000 個位元組（兩個語言一致）。

應用場景：
    - 生物資訊：基因序列的最長公共子串、重複片段偵測
    - 抄襲 / 重複內容檢測：兩篇文章的最長公共子串、不同 n-gram 計數
    - 搜尋引擎與 IDE 的子串索引、拼寫校正的候選生成
    - 壓縮演算法（LZ 系列）的重複子串查找、資料去重
    - 編譯器的字串常量池、惡意程式碼特徵碼比對

複雜度：
    構造            O(n · α)，α 為字元集相關常數（map 為 O(log α)）
    不同子串數      O(狀態數)
    子串判定 / 出現次數  O(|t|)
    最長公共子串    O(|t|)
    第 k 小子串     O(答案長度 · log α)
    空間            O(n)

輸入格式（stdin，以空白分隔）：
    s      基底字串（第 1 個 token；缺失視為空字串）
    t      查詢字串（第 2 個 token；缺失視為空字串）
    k      整數，第 k 小子串（第 3 個 token；非數字或缺失視為 0）

    兩個字串都截斷到前 3000 個位元組。位元組以 latin-1 解碼，
    因此任意位元組都合法（含 UTF-8 中文字串，會被視為一串位元組）。
輸出格式（stdout，共 6 行）：
    第 1 行  s 的不同子串數
    第 2 行  t 是否為 s 的子串（1 / 0）；空字串恆為 1
    第 3 行  t 在 s 中的出現次數（重疊計入；空字串為 |s|）
    第 4 行  s 與 t 的最長公共子串長度
    第 5 行  最長公共子串本身（長度 0 時為空行；並列時取 t 中最早出現者）
    第 6 行  字典序第 k 小子串（k 越界時為空行）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import Dict, List, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")
MAX_LEN = 3000                      # 輸入字串截斷長度（兩個語言必須一致）


class SuffixAutomaton:
    """後綴自動機。

    屬性：
        nxt[v]     : {字元 -> 目標狀態}
        link[v]    : 後綴鏈（根為 -1）
        length[v]  : 該 endpos 等價類最長串的長度
        cnt[v]     : 該等價類每個子串在 s 中的出現次數
        dp[v]      : 從 v 出發可讀出的不同非空字串數
    """

    __slots__ = ("nxt", "link", "length", "created", "last", "cnt", "dp")

    def __init__(self, s: str = "") -> None:
        self.nxt: List[Dict[str, int]] = [{}]
        self.link: List[int] = [-1]
        self.length: List[int] = [0]
        self.created: List[int] = []
        self.last = 0
        for ch in s:
            self.extend(ch)
        self._build_cnt()
        self._build_dp()

    # ---------------------------------------------------------------- 構造
    def extend(self, ch: str) -> None:
        """在線加入一個字元（攤還 O(1)）。"""
        cur = len(self.nxt)
        self.nxt.append({})
        self.link.append(0)
        self.length.append(self.length[self.last] + 1)
        p = self.last
        while p != -1 and ch not in self.nxt[p]:
            self.nxt[p][ch] = cur
            p = self.link[p]
        if p == -1:
            self.link[cur] = 0
        else:
            q = self.nxt[p][ch]
            if self.length[p] + 1 == self.length[q]:
                self.link[cur] = q
            else:
                # q 代表的串長短不一，克隆出只代表「短的那些」的狀態
                clone = len(self.nxt)
                self.nxt.append(dict(self.nxt[q]))
                self.length.append(self.length[p] + 1)
                self.link.append(self.link[q])
                while p != -1:
                    it = self.nxt[p].get(ch)
                    if it is None or it != q:
                        break
                    self.nxt[p][ch] = clone
                    p = self.link[p]
                self.link[q] = clone
                self.link[cur] = clone
        self.last = cur
        self.created.append(cur)

    def _order(self) -> List[int]:
        """狀態依 length 由大到小 —— 後綴鏈與轉移都指向 length 更大的狀態，故為拓撲序。"""
        return sorted(range(len(self.nxt)), key=lambda v: -self.length[v])

    def _build_cnt(self) -> None:
        self.cnt = [0] * len(self.nxt)
        for v in self.created:
            self.cnt[v] = 1
        for v in self._order():
            p = self.link[v]
            if p >= 0:
                self.cnt[p] += self.cnt[v]

    def _build_dp(self) -> None:
        self.dp = [0] * len(self.nxt)
        for v in self._order():
            total = 0
            for u in self.nxt[v].values():
                total += 1 + self.dp[u]
            self.dp[v] = total

    # ---------------------------------------------------------------- 查詢
    def walk(self, t: str) -> int:
        """沿自動機走 t；走不通回傳 -1，否則回傳落點狀態。"""
        v = 0
        for ch in t:
            u = self.nxt[v].get(ch)
            if u is None:
                return -1
            v = u
        return v

    def distinct_count(self) -> int:
        """不同非空子串數。兩種算法互相對拍（路徑計數 vs 後綴鏈長度差）。"""
        by_link = sum(self.length[v] - self.length[self.link[v]]
                      for v in range(1, len(self.nxt)))
        assert by_link == self.dp[0]
        return self.dp[0]

    def is_substring(self, t: str) -> bool:
        return self.walk(t) != -1

    def count_occurrences(self, t: str) -> int:
        v = self.walk(t)
        return 0 if v < 0 else self.cnt[v]

    def longest_common_substring(self, t: str) -> Tuple[int, str]:
        """s 與 t 的最長公共子串：長度與字串本身（並列時取 t 中最早出現者）。"""
        v = 0
        matched = 0
        best = 0
        best_pos = -1
        for i, ch in enumerate(t):
            while v != 0 and ch not in self.nxt[v]:
                v = self.link[v]                 # 沿後綴鏈退，縮短目前匹配
                matched = self.length[v]
            u = self.nxt[v].get(ch)
            if u is None:
                v = 0
                matched = 0
            else:
                v = u
                matched += 1
            if matched > best:
                best = matched
                best_pos = i
        return best, ("" if best == 0 else t[best_pos - best + 1:best_pos + 1])

    def kth_substring(self, k: int) -> str:
        """字典序第 k 小子串（1-indexed）；k 越界回傳空字串。"""
        if k < 1 or k > self.dp[0]:
            return ""
        out: List[str] = []
        v = 0
        while k > 0:
            nxt_v = self.nxt[v]
            for ch in sorted(nxt_v):
                u = nxt_v[ch]
                block = 1 + self.dp[u]
                if k <= block:
                    out.append(ch)
                    k -= 1
                    v = u
                    break
                k -= block
            else:                                # pragma: no cover - 理論上不會到這裡
                break
        return "".join(out)


# ---------------------------------------------------------------- 對拍用暴力


def bf_distinct(s: str) -> int:
    return len({s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)})


def bf_occurrences(s: str, t: str) -> int:
    if t == "":
        return len(s)
    ans = 0
    i = s.find(t)
    while i != -1:
        ans += 1
        i = s.find(t, i + 1)
    return ans


def bf_lcs_len(s: str, t: str) -> int:
    subs = {s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)}
    best = 0
    for i in range(len(t)):
        for j in range(i + 1, len(t) + 1):
            if t[i:j] in subs and j - i > best:
                best = j - i
    return best


def bf_kth(s: str, k: int) -> str:
    subs = sorted({s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)})
    return subs[k - 1] if 1 <= k <= len(subs) else ""


# ---------------------------------------------------------------- 輸出


def run(s: str, t: str, k: int) -> str:
    sam = SuffixAutomaton(s)
    length, text = sam.longest_common_substring(t)
    lines = [
        str(sam.distinct_count()),
        "1" if sam.is_substring(t) else "0",
        str(sam.count_occurrences(t)),
        str(length),
        text,
        sam.kth_substring(k),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check(s: str, t: str, k: int) -> None:
    sam = SuffixAutomaton(s)
    assert sam.distinct_count() == bf_distinct(s), s
    assert sam.is_substring(t) == (t in s), (s, t)
    assert sam.count_occurrences(t) == bf_occurrences(s, t), (s, t)
    length, text = sam.longest_common_substring(t)
    assert length == bf_lcs_len(s, t), (s, t)
    if length > 0:
        assert text in s and text in t, (s, t, text)
    assert sam.kth_substring(k) == bf_kth(s, k), (s, k)


def run_tests() -> None:
    rng = random.Random(20261013)

    # ---- 固定用例 ----
    assert run("aba", "babab", 4) == "5\n0\n0\n3\naba\nb"
    assert run("aaaa", "aa", 3) == "4\n1\n3\n2\naa\naaa"
    assert run("", "", 1) == "0\n1\n0\n0\n\n"
    assert run("abc", "", 1) == "6\n1\n3\n0\n\na"
    assert run("a", "a", 1) == "1\n1\n1\n1\na\na"
    assert run("abab", "baba", 5) == "7\n0\n0\n3\nbab\nb"

    check("aba", "babab", 4)
    check("aaaa", "aa", 3)
    check("", "", 1)
    check("abc", "", 1)
    check("abab", "baba", 5)
    check("mississippi", "issipp", 20)

    # ---- 隨機對拍（小字串，暴力列舉）----
    for _ in range(500):
        alpha = rng.choice(["ab", "abc", "abcd"])
        s = "".join(rng.choice(alpha) for _ in range(rng.randint(0, 12)))
        t = "".join(rng.choice(alpha) for _ in range(rng.randint(0, 8)))
        k = rng.randint(-2, 30)
        check(s, t, k)

    # ---- 隨機對拍（大字串，只與暴力比「不同子串數」與「出現次數」）----
    for _ in range(60):
        s = "".join(rng.choice("abcde") for _ in range(rng.randint(60, 200)))
        sam = SuffixAutomaton(s)
        assert sam.distinct_count() == bf_distinct(s)
        # 字典序第 k 小對整條 list 抽樣比對
        subs = sorted({s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)})
        for _ in range(15):
            k = rng.randint(1, len(subs))
            assert sam.kth_substring(k) == subs[k - 1], (s, k)
        for _ in range(15):
            i = rng.randrange(len(s))
            j = rng.randint(i + 1, len(s))
            t = s[i:j]
            assert sam.count_occurrences(t) == bf_occurrences(s, t)
            assert sam.is_substring(t)
        length, text = sam.longest_common_substring(s)
        assert length == len(s) and text == s

    # ---- 邊界：單字元重複、全相異、長度上限 ----
    assert SuffixAutomaton("a" * 300).distinct_count() == 300
    all_diff = "".join(chr(i) for i in range(1, 201))     # 200 個相異位元組
    assert SuffixAutomaton(all_diff).distinct_count() == 200 * 201 // 2
    # 長度上限：3000 只跑 SAM 本身（暴力是 O(n^2) 個子串，3000 會吃掉數 GB），
    # 暴力對拍改用 400：子串總量約 8×10^4，兩個語言都在一秒內結束。
    big = "".join(rng.choice("ab") for _ in range(MAX_LEN))
    sam_big = SuffixAutomaton(big)
    assert 1 < sam_big.distinct_count() <= MAX_LEN * (MAX_LEN + 1) // 2
    mid = "".join(rng.choice("abc") for _ in range(400))
    assert SuffixAutomaton(mid).distinct_count() == bf_distinct(mid)

    # ---- 後綴鏈與狀態數的結構性質 ----
    for _ in range(200):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(1, 40)))
        sam = SuffixAutomaton(s)
        assert len(sam.nxt) <= 2 * len(s), s          # 狀態數 ≤ 2n-1（n ≥ 1）
        for v in range(1, len(sam.nxt)):
            assert 0 <= sam.link[v] < len(sam.nxt)
            assert sam.length[sam.link[v]] < sam.length[v]
        assert sum(sam.cnt[v] for v in range(len(sam.nxt)) if sam.link[v] == -1) == len(s)


def run_io(raw: str) -> None:
    """stdin 解析：前兩個 token 是字串，第三個必須是嚴格整數否則視為 0。"""
    tokens = raw.split()
    s = tokens[0][:MAX_LEN] if len(tokens) >= 1 else ""
    t = tokens[1][:MAX_LEN] if len(tokens) >= 2 else ""
    k = 0
    if len(tokens) >= 3 and NUM_RE.match(tokens[2]):
        k = int(tokens[2])

    out = run(s, t, k) + "\n"
    emit(out)


def emit(text: str) -> None:
    """以二進位輸出（latin-1），避免 Windows 把 \n 換成 \r\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
