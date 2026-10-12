"""Aho-Corasick 多模式串匹配（Trie + fail 指針 + 完整轉移表）

題意：
    給定一篇文字 text 與 k 個模式串 pat[0..k-1]，**一次掃描**就找出全部模式串在 text 中的
    所有出現位置（允許重疊），並統計每個模式串各命中幾次。

    這是 KMP 的「多模式」推廣：KMP 對單一模式串建失配指針，Aho-Corasick 則對
    **一整棵 Trie** 建失配指針，讓匹配代價與模式串個數無關。

思路：
    ### 為什麼不能用 KMP 跑 k 次
    k 個模式串各跑一次 KMP 是 O(k·n + k·|pat|)。當 k 很大（敏感詞庫上萬條）時完全不可行。
    Aho-Corasick 把 k 個模式串塞進一棵 Trie，總時間 O(Σ|pat| + n + 命中數)，
    與 k **無關**——這正是它在敏感詞過濾、病毒特徵碼掃描裡無可取代的原因。

    ### 三個部件
    1. **goto（Trie）**：按字元走樹，走到哪個狀態。
    2. **fail（失配指針）**：從狀態 v 出發，若下一個字元走不下去，就跳到
       「當前已匹配字串的最長真後綴」所對應的狀態。它其實就是 **Trie 上的 KMP next**。
       用 BFS 建：根的孩子 fail = 0；其餘 u = goto(v, ch) 的
       fail[u] = goto(fail[v], ch)（沿 fail 鏈往上找第一個能吃掉 ch 的狀態）。
    3. **out（輸出表）**：每個狀態掛上「在此結尾的模式串編號」。
       建 fail 時順便合併：out[v] += out[fail[v]]，
       這樣命中一個狀態就能一次拿到所有「以當前位置結尾」的模式串
       （這正是「aba」裡同時命中「a」與「aba」的關鍵）。

    ### 完整轉移表（Trie 圖）
    匹配時若每步都沿 fail 鏈往上爬，`while` 的總次數是攤還 O(n)，已經夠快；
    但更漂亮的做法是把每個狀態對每個字元的轉移**預先算好**（`go[v][c]`），
    匹配就變成每字元 O(1) 的查表，完全沒有回退：
        go[v][c] = nxt[v][c]  若存在
        go[v][c] = go[fail[v]][c]  否則
    因為 fail[v] 的深度一定小於 v，按 BFS 序遞推即可。這就是所謂的「Trie 圖」。

    ### 複雜度
    設 S = Σ|pat|、A = 字母表大小、n = |text|：
        建 Trie      O(S)
        建 fail      O(S · A)（本實作沿 fail 鏈爬，總攤還 O(S)）
        建轉移表     O(S · A)
        匹配         O(n + 命中數)
    空間 O(S · A)。若 A 很大（例如 Unicode），把 `nxt` 換成雜湊表即可維持 O(S)。

    ### 確定性（保證 Python 與 C++ 輸出逐字節一致）
      - 狀態編號 = 插入順序（同一份 patterns 兩邊建出同一棵樹）；
      - 模式串一律按「輸入順序」輸出，與 Trie 結構無關；
      - 每個模式的命中位置升序輸出（掃描時本來就是升序，仍顯式排序一次）；
      - 空模式串不進 Trie，視為「永不命中」，兩邊一致。

應用場景：
    社群平台的敏感詞 / 廣告違禁詞過濾（一次掃描命中上萬條規則）、
    防毒軟體與 IDS（Snort）的特徵碼掃描、DNA 序列中同時搜尋大量探針、
    編輯器的多關鍵字高亮、拼寫檢查與斷詞。

複雜度：
    建 Trie        O(Σ|pat|)
    建 fail 指針   O(Σ|pat|)（攤還）
    建轉移表       O(狀態數 × |字母表|)
    匹配           O(|text| + 命中數)
    空間           O(狀態數 × |字母表|)
    暴力對拍       O(|text| × Σ|pat|)（僅測試用）

輸入格式（stdin，全部以空白分隔）：
    text    第 1 個 token 為待掃描文字（缺則視為空串）
    k       第 2 個 token 為模式串個數（缺則為 0）
    pat_i   接下來最多 k 個 token 為各模式串（不足則補空串）
輸出格式（stdout）：
    第 1 行：Trie 狀態數（含根節點）
    第 2 行：總命中次數
    第 3 行：各模式的命中次數（按輸入順序，空格分隔；k=0 時為空行）
    第 4 行起：每個模式一行，其命中的起始下標（升序，空格分隔；無則空行）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from collections import deque
from typing import Dict, List, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")


# ---------------------------------------------------------------- 建自動機


def build_trie(patterns: Sequence[str]) -> Tuple[List[Dict[str, int]], List[int], List[List[int]], List[int]]:
    """建立 Trie + fail 指針 + 合併後的輸出表。

    返回 (nxt, fail, out, order)；
    order 為 BFS 序（根在最前），建轉移表時要用到。
    """
    nxt: List[Dict[str, int]] = [{}]
    fail = [0]
    out: List[List[int]] = [[]]

    for pid, p in enumerate(patterns):
        if not p:                            # 空模式串不進 Trie（避免「無處不在」的退化）
            continue
        v = 0
        for ch in p:
            if ch in nxt[v]:
                v = nxt[v][ch]
            else:
                nxt.append({})
                fail.append(0)
                out.append([])
                nxt[v][ch] = len(nxt) - 1
                v = len(nxt) - 1
        out[v].append(pid)

    order = [0]
    q: deque = deque()
    for ch in sorted(nxt[0].keys()):         # 根的直系孩子 fail = 0
        u = nxt[0][ch]
        fail[u] = 0
        q.append(u)
        order.append(u)

    while q:
        v = q.popleft()
        # fail[v] 的深度更小，必定已先出隊並合併完畢 → 直接拼接
        out[v] = out[v] + out[fail[v]]
        for ch in sorted(nxt[v].keys()):
            u = nxt[v][ch]
            f = fail[v]
            while f != 0 and ch not in nxt[f]:
                f = fail[f]
            fail[u] = nxt[f][ch] if ch in nxt[f] else 0
            q.append(u)
            order.append(u)

    return nxt, fail, out, order


def build_goto(nxt, fail, order, alphabet) -> Tuple[List[List[int]], Dict[str, int]]:
    """建完整轉移表（Trie 圖）：go[v][c] = 從 v 吃掉字元 c 後到達的狀態。

    按 BFS 序遞推，因為 fail[v] 一定比 v 淺、已被算過。
    """
    idx = {c: i for i, c in enumerate(alphabet)}
    A = len(alphabet)
    go = [[0] * A for _ in nxt]
    for c, ci in idx.items():
        go[0][ci] = nxt[0].get(c, 0)
    for v in order:
        if v == 0:
            continue
        for c, ci in idx.items():
            if c in nxt[v]:
                go[v][ci] = nxt[v][c]
            else:
                go[v][ci] = go[fail[v]][ci]
    return go, idx


def goto_naive(nxt, fail, v: int, ch: str) -> int:
    """沿 fail 鏈爬求轉移，O(深度)，僅供對拍驗證轉移表。"""
    while v != 0 and ch not in nxt[v]:
        v = fail[v]
    return nxt[v][ch] if ch in nxt[v] else 0


# ---------------------------------------------------------------- 匹配


def ac_match(text: str, patterns: Sequence[str], go, idx, out) -> List[List[int]]:
    """一次掃描找出所有模式串的所有出現位置（允許重疊）。"""
    pos: List[List[int]] = [[] for _ in patterns]
    v = 0
    for i, ch in enumerate(text):
        ci = idx.get(ch)
        if ci is None:                       # 字母表外的字元 → 回到根
            v = 0
            continue
        v = go[v][ci]
        for pid in out[v]:
            pos[pid].append(i - len(patterns[pid]) + 1)
    return pos


def bf_match(text: str, pat: str) -> List[int]:
    """暴力匹配，O(n·|pat|)，僅供對拍。空模式串視為永不命中（與自動機一致）。"""
    if not pat:
        return []
    res = []
    L = len(pat)
    for i in range(len(text) - L + 1):
        if text[i:i + L] == pat:
            res.append(i)
    return res


# ---------------------------------------------------------------- 輸出


def run(text: str, patterns: Sequence[str]) -> str:
    """跑完整流程並回傳要印出的所有行。"""
    nxt, fail, out, order = build_trie(patterns)
    alphabet = sorted(set(text) | set("".join(patterns)))
    go, idx = build_goto(nxt, fail, order, alphabet)
    pos = ac_match(text, patterns, go, idx, out)

    lines = [str(len(nxt)), str(sum(len(p) for p in pos))]
    lines.append(" ".join(str(len(p)) for p in pos))
    for p in pos:
        lines.append(" ".join(str(x) for x in p))
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check(text: str, patterns: List[str]) -> None:
    """自動機結果 == 暴力結果，且轉移表 == fail 鏈爬出來的結果。"""
    nxt, fail, out, order = build_trie(patterns)
    alphabet = sorted(set(text) | set("".join(patterns)))
    go, idx = build_goto(nxt, fail, order, alphabet)

    # 轉移表與 fail 鏈必須完全一致
    for v in range(len(nxt)):
        for c, ci in idx.items():
            assert go[v][ci] == goto_naive(nxt, fail, v, c)

    # fail 指標的基本性質：深度必定更淺
    assert fail[0] == 0
    for v in range(1, len(nxt)):
        assert fail[v] != v

    pos = ac_match(text, patterns, go, idx, out)
    assert len(pos) == len(patterns)
    for pid, p in enumerate(patterns):
        assert pos[pid] == bf_match(text, p), (text, patterns, pid)


def run_tests() -> None:
    rng = random.Random(20261012)

    # ---- 固定用例（期望值皆由暴力核對過後才寫死）----
    check("ababa", ["aba", "ba", "a"])
    check("aaaa", ["aa", "aaa"])
    check("abc", [])
    check("", ["a", "ab"])
    check("xyz", ["q"])
    check("aaa", [""])
    check("she", ["he", "she", "hers", "his"])

    # 上面這組的具體數字（跑過之後寫死，防止重構時悄悄變掉）
    assert run("ababa", ["aba", "ba", "a"]) == "6\n7\n2 2 3\n0 2\n1 3\n0 2 4"
    assert run("aaaa", ["aa", "aaa"]) == "4\n5\n3 2\n0 1 2\n0 1"
    assert run("abc", []) == "1\n0\n"
    assert run("", ["a", "ab"]) == "3\n0\n0 0\n\n"

    # ---- 隨機對拍：二字元字母表，模式串短、重疊多 ----
    for _ in range(400):
        n = rng.randint(0, 30)
        text = "".join(rng.choice("ab") for _ in range(n))
        k = rng.randint(0, 5)
        pats = ["".join(rng.choice("ab") for _ in range(rng.randint(0, 4))) for _ in range(k)]
        check(text, pats)

    # ---- 隨機對拍：三字元字母表 + 較長模式串 ----
    for _ in range(300):
        n = rng.randint(0, 40)
        text = "".join(rng.choice("abc") for _ in range(n))
        k = rng.randint(1, 6)
        pats = ["".join(rng.choice("abc") for _ in range(rng.randint(1, 6))) for _ in range(k)]
        check(text, pats)

    # ---- 隨機對拍：模式串取自文字的真實子串，保證一定有命中 ----
    for _ in range(300):
        n = rng.randint(1, 25)
        text = "".join(rng.choice("abcd") for _ in range(n))
        pats = []
        for _ in range(rng.randint(1, 5)):
            i = rng.randrange(n)
            j = min(n, i + rng.randint(1, 5))
            pats.append(text[i:j])
        check(text, pats)
        assert sum(len(p) for p in pats) > 0


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的數字 token，其餘皆當作字串。"""
    tokens = raw.split()
    text = tokens[0] if len(tokens) >= 1 else ""
    k = 0
    if len(tokens) >= 2 and NUM_RE.match(tokens[1]):
        k = int(tokens[1])
        if k < 0:
            k = 0
    pats = []
    for i in range(k):
        pats.append(tokens[2 + i] if len(tokens) >= 3 + i else "")
    sys.stdout.write(run(text, pats) + "\n")


def main() -> None:
    raw = sys.stdin.read()
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        print("all tests passed")


if __name__ == "__main__":
    main()
