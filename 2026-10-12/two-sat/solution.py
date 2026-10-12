"""2-SAT（2-可滿足性）：蘊含圖 + 強連通分量 + 構造解

題意：
    n 個布爾變數 x_0..x_{n-1}，給 m 個子句，每個子句形如 (a ∨ b)，
    其中 a、b 是**文字**（literal）：正數 +(i+1) 表示 x_i 為真，負數 -(i+1) 表示 x_i 為假。
    問是否存在一組取值讓所有子句同時成立；若存在，輸出一組解。

    本篇同時收錄兩個獨立的 SCC 實作（Kosaraju 迭代版 / Tarjan 迭代版）互相對拍，
    以及一個真實建模範例：**每個任務兩個候選時段的排程問題**。

思路：
    ### 為什麼 2-SAT 是多項式時間，而 3-SAT 是 NP-complete
    子句 (a ∨ b) 等價於兩個蘊含式：
        ¬a → b      （a 不成立，那 b 必須成立）
        ¬b → a
    把每個「文字」當成圖上的一個點，每個蘊含式當成一條有向邊，
    就得到有 2n 個點的**蘊含圖**。

    沿著邊走代表「被迫成立」。若某個變數的 x 與 ¬x 在**同一個強連通分量**裡，
    就代表 x ⟹ ¬x 且 ¬x ⟹ x，自相矛盾 → 無解。
    反過來也成立：只要沒有變數與其否定同分量，就一定有解。

    ### 構造解：關鍵在「分量的拓撲序」
    把蘊含圖縮成 DAG（每個 SCC 縮成一個點）。若在 DAG 上 A → B，
    則「A 成立 ⟹ B 成立」。取拓撲序 rank（邊一定從小 rank 指向大 rank），
    令 **x 為真 ⟺ rank(comp(x)) > rank(comp(¬x))**。

    為什麼這樣對？反證：若某子句 (a ∨ b) 不成立，則 a、b 皆假，即
        rank(comp(a)) ≤ rank(comp(¬a))、rank(comp(b)) ≤ rank(comp(¬b))
    而子句給出邊 ¬a → b、¬b → a，故
        rank(comp(¬a)) ≤ rank(comp(b)) ≤ rank(comp(¬b)) ≤ rank(comp(a)) ≤ rank(comp(¬a))
    全部相等 → 四個文字同一個 SCC → 與「無解已被排除」矛盾。

    ### 拓撲序的正規化（保證 Python 與 C++ 輸出逐字節一致）
    SCC 的**編號方式**取決於實作細節（Kosaraju 與 Tarjan 給出的編號就不同），
    若直接拿編號比較會讓兩個語言算出不同的解。因此：
      1. 輸出時一律用「所在分量的最小點號」這個**與編號無關**的正規形式；
      2. 求拓撲序時用 priority queue，以「分量的最小點號」為鍵。
    因為每個分量的最小點號唯一，整個過程與 SCC 編號完全無關，兩邊必然一致。

    ### 應用：兩時段排程
    每個任務有兩個候選時段（半開區間 [a₁,b₁) 與 [a₂,b₂)），必須選一個；
    不同任務的時段只要重疊就違法。這正是 2-SAT 的經典模型：
    對任意兩個任務 i<j 與任意一對時段組合 (Xᵢ, Yⱼ)，若兩段重疊，
    就加上子句「不能同時選 Xᵢ 與 Yⱼ」，四種組合各一個子句。

    ### 複雜度
    點數 2n、邊數 2m，Kosaraju / Tarjan 都是 O(n + m)，構造解也是 O(n + m)。
    相較之下暴力枚舉是 O(2ⁿ·m)，n=30 就爆了。

應用場景：
    EDA 的電路可滿足性檢查、排班與排課的「二選一」約束、
    編譯器的暫存器配置與指令排程、遊戲關卡的連鎖謎題求解、
    自動規劃（planning）中的互斥資源約束、封裝佈局的通道繞線。

複雜度：
    Kosaraju 迭代版   O(V + E) 時間、O(V + E) 空間
    Tarjan 迭代版     O(V + E) 時間、O(V + E) 空間
    2-SAT 判定與構造  O(n + m)
    兩時段排程建模    O(k²) 個子句，k 為任務數
    暴力枚舉對拍      O(2ⁿ·m)（僅測試用，n ≤ 12）

輸入格式（stdin，全部以空白分隔）：
    n   變數個數
    m   子句個數
    接下來 m 對整數 (a, b)，每對為一個子句 (a ∨ b)；
    文字為非零整數，+(i+1) 表示 x_i 真、-(i+1) 表示 x_i 假。
    非法文字（0 或絕對值 > n）所在的子句直接忽略；token 不足則提前結束。
輸出格式（stdout）：
    第 1 行：SATISFIABLE 或 UNSATISFIABLE
    第 2 行：n 個變數的取值（0/1，空格分隔；UNSAT 時為空行）
    第 3 行：蘊含圖的 SCC 數量
    第 4 行：2n 個點各自所在 SCC 的「最小點號」（空格分隔；正規化表示）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import heapq
import random
import re
import sys
from typing import List, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")


# ---------------------------------------------------------------- SCC 兩種寫法


def kosaraju_scc(n: int, adj: Sequence[Sequence[int]]) -> List[int]:
    """Kosaraju 迭代版：兩次 DFS，第二次在反圖上按完成時間由晚到早走。

    返回 comp[v] = 點 v 所屬的分量編號。
    """
    radj: List[List[int]] = [[] for _ in range(n)]
    for u in range(n):
        for v in adj[u]:
            radj[v].append(u)

    visited = [False] * n
    order: List[int] = []
    for s in range(n):
        if visited[s]:
            continue
        visited[s] = True
        stk = [(s, 0)]
        while stk:
            u, i = stk[-1]
            if i < len(adj[u]):
                stk[-1] = (u, i + 1)          # 顯式棧不持有引用，避免 vector/list 重配置
                v = adj[u][i]
                if not visited[v]:
                    visited[v] = True
                    stk.append((v, 0))
            else:
                order.append(u)
                stk.pop()

    comp = [-1] * n
    c = 0
    for s in reversed(order):
        if comp[s] != -1:
            continue
        comp[s] = c
        stk = [s]
        while stk:
            u = stk.pop()
            for v in radj[u]:
                if comp[v] == -1:
                    comp[v] = c
                    stk.append(v)
        c += 1
    return comp


def tarjan_scc(n: int, adj: Sequence[Sequence[int]]) -> List[int]:
    """Tarjan 迭代版：dfn / low + 顯式棧，返回 comp[v]。"""
    dfn = [0] * n                 # 0 表示未訪問，否則為發現序號（從 1 起）
    low = [0] * n
    on = [False] * n
    comp = [-1] * n
    stk: List[int] = []
    timer = 1
    nc = 0

    for s in range(n):
        if dfn[s]:
            continue
        dfn[s] = low[s] = timer
        timer += 1
        stk.append(s)
        on[s] = True
        work = [(s, 0)]
        while work:
            u, i = work[-1]
            if i < len(adj[u]):
                work[-1] = (u, i + 1)
                v = adj[u][i]
                if not dfn[v]:
                    dfn[v] = low[v] = timer
                    timer += 1
                    stk.append(v)
                    on[v] = True
                    work.append((v, 0))
                elif on[v]:
                    low[u] = min(low[u], dfn[v])
            else:
                work.pop()
                if low[u] == dfn[u]:          # u 是一個分量的根，彈到 u 為止
                    while True:
                        w = stk.pop()
                        on[w] = False
                        comp[w] = nc
                        if w == u:
                            break
                    nc += 1
                if work:                      # 回到父節點，傳遞 low
                    p = work[-1][0]
                    low[p] = min(low[p], low[u])
    return comp


def comp_canonical(comp: Sequence[int]) -> List[int]:
    """把分量編號換成「該分量中最小點號」——與編號方式無關的正規形式。"""
    smallest = {}
    for v, c in enumerate(comp):
        if c not in smallest or v < smallest[c]:
            smallest[c] = v
    return [smallest[c] for c in comp]


# ---------------------------------------------------------------- 2-SAT


def lit_node(lit: int) -> int:
    """文字 → 點號：+(i+1) → 2i（真），-(i+1) → 2i+1（假）。"""
    return 2 * (abs(lit) - 1) + (0 if lit > 0 else 1)


def build_implication(n: int, clauses: Sequence[Tuple[int, int]]) -> List[List[int]]:
    """建蘊含圖：(a ∨ b) ⟹ 邊 ¬a→b 與 ¬b→a。"""
    adj: List[List[int]] = [[] for _ in range(2 * n)]
    for a, b in clauses:
        na, nb = lit_node(a), lit_node(b)
        adj[na ^ 1].append(nb)
        adj[nb ^ 1].append(na)
    return adj


def condensation_rank(n: int, adj: Sequence[Sequence[int]], comp: Sequence[int], ncomp: int) -> List[int]:
    """縮點 DAG 的拓撲序排名。

    以「分量的最小點號」為 priority queue 的鍵 —— 鍵唯一，因此結果與 SCC 編號無關。
    邊 u→v 保證 rank[comp[u]] < rank[comp[v]]。
    """
    cmin = [n] * ncomp
    for v, c in enumerate(comp):
        if v < cmin[c]:
            cmin[c] = v

    cedge: List[set] = [set() for _ in range(ncomp)]
    indeg = [0] * ncomp
    for u in range(n):
        for v in adj[u]:
            cu, cv = comp[u], comp[v]
            if cu != cv and cv not in cedge[cu]:
                cedge[cu].add(cv)
                indeg[cv] += 1

    heap = [(cmin[c], c) for c in range(ncomp) if indeg[c] == 0]
    heapq.heapify(heap)
    rank = [0] * ncomp
    t = 0
    while heap:
        _, c = heapq.heappop(heap)
        rank[c] = t
        t += 1
        for d in sorted(cedge[c]):
            indeg[d] -= 1
            if indeg[d] == 0:
                heapq.heappush(heap, (cmin[d], d))
    assert t == ncomp                          # 縮點後必為 DAG，不會有環
    return rank


def solve_2sat(n: int, clauses: Sequence[Tuple[int, int]]) -> Tuple[bool, List[int]]:
    """解 2-SAT，返回 (是否有解, 一組解)。無解時解為空串列。"""
    adj = build_implication(n, clauses)
    comp = kosaraju_scc(2 * n, adj)
    for v in range(n):
        if comp[2 * v] == comp[2 * v + 1]:     # x 與 ¬x 互通 → 矛盾
            return False, []
    ncomp = max(comp) + 1 if comp else 0
    rank = condensation_rank(2 * n, adj, comp, ncomp)
    # x 為真 ⟺ 真節點所在分量的拓撲序比較晚
    assign = [1 if rank[comp[2 * v]] > rank[comp[2 * v + 1]] else 0 for v in range(n)]
    return True, assign


# ---------------------------------------------------------------- 對拍用暴力


def lit_value(lit: int, assign: Sequence[int]) -> bool:
    """文字在給定取值下是否為真。"""
    v = abs(lit) - 1
    return assign[v] == (1 if lit > 0 else 0)


def check_assign(n: int, clauses: Sequence[Tuple[int, int]], assign: Sequence[int]) -> bool:
    """驗證一組取值是否滿足所有子句。"""
    return all(lit_value(a, assign) or lit_value(b, assign) for a, b in clauses)


def bf_2sat(n: int, clauses: Sequence[Tuple[int, int]]) -> bool:
    """暴力枚舉 2ⁿ 種取值，僅供對拍（n ≤ 12）。"""
    for mask in range(1 << n):
        assign = [(mask >> i) & 1 for i in range(n)]
        if all(lit_value(a, assign) or lit_value(b, assign) for a, b in clauses):
            return True
    return False


def bf_2sat_solution(n: int, clauses: Sequence[Tuple[int, int]]):
    """暴力找一組解（字典序最小者），僅供對拍。"""
    for mask in range(1 << n):
        assign = [(mask >> i) & 1 for i in range(n)]
        if check_assign(n, clauses, assign):
            return assign
    return None


# ---------------------------------------------------------------- 應用：兩時段排程


def overlap(p: Tuple[int, int], q: Tuple[int, int]) -> bool:
    """兩個半開區間 [s,e) 是否重疊（長度 0 的區間不與任何區間重疊）。"""
    return max(p[0], q[0]) < min(p[1], q[1])


def schedule_two_slots(slots: Sequence[Tuple[int, int, int, int]]) -> Tuple[bool, List[int]]:
    """每個任務兩個候選時段（半開區間），任務之間不得重疊，是否可排？

    slots[i] = (a1, b1, a2, b2)：任務 i 可選 [a1,b1) 或 [a2,b2)。
    返回 (是否可排, choice)，choice[i] = 0 選第一個、1 選第二個。
    """
    n = len(slots)
    clauses: List[Tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            ai1, bi1, ai2, bi2 = slots[i]
            aj1, bj1, aj2, bj2 = slots[j]
            if overlap((ai1, bi1), (aj1, bj1)):
                clauses.append((-(i + 1), -(j + 1)))        # 不能同時選 A
            if overlap((ai1, bi1), (aj2, bj2)):
                clauses.append((-(i + 1), (j + 1)))         # 不能 i 選 A 且 j 選 B
            if overlap((ai2, bi2), (aj1, bj1)):
                clauses.append(((i + 1), -(j + 1)))         # 不能 i 選 B 且 j 選 A
            if overlap((ai2, bi2), (aj2, bj2)):
                clauses.append(((i + 1), (j + 1)))          # 不能同時選 B
    sat, assign = solve_2sat(n, clauses)
    if not sat:
        return False, []
    return True, [0 if assign[i] == 1 else 1 for i in range(n)]


def bf_schedule(slots: Sequence[Tuple[int, int, int, int]]) -> Tuple[bool, List[int]]:
    """暴力枚舉所有選擇，僅供對拍（任務數 ≤ 10）。"""
    n = len(slots)
    for mask in range(1 << n):
        chosen = []
        for i in range(n):
            chosen.append((slots[i][0], slots[i][1]) if ((mask >> i) & 1) == 0
                          else (slots[i][2], slots[i][3]))
        ok = True
        for i in range(n):
            for j in range(i + 1, n):
                if overlap(chosen[i], chosen[j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True, [(mask >> i) & 1 for i in range(n)]
    return False, []


# ---------------------------------------------------------------- 輸出


def run(n: int, clauses: Sequence[Tuple[int, int]]) -> str:
    sat, assign = solve_2sat(n, clauses)
    adj = build_implication(n, clauses)
    comp = kosaraju_scc(2 * n, adj)
    ncomp = max(comp) + 1 if comp else 0
    canon = comp_canonical(comp)
    lines = ["SATISFIABLE" if sat else "UNSATISFIABLE"]
    lines.append(" ".join(str(x) for x in assign) if sat else "")
    lines.append(str(ncomp))
    lines.append(" ".join(str(x) for x in canon))
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check_2sat(n: int, clauses: List[Tuple[int, int]]) -> None:
    """Kosaraju 與 Tarjan 分量一致；判定與暴力一致；構造解確實可行。"""
    adj = build_implication(n, clauses)
    c1 = kosaraju_scc(2 * n, adj)
    c2 = tarjan_scc(2 * n, adj)
    assert comp_canonical(c1) == comp_canonical(c2), (n, clauses)
    n1 = (max(c1) + 1) if c1 else 0
    n2 = (max(c2) + 1) if c2 else 0
    assert n1 == n2, (n, clauses)

    sat, assign = solve_2sat(n, clauses)
    assert sat == bf_2sat(n, clauses), (n, clauses)
    if sat:
        assert len(assign) == n
        assert check_assign(n, clauses, assign), (n, clauses, assign)


def run_tests() -> None:
    rng = random.Random(20261012)

    # ---- 固定用例 ----
    # (x0 ∨ x1) ∧ (¬x0 ∨ x1) ∧ (x0 ∨ ¬x1) → 唯一解 x0=1, x1=1
    check_2sat(2, [(1, 2), (-1, 2), (1, -2)])
    assert solve_2sat(2, [(1, 2), (-1, 2), (1, -2)]) == (True, [1, 1])

    # (x0 ∨ x1) ∧ (¬x0 ∨ x1) ∧ (x0 ∨ ¬x1) ∧ (¬x0 ∨ ¬x1) → 無解
    check_2sat(2, [(1, 2), (-1, 2), (1, -2), (-1, -2)])
    assert solve_2sat(2, [(1, 2), (-1, 2), (1, -2), (-1, -2)])[0] is False

    # 空公式：永遠有解
    check_2sat(0, [])
    check_2sat(3, [])
    assert run(0, []) == "SATISFIABLE\n\n0\n"

    # 單位子句 (x ∨ x) 強制 x 為真；矛盾對 (x)∧(¬x) 無解
    check_2sat(1, [(1, 1)])
    assert solve_2sat(1, [(1, 1)]) == (True, [1])
    check_2sat(1, [(1, 1), (-1, -1)])
    assert solve_2sat(1, [(1, 1), (-1, -1)])[0] is False

    # 上面這組的具體輸出（跑過之後寫死）
    assert run(2, [(1, 2), (-1, 2), (1, -2)]) == "SATISFIABLE\n1 1\n2\n0 1 0 1"
    assert run(1, [(1, 1)]) == "SATISFIABLE\n1\n2\n0 1"

    # ---- 隨機對拍：小規模全枚舉 ----
    for _ in range(600):
        n = rng.randint(1, 7)
        m = rng.randint(0, 10)
        clauses = []
        for _ in range(m):
            a = rng.choice([-1, 1]) * rng.randint(1, n)
            b = rng.choice([-1, 1]) * rng.randint(1, n)
            clauses.append((a, b))
        check_2sat(n, clauses)

    # ---- 隨機對拍：限制子句數，讓大 n 也有解/無解兩邊都出現 ----
    for _ in range(300):
        n = rng.randint(8, 12)
        m = rng.randint(0, 6)
        clauses = []
        for _ in range(m):
            a = rng.choice([-1, 1]) * rng.randint(1, n)
            b = rng.choice([-1, 1]) * rng.randint(1, n)
            clauses.append((a, b))
        sat, assign = solve_2sat(n, clauses)
        assert sat == bf_2sat(n, clauses)
        if sat:
            assert check_assign(n, clauses, assign)

    # ---- 隨機圖：單純驗證兩種 SCC 寫法一致 ----
    for _ in range(300):
        v = rng.randint(1, 20)
        e = rng.randint(0, 40)
        adj: List[List[int]] = [[] for _ in range(v)]
        for _ in range(e):
            a, b = rng.randrange(v), rng.randrange(v)
            adj[a].append(b)
        assert comp_canonical(kosaraju_scc(v, adj)) == comp_canonical(tarjan_scc(v, adj))

    # ---- 應用：兩時段排程對拍 ----
    for _ in range(400):
        k = rng.randint(0, 8)
        slots = []
        for _ in range(k):
            a1 = rng.randint(0, 8)
            b1 = a1 + rng.randint(0, 3)
            a2 = rng.randint(0, 8)
            b2 = a2 + rng.randint(0, 3)
            slots.append((a1, b1, a2, b2))
        ok1, ch1 = schedule_two_slots(slots)
        ok2, ch2 = bf_schedule(slots)
        assert ok1 == ok2, slots
        if ok1:
            chosen = [(s[0], s[1]) if ch1[i] == 0 else (s[2], s[3]) for i, s in enumerate(slots)]
            for i in range(k):
                for j in range(i + 1, k):
                    assert not overlap(chosen[i], chosen[j]), (slots, ch1)

    # 經典不可排範例：三個任務，兩兩時段組合都撞
    assert schedule_two_slots([(0, 2, 0, 2), (0, 2, 0, 2), (0, 2, 0, 2)])[0] is False
    # 可排：錯開即可
    assert schedule_two_slots([(0, 1, 5, 6), (1, 2, 6, 7)])[0] is True


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。"""
    tokens = raw.split()

    def read_int(i: int):
        return int(tokens[i]) if i < len(tokens) and NUM_RE.match(tokens[i]) else None

    n = read_int(0)
    m = read_int(1)
    if n is None:
        n = 0
    if m is None:
        m = 0
    if n < 0:
        n = 0
    if m < 0:
        m = 0

    clauses: List[Tuple[int, int]] = []
    for t in range(m):
        a = read_int(2 + 2 * t)
        b = read_int(3 + 2 * t)
        if a is None or b is None:
            break
        if a == 0 or b == 0 or abs(a) > n or abs(b) > n:
            continue                             # 非法文字：整條子句忽略
        clauses.append((a, b))

    sys.stdout.write(run(n, clauses) + "\n")


def main() -> None:
    raw = sys.stdin.read()
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        print("all tests passed")


if __name__ == "__main__":
    main()
