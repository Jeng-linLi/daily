"""Tarjan 強連通分量：SCC 分解 / 縮點 DAG / 使全圖強連通的最少加邊數

題意：
    給定一張 n 個點、m 條邊的有向圖，要求：
      1. 求出所有**強連通分量**（SCC：分量內任意兩點互相可達）；
      2. 輸出每個點所屬的分量編號（按分量內最小點號遞增重新編號，保證輸出唯一）；
      3. 統計「縮點」後 DAG 的邊數（兩點間重複邊只算一次）；
      4. 最大的強連通分量有多少個點；
      5. 縮點 DAG 的源點數 / 匯點數，以及**最少再加幾條有向邊能讓整張圖強連通**。

思路：
    有向圖的強連通分量把圖「壓縮」成一棵 DAG：分量內互相可達，分量之間單向連通。
    這是很多圖論題的預處理步驟——縮點之後，帶環的問題往往就退化成 DAG 上的問題。

    **Tarjan 算法**做一次 DFS，給每個點記兩個量：
      - `dfn[v]`：第一次訪問到 v 的時間戳；
      - `low[v]`：從 v 出發，沿著「樹邊」往下走、**最多再走一條回邊 / 橫叉邊**，
        能到達的節點中最小的 dfn。
    遞迴過程中把訪問到的點壓進一個棧。當某點 v 滿足 `low[v] == dfn[v]` 時，
    說明 v 的子樹裡沒有任何邊能連到 v 的**祖先**（否則 low[v] 會被拉小），
    於是以 v 為根的子樹中「還留在棧裡」的點剛好構成一個 SCC，把它們全部彈出。
    每個點進棧出棧各一次，複雜度 O(n + m)。

    與 Kosaraju（兩次 DFS，需要反圖）相比，Tarjan 只要一次 DFS、不需要反圖，
    常數更小；本實作把 Tarjan 寫成**迭代版**（顯式棧），
    避免圖退化成一條鏈時 Python 遞迴深度爆掉，C++ 版同樣用迭代寫法，兩邊行為一致。

    第 5 問的結論：縮點 DAG 上有 `src` 個入度為 0 的點、`snk` 個出度為 0 的點，
    則最少需要加 `max(src, snk)` 條邊（k = 1 時為 0）——每次加邊都能把一個源點和一個匯點接上。

輸入格式（stdin，數字按空白分隔即可）：
    n m
    u1 v1
    u2 v2
    …… （共 m 行，1 ≤ u, v ≤ n）
輸出格式（stdout）：
    第 1 行：強連通分量個數 k
    第 2 行：n 個數，第 i 個是點 i 所屬的分量編號（0 … k-1），空格分隔（n = 0 時空行）
    第 3 行：縮點後 DAG 的邊數（去重）
    第 4 行：最大分量包含的點數（k = 0 時為 0）
    第 5 行：源點數 匯點數 最少加邊數（三個數字空格分隔）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from collections import deque
from typing import List, Tuple


def tarjan_scc(n: int, adj: List[List[int]]) -> Tuple[List[int], int]:
    """Tarjan 求強連通分量（迭代版）。

    回傳 (comp, k)：comp[v] 是點 v（1 … n）所屬的原始分量編號，k 為分量個數。
    原始編號的順序取決於 DFS，需要穩定輸出時請再用 relabel 重新編號。
    """
    dfn = [0] * (n + 1)
    low = [0] * (n + 1)
    on_stack = [False] * (n + 1)
    comp = [-1] * (n + 1)
    stk: List[int] = []          # Tarjan 的「待彈出」棧
    timer = 0
    ncomp = 0

    for root in range(1, n + 1):
        if dfn[root]:
            continue
        # 顯式棧的每個元素 = (當前節點, 下一條要處理的邊的下標)
        work: List[Tuple[int, int]] = [(root, 0)]
        while work:
            v, ei = work[-1]
            if ei == 0:                       # 第一次進入 v
                timer += 1
                dfn[v] = low[v] = timer
                stk.append(v)
                on_stack[v] = True
            descended = False
            while ei < len(adj[v]):
                w = adj[v][ei]
                ei += 1
                if dfn[w] == 0:               # 樹邊：先下去處理 w
                    work[-1] = (v, ei)
                    work.append((w, 0))
                    descended = True
                    break
                elif on_stack[w]:             # 回邊 / 橫叉邊：用 dfn 更新 low
                    low[v] = min(low[v], dfn[w])
            if descended:
                continue
            # v 的所有邊都處理完了
            if low[v] == dfn[v]:              # v 是一個 SCC 的根，彈出整個分量
                while True:
                    x = stk.pop()
                    on_stack[x] = False
                    comp[x] = ncomp
                    if x == v:
                        break
                ncomp += 1
            work.pop()
            if work:                          # 把 low[v] 回傳給父節點
                u = work[-1][0]
                low[u] = min(low[u], low[v])
    return comp, ncomp


def relabel(comp: List[int], k: int, n: int) -> List[int]:
    """把分量編號重新排過：按「分量內最小點號」遞增編號，讓輸出與 DFS 順序無關。

    回傳長度 n+1 的陣列，下標 1 … n 對應節點（下標 0 未使用，恆為 -1）。
    """
    first = [n + 1] * k
    for v in range(1, n + 1):
        c = comp[v]
        if v < first[c]:
            first[c] = v
    order = sorted(range(k), key=lambda c: first[c])
    new_id = [-1] * k
    for i, c in enumerate(order):
        new_id[c] = i
    res = [-1] * (n + 1)
    for v in range(1, n + 1):
        res[v] = new_id[comp[v]]
    return res


def kosaraju_scc(n: int, adj: List[List[int]]) -> Tuple[List[int], int]:
    """Kosaraju 求 SCC（第二次 DFS 在反圖上做），測試裡當 Tarjan 的對拍基準。"""
    radj: List[List[int]] = [[] for _ in range(n + 1)]
    for v in range(1, n + 1):
        for w in adj[v]:
            radj[w].append(v)

    order: List[int] = []
    visited = [False] * (n + 1)
    for s in range(1, n + 1):
        if visited[s]:
            continue
        stack = [(s, 0)]
        visited[s] = True
        while stack:
            v, ei = stack[-1]
            if ei < len(adj[v]):
                stack[-1] = (v, ei + 1)
                w = adj[v][ei]
                if not visited[w]:
                    visited[w] = True
                    stack.append((w, 0))
            else:
                order.append(v)
                stack.pop()

    comp = [-1] * (n + 1)
    k = 0
    for s in reversed(order):
        if comp[s] != -1:
            continue
        stack = [s]
        comp[s] = k
        while stack:
            v = stack.pop()
            for w in radj[v]:
                if comp[w] == -1:
                    comp[w] = k
                    stack.append(w)
        k += 1
    return comp, k


def brute_force_scc(n: int, adj: List[List[int]]) -> List[int]:
    """暴力法：對每個點做 BFS，兩點互相可達則同屬一個分量（O(n·(n+m))，只用於測試）。"""
    reach = [[False] * (n + 1) for _ in range(n + 1)]
    for s in range(1, n + 1):
        reach[s][s] = True
        q = deque([s])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if not reach[s][w]:
                    reach[s][w] = True
                    q.append(w)
    comp = [-1] * (n + 1)
    k = 0
    for v in range(1, n + 1):
        if comp[v] != -1:
            continue
        for w in range(v, n + 1):
            if reach[v][w] and reach[w][v]:
                comp[w] = k
        k += 1
    return comp


def condensation_info(n: int, adj: List[List[int]], comp: List[int], k: int) -> Tuple[int, int, int, int]:
    """縮點後 DAG 的統計：回傳 (邊數, 最大分量大小, 源點數, 匯點數)。"""
    size = [0] * k
    for v in range(1, n + 1):
        size[comp[v]] += 1
    edges = set()
    in_deg = [0] * k
    out_deg = [0] * k
    for v in range(1, n + 1):
        for w in adj[v]:
            a, b = comp[v], comp[w]
            if a != b and (a, b) not in edges:
                edges.add((a, b))
                out_deg[a] += 1
                in_deg[b] += 1
    src = sum(1 for i in range(k) if in_deg[i] == 0)
    snk = sum(1 for i in range(k) if out_deg[i] == 0)
    return len(edges), (max(size) if k else 0), src, snk


def solve(n: int, m: int, edges: List[Tuple[int, int]]) -> List[str]:
    """按題目格式算出五行的輸出。"""
    adj: List[List[int]] = [[] for _ in range(n + 1)]
    for u, v in edges:
        adj[u].append(v)

    comp, k = tarjan_scc(n, adj)
    comp = relabel(comp, k, n)                     # 重編號，保證輸出唯一
    dag_edges, biggest, src, snk = condensation_info(n, adj, comp, k)
    add = 0 if k <= 1 else max(src, snk)           # 使全圖強連通的最少加邊數
    return [
        str(k),
        " ".join(str(c) for c in comp[1:]),      # 下標 0 是佔位，不輸出
        str(dag_edges),
        str(biggest),
        f"{src} {snk} {add}",
    ]


# ---------------------------------------------------------------- 輸入輸出


def run_io(data: str) -> None:
    """按題目格式解析 stdin 並輸出結果。"""
    tokens = list(map(int, data.split()))
    if not tokens:
        n = 0
        rest: List[int] = []
    else:
        n, rest = tokens[0], tokens[1:]
    m = rest[0] if rest else 0
    nums = rest[1:]

    edges: List[Tuple[int, int]] = []
    for i in range(min(m, len(nums) // 2)):
        u, v = nums[2 * i], nums[2 * i + 1]
        if 1 <= u <= n and 1 <= v <= n:
            edges.append((u, v))

    sys.stdout.write("\n".join(solve(n, m, edges)) + "\n")


# ---------------------------------------------------------------- 內置測試


def build(n: int, edges: List[Tuple[int, int]]) -> List[List[int]]:
    adj: List[List[int]] = [[] for _ in range(n + 1)]
    for u, v in edges:
        adj[u].append(v)
    return adj


def run_tests() -> None:
    # ---- 固定用例 ----
    # 1→2→3→1 構成一個環，4 單獨，5 被 3 指向
    e1 = [(1, 2), (2, 3), (3, 1), (3, 5)]
    n1 = 5
    adj = build(n1, e1)
    comp, k = tarjan_scc(n1, adj)
    assert k == 3
    comp2, k2 = kosaraju_scc(n1, adj)
    assert k2 == k
    # 與暴力法對拍：同屬一個分量 ⟺ 互相可達
    bf = brute_force_scc(n1, adj)
    for v in range(1, n1 + 1):
        for w in range(1, n1 + 1):
            assert (comp[v] == comp[w]) == (bf[v] == bf[w])
    labeled = relabel(comp, k, n1)
    assert labeled[1:] == [0, 0, 0, 1, 2]
    de, big, src, snk = condensation_info(n1, adj, labeled, k)
    assert (de, big, src, snk) == (1, 3, 2, 2)     # 只有 0→2 一條縮點邊；源點 {0,1}，匯點 {1,2}
    assert max(src, snk) == 2                      # 加 2 條邊可讓全圖強連通

    # 空圖
    assert tarjan_scc(0, [[]]) == ([-1], 0)
    assert solve(0, 0, []) == ["0", "", "0", "0", "0 0 0"]

    # 單點無邊
    assert solve(1, 0, []) == ["1", "0", "0", "1", "1 1 0"]

    # 整張圖強連通
    assert solve(3, 3, [(1, 2), (2, 3), (3, 1)]) == ["1", "0 0 0", "0", "3", "1 1 0"]

    # 一條鏈：每個點自成一個分量，DAG 就是這條鏈
    assert solve(4, 3, [(1, 2), (2, 3), (3, 4)]) == ["4", "0 1 2 3", "3", "1", "1 1 1"]

    random.seed(20260930)
    for _ in range(500):
        n = random.randint(0, 8)
        m = random.randint(0, 12)
        edges = [(random.randint(1, n), random.randint(1, n)) for _ in range(m)] if n else []
        adj = build(n, edges)

        comp, k = tarjan_scc(n, adj)
        comp_k, kk = kosaraju_scc(n, adj)
        assert kk == k                              # Tarjan vs Kosaraju：分量數一致

        bf = brute_force_scc(n, adj)                # Tarjan vs 暴力：分割完全一致
        for v in range(1, n + 1):
            for w in range(1, n + 1):
                assert (comp[v] == comp[w]) == (bf[v] == bf[w])
                assert (comp_k[v] == comp_k[w]) == (bf[v] == bf[w])

        lab = relabel(comp, k, n)
        lab_k = relabel(comp_k, k, n)
        assert lab == lab_k                         # 重編號後兩個算法應給出同一個陣列
        # 分量編號必須按「分量內最小點號」遞增
        first = {}
        for v in range(1, n + 1):
            first.setdefault(lab[v], v)
        assert [first[i] for i in range(k)] == sorted(first[i] for i in range(k))
        assert (sorted(set(lab[1:])) == list(range(k))) if k else lab[1:] == []

        de, big, src, snk = condensation_info(n, adj, lab, k)
        assert big == max([sum(1 for v in range(1, n + 1) if lab[v] == i) for i in range(k)], default=0)
        assert 0 <= de <= k * (k - 1)
        # 縮點 DAG 無環 → 至少有源點也至少有匯點（k ≥ 1）
        if k > 0:
            assert src >= 1 and snk >= 1

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
