"""樹形 DP 與換根 DP（最大獨立集 / 直徑 / 重心 / 各點最遠距離 / 各點距離和）

題意：
    給定一棵樹（或森林，無向無環圖，n 個點 m 條邊），求：
      1. **最大獨立集**：最大的點集，使集合中任意兩點之間沒有邊；
      2. **樹的直徑**：最長路徑的長度（邊數）與具體路徑；
      3. **樹的重心**：刪掉該點後，剩下的最大連通塊大小 ≤ n/2 的點（可能 1 或 2 個）；
      4. **換根 DP（rerooting）之 maxdist**：對每個點 u，求 u 到「所在連通分量內最遠點」的距離；
      5. **換根 DP 之 sumdist**：對每個點 u，求 u 到「所在連通分量內所有點」的距離之和。

思路：
    ### 樹形 DP 的基本套路
    任取一點為根，把「無根樹」轉成「有根樹」，於是可以自底向上（後序）做 DP：
    每個點 u 的狀態只依賴它的子節點。`dp[u] = 合併(所有子節點 v 的 dp[v])`。
    複雜度通常是 O(n)，因為每條邊只被處理常數次。

    ### 最大獨立集
    每個點有「選 / 不選」兩種狀態：
        dp1[u] = 1 + Σ dp0[v]                （選 u ⟹ 所有子節點都不能選）
        dp0[u] = Σ max(dp0[v], dp1[v])       （不選 u ⟹ 子節點隨意）
    答案 max(dp0[root], dp1[root])。這是「樹上背包 / 狀態機 DP」最簡單的一例。
    方案重建：父節點被選則子節點強制不選；否則取 dp1 ≥ dp0 者（平手優先選）。

    ### 樹的直徑（兩次 BFS / DFS）
    從任意點 s 出發走到最遠的點 a，再從 a 出發走到最遠的點 b，
    則 a–b 必定是一條直徑。對森林則對每個連通分量各做一次並取最大。
    （這個性質對**樹**成立；一般圖不成立，一般圖要用 Floyd 或 Johnson。）

    ### 樹的重心
    以任意點為根，記 sz[u] 為 u 的子樹大小。刪掉 u 後會分成若干塊：
    每個子節點 v 對應一塊大小 sz[v]，父親方向對應一塊大小 compSize - sz[u]。
    取這些塊的最大值，若 ≤ compSize/2 則 u 是重心。重心必存在且最多 2 個。

    ### 換根 DP（本題重點）
    上面那些量若「以每個點為根各算一次」是 O(n^2)。換根 DP 把它降到 O(n)：
    先做一次自底向上（bottom-up）拿到「向下」的資訊，再做一次自頂向下（top-down）
    把「向上」的資訊補給每個子節點，於是每個點都同時擁有「往下看」和「往上看」的視角。

    - **maxdist**：對每個 u 維護向下最長的兩條路徑 `best1 / best2`（分別記是從哪個子節點來的）。
      對子節點 v，其「向上的最遠距離」為 `1 + max(up[u], best_excluding_v[u])`，
      其中 best_excluding_v 就是「不從 v 這個分支走」的最佳值（who1[u] == v 時取 best2，否則取 best1）。
      答案 `maxdist[u] = max(best1[u], up[u])`。
    - **sumdist**：記 `sz[u]` 與 `sub[u]`（u 到子樹內所有點的距離和）。
      經典公式：把根從 u 換到子節點 v 時，v 這一側的 sz[v] 個點距離各減 1，
      其餘 compSize - sz[v] 個點距離各加 1，所以
          sumdist[v] = sumdist[u] - sz[v] + (compSize - sz[v])
      O(1) 轉移，一趟 top-down 全部算完。

    森林的情況：對每個連通分量獨立計算，compSize 用該分量的大小。

應用場景：
    社交網路中「影響力最大」的節點（距離和最小 = 最接近所有人的點，即重心類指標）、
    伺服器 / 倉庫選址（樹上最小化最大延遲 = 直徑中心）、
    依賴樹上的任務調度（最大獨立集 = 互不衝突的最大任務集）、
    文件目錄樹與組織架構的聚合統計。

複雜度：
    最大獨立集 / 重心 / 換根 DP    O(n) 時間、O(n) 空間
    直徑（兩次 BFS，對每個分量）   O(n) 時間（森林亦然）、O(n) 空間
    暴力基準（從每個點 BFS）       O(n · (n + m))

輸入格式（stdin，全部以空白分隔）：
    n m
    m 行：u v        （0-indexed 無向邊；重邊與自環會被忽略；越界邊忽略）
    n = 0 或輸入不足時視為空圖。遇到非整數 token 視為輸入結束。
輸出格式（stdout）：
    第 1 行：最大獨立集大小
    第 2 行：獨立集節點（升序，空格分隔）
    第 3 行：直徑長度（邊數；空圖為 0）
    第 4 行：直徑路徑（空格分隔；空圖輸出空行）
    第 5 行：重心個數
    第 6 行：重心（升序，空格分隔）
    第 7 行：maxdist[0..n-1]（空格分隔）
    第 8 行：sumdist[0..n-1]（空格分隔）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
from collections import deque
import sys
from typing import List, Optional, Tuple

_INT_RE = re.compile(r"^[+-]?[0-9]+$")   # 嚴格整數規則，與 C++ 的 tryLL 完全一致


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；非十進制整數則返回 None（視為輸入到此為止）。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


def build_adj(n: int, edges: List[Tuple[int, int]]) -> List[List[int]]:
    """建立鄰接表：去重、去自環、去越界邊，並將鄰點升序排序（保證確定性）。"""
    st = set()
    for u, v in edges:
        if u == v:
            continue
        if u < 0 or u >= n or v < 0 or v >= n:
            continue
        if u > v:
            u, v = v, u
        st.add((u, v))
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v in st:
        adj[u].append(v)
        adj[v].append(u)
    for i in range(n):
        adj[i].sort()
    return adj


def forest_data(n: int, adj: List[List[int]]):
    """以每個分量的最小編號點為根做迭代 DFS，返回 (parent, root, order)。

    parent[root] = -2 表示該點是根；order 為先序（父必在子之前）。
    """
    parent = [-1] * n
    root = list(range(n))
    order: List[int] = []
    for s in range(n):
        if parent[s] != -1:
            continue
        parent[s] = -2
        stack = [s]
        while stack:
            u = stack.pop()
            order.append(u)
            for v in adj[u]:
                if parent[v] == -1:
                    parent[v] = u
                    root[v] = s
                    stack.append(v)
    return parent, root, order


# ---------------------------------------------------------------- 最大獨立集


def max_independent_set(n: int, adj: List[List[int]],
                        parent: List[int], order: List[int]) -> Tuple[int, List[int]]:
    """樹上最大獨立集：返回 (大小, 節點升序列表)。"""
    dp0 = [0] * n
    dp1 = [0] * n
    for u in reversed(order):
        dp1[u] = 1
        for v in adj[u]:
            if parent[v] == u:
                dp1[u] += dp0[v]
                dp0[u] += dp0[v] if dp0[v] >= dp1[v] else dp1[v]
    chosen = [False] * n
    for u in order:
        p = parent[u]
        blocked = p >= 0 and chosen[p]
        chosen[u] = (not blocked) and dp1[u] >= dp0[u]
    nodes = [u for u in range(n) if chosen[u]]
    return len(nodes), nodes


def mis_bruteforce(n: int, adj: List[List[int]]) -> int:
    """最大獨立集暴力版：枚舉所有子集（n ≤ 14）。"""
    best = 0
    for mask in range(1 << n):
        ok = True
        for u in range(n):
            if not (mask >> u & 1):
                continue
            for v in adj[u]:
                if mask >> v & 1:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            best = max(best, bin(mask).count("1"))
    return best


# ---------------------------------------------------------------- 直徑 / 重心


def bfs_far(n: int, adj: List[List[int]], s: int):
    """從 s 出發 BFS，返回 (最遠點編號最小的那個, dist, parent)。"""
    dist = [-1] * n
    par = [-1] * n
    dist[s] = 0
    q = deque([s])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                par[v] = u
                q.append(v)
    best = s
    for i in range(n):
        if dist[i] > dist[best]:
            best = i
    return best, dist, par


def tree_diameter(n: int, adj: List[List[int]], root: List[int]) -> Tuple[int, List[int]]:
    """森林中最長路徑：返回 (長度, 路徑)。平手時取根編號最小的分量。"""
    best_len = -1                                   # 用 -1 當哨兵，確保第一個分量必定被採用
    best_path: List[int] = []
    for s in range(n):
        if root[s] != s:
            continue
        a, _, _ = bfs_far(n, adj, s)
        b, dist, par = bfs_far(n, adj, a)
        length = dist[b]
        path = []
        cur = b
        while cur != -1:                            # 沿 parent 從 b 一路走回 a
            path.append(cur)
            if cur == a:
                break
            cur = par[cur]
        if length > best_len:
            best_len = length
            best_path = path
    return (0 if best_len < 0 else best_len), best_path


def centroids(n: int, adj: List[List[int]], parent: List[int], root: List[int],
              sz: List[int], comp: List[int]) -> List[int]:
    """所有重心：刪除後最大連通塊 ≤ compSize/2。"""
    res = []
    for u in range(n):
        mx = comp[root[u]] - sz[u]
        for v in adj[u]:
            if parent[v] == u:
                if sz[v] > mx:
                    mx = sz[v]
        if mx * 2 <= comp[root[u]]:
            res.append(u)
    return res


def centroids_bruteforce(n: int, adj: List[List[int]], comp: List[int],
                         root: List[int]) -> List[int]:
    """重心暴力版：真的把 u 刪掉，再數 u 所在分量內各連通塊的大小。

    注意森林要**逐分量**判定：只統計與 u 同分量的塊，且與 comp[root[u]] 比較。
    若把其它分量也算進來，塊大小會被無關的分量撐大，定義就錯了。
    """
    res = []
    for u in range(n):
        r = root[u]
        seen = [False] * n
        seen[u] = True
        mx = 0
        for s in range(n):
            if seen[s] or root[s] != r or s == u:
                continue
            cnt = 0
            q = deque([s])
            seen[s] = True
            while q:
                x = q.popleft()
                cnt += 1
                for y in adj[x]:
                    if not seen[y] and y != u and root[y] == r:
                        seen[y] = True
                        q.append(y)
            if cnt > mx:
                mx = cnt
        if mx * 2 <= comp[r]:
            res.append(u)
    return res


# ---------------------------------------------------------------- 換根 DP


def reroot(n: int, adj: List[List[int]], parent: List[int], root: List[int],
           order: List[int]):
    """換根 DP：返回 (maxdist, sumdist, sz, comp)。

    maxdist[u] = u 到所在分量內最遠點的距離；
    sumdist[u] = u 到所在分量內所有點的距離和。
    """
    comp = [0] * n
    for u in range(n):
        comp[root[u]] += 1

    best1 = [0] * n          # 向下的最長路徑
    best2 = [0] * n          # 向下的次長路徑（來自不同子節點）
    who1 = [-1] * n
    sz = [1] * n
    sub = [0] * n            # u 到子樹內所有點的距離和

    for u in reversed(order):
        for v in adj[u]:
            if parent[v] == u:
                sz[u] += sz[v]
                sub[u] += sub[v] + sz[v]
                val = best1[v] + 1
                if val > best1[u]:
                    best2[u] = best1[u]
                    best1[u] = val
                    who1[u] = v
                elif val > best2[u]:
                    best2[u] = val

    up = [0] * n
    maxdist = [0] * n
    sumdist = [0] * n
    for u in order:
        maxdist[u] = best1[u] if best1[u] > up[u] else up[u]
        if parent[u] == -2:
            sumdist[u] = sub[u]
        for v in adj[u]:
            if parent[v] == u:
                excl = best2[u] if who1[u] == v else best1[u]
                via_up = up[u] if up[u] > excl else excl
                up[v] = 1 + via_up
                sumdist[v] = sumdist[u] - sz[v] + (comp[root[u]] - sz[v])
    return maxdist, sumdist, sz, comp


def brute_distances(n: int, adj: List[List[int]]) -> Tuple[List[int], List[int]]:
    """從每個點各跑一次 BFS，得到 maxdist / sumdist 的暴力基準。"""
    md = [0] * n
    sd = [0] * n
    for s in range(n):
        _, dist, _ = bfs_far(n, adj, s)
        mx = 0
        total = 0
        for i in range(n):
            if dist[i] > 0:
                if dist[i] > mx:
                    mx = dist[i]
                total += dist[i]
        md[s] = mx
        sd[s] = total
    return md, sd


# ---------------------------------------------------------------- IO 與測試


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> int:
        nonlocal pos
        v = parse_int(toks[pos]) if pos < len(toks) else None
        pos += 1
        return 0 if v is None else v

    n = nxt()
    m = nxt()
    if n < 0:
        n = 0
    if m < 0:
        m = 0
    edges = [(nxt(), nxt()) for _ in range(m)]

    adj = build_adj(n, edges)
    parent, root, order = forest_data(n, adj)
    mis_size, mis_nodes = max_independent_set(n, adj, parent, order)
    dia_len, dia_path = tree_diameter(n, adj, root)
    maxdist, sumdist, sz, comp = reroot(n, adj, parent, root, order)
    cens = centroids(n, adj, parent, root, sz, comp)

    out = [
        str(mis_size),
        " ".join(str(x) for x in mis_nodes),
        str(dia_len),
        " ".join(str(x) for x in dia_path),
        str(len(cens)),
        " ".join(str(x) for x in cens),
        " ".join(str(x) for x in maxdist),
        " ".join(str(x) for x in sumdist),
    ]
    sys.stdout.write("\n".join(out) + "\n")


def run_tests() -> None:
    # ---- 固定用例：鏈 0-1-2-3 ----
    chain = [(0, 1), (1, 2), (2, 3)]
    adj = build_adj(4, chain)
    parent, root, order = forest_data(4, adj)
    size, nodes = max_independent_set(4, adj, parent, order)
    assert size == 2 and nodes == [0, 2] or size == 2 and nodes == [1, 3]
    assert mis_bruteforce(4, adj) == 2
    assert tree_diameter(4, adj, root) == (3, [0, 1, 2, 3])
    md, sd, sz, comp = reroot(4, adj, parent, root, order)
    assert md == [3, 2, 2, 3]
    assert sd == [6, 4, 4, 6]
    assert brute_distances(4, adj) == (md, sd)
    assert centroids(4, adj, parent, root, sz, comp) == [1, 2]
    assert centroids_bruteforce(4, adj, comp, root) == [1, 2]

    # ---- 固定用例：星形（中心 0，葉 1..4）----
    star = [(0, i) for i in range(1, 5)]
    adj = build_adj(5, star)
    parent, root, order = forest_data(5, adj)
    size, nodes = max_independent_set(5, adj, parent, order)
    assert size == 4 and nodes == [1, 2, 3, 4]
    assert tree_diameter(5, adj, root) == (2, [2, 0, 1])
    md, sd, sz, comp = reroot(5, adj, parent, root, order)
    assert md == [1, 2, 2, 2, 2]
    assert sd == [4, 7, 7, 7, 7]
    assert brute_distances(5, adj) == (md, sd)
    assert centroids(5, adj, parent, root, sz, comp) == [0]

    # ---- 邊界：空圖 / 單點 / 森林 ----
    adj = build_adj(0, [])
    parent, root, order = forest_data(0, adj)
    assert max_independent_set(0, adj, parent, order) == (0, [])
    assert tree_diameter(0, adj, root) == (0, [])
    assert reroot(0, adj, parent, root, order)[0] == []

    adj = build_adj(1, [])
    parent, root, order = forest_data(1, adj)
    assert max_independent_set(1, adj, parent, order) == (1, [0])
    assert tree_diameter(1, adj, root) == (0, [0])
    md, sd, sz, comp = reroot(1, adj, parent, root, order)
    assert md == [0] and sd == [0]
    assert centroids(1, adj, parent, root, sz, comp) == [0]

    # 森林：兩個不相連的鏈
    forest = [(0, 1), (2, 3)]
    adj = build_adj(4, forest)
    parent, root, order = forest_data(4, adj)
    assert max_independent_set(4, adj, parent, order)[0] == 2   # 兩條邊各取一端
    md, sd, sz, comp = reroot(4, adj, parent, root, order)
    assert md == [1, 1, 1, 1]
    assert sd == [1, 1, 1, 1]
    assert brute_distances(4, adj) == (md, sd)

    # ---- 隨機對拍：隨機樹（n ≤ 12）----
    random.seed(20261007)
    for _ in range(400):
        n = random.randint(1, 12)
        perm = list(range(n))
        random.shuffle(perm)
        edges = []
        for i in range(1, n):                       # 隨機父節點 → 保證是一棵樹
            j = random.randint(0, i - 1)
            edges.append((perm[i], perm[j]))
        adj = build_adj(n, edges)
        parent, root, order = forest_data(n, adj)

        # 最大獨立集 vs 位掩碼暴力
        size, nodes = max_independent_set(n, adj, parent, order)
        assert size == mis_bruteforce(n, adj)
        assert len(set(nodes)) == size
        for u in nodes:                             # 驗證確實是獨立集
            for v in adj[u]:
                assert v not in nodes

        # 換根 DP vs 從每點 BFS
        md, sd, sz, comp = reroot(n, adj, parent, root, order)
        bmd, bsd = brute_distances(n, adj)
        assert md == bmd and sd == bsd, (edges, md, bmd, sd, bsd)

        # 重心 vs 暴力刪點
        assert centroids(n, adj, parent, root, sz, comp) == \
            centroids_bruteforce(n, adj, comp, root)

        # 直徑：長度合理且路徑是真的簡單路徑
        dl, dp = tree_diameter(n, adj, root)
        assert dl == len(dp) - 1
        assert len(set(dp)) == len(dp)
        for i in range(len(dp) - 1):
            assert dp[i + 1] in adj[dp[i]]
        assert dl == max(bmd)

    # ---- 隨機對拍：隨機森林（含孤立點）----
    for _ in range(200):
        n = random.randint(1, 12)
        m = random.randint(0, n)
        edges = [(random.randint(0, n - 1), random.randint(0, n - 1)) for _ in range(m)]
        edges = [(u, v) for u, v in edges if u != v]
        adj = build_adj(n, edges)
        # 只保留無環的情形（隨機邊可能成環，這裡用並查集過濾）
        p = list(range(n))

        def find(x: int) -> int:
            while p[x] != x:
                p[x] = p[p[x]]
                x = p[x]
            return x

        keep = []
        for u, v in edges:
            a, b = find(u), find(v)
            if a != b:
                p[a] = b
                keep.append((u, v))
        adj = build_adj(n, keep)
        parent, root, order = forest_data(n, adj)
        md, sd, sz, comp = reroot(n, adj, parent, root, order)
        bmd, bsd = brute_distances(n, adj)
        assert md == bmd and sd == bsd
        assert max_independent_set(n, adj, parent, order)[0] == mis_bruteforce(n, adj)
        assert centroids(n, adj, parent, root, sz, comp) == \
            centroids_bruteforce(n, adj, comp, root)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
