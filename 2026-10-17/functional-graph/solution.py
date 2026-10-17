"""基環樹與函數圖（Functional Graph / Pseudo-Forest）

題意：
    本題分兩部分，兩部分都同時給出「高效演算法」與「獨立實作的暴力版本」並互相對拍。

    A. 有向函數圖：n 個點，每個點 i 恰有一條出邊 nxt[i]（0 ≤ nxt[i] < n）。
       1. 找出所有環，每個環正規化成「從環上最小點號開始」的順序
       2. belong[i]：i 沿出邊一直走，最終進入的那個環的代表（環上最小點號）
       3. dist_to_cycle[i]：i 走到環上所需步數（環上點為 0）
       4. cycle_len[i]：i 最終進入的那個環的長度
       5. jump(i, K)：從 i 出發走 K 步後落在哪個點（K 可達 10^18，用二進制跳躍）

    B. 無向基環樹森林：n2 個點、m2 條邊，每個連通分量恰有一個環（退化成樹時環長為 0）。
       對每個分量輸出：最小點號、環長、直徑（任意兩點間最短路徑的最大值）。

思路：
    ### A1. 找環：沿著唯一出邊走，一定會「走進一條尾巴，再進入一個環」
    函數圖的每個連通分量恰有一個環，形狀是「很多棵內向樹掛在環上」。
    用 0/1/2 三色標記（0 未訪問、1 在當前這條鏈上、2 已處理完）：
      - 從未訪問點 s 出發一直走，把沿途節點記進 path 與 pos（節點 → 在 path 中的下標）。
      - 若撞到「顏色 1」的點 v：path[pos[v]:] 就是那個環（因為 v 在當前鏈上出現了兩次）。
      - 若撞到「顏色 2」的點 v：說明這條尾巴接到了別人已經處理好的環上。
    環上的點距離為 0；尾巴上的點用**倒序**傳遞：
        dist_to_cycle[u] = dist_to_cycle[nxt[u]] + 1
    倒著走 path 就能保證兒子先算完，一次搞定，總共 O(n)。

    ### A2. 走 K 步：二進制跳躍（倍增）
    up[j][i] 表示從 i 出發走 2^j 步到達的點；up[0] = nxt，
    up[j][i] = up[j-1][ up[j-1][i] ]。把 K 拆成二進位，從低位往高位跳即可。
    複雜度 O(n log K) 預處理、O(log K) 單次查詢。

    ### B1. 找環：拓撲剝葉（砍掉所有度數 ≤ 1 的點）
    把度數 ≤ 1 的點丟進佇列，反覆取出、從圖上刪掉、並把它鄰居的度數減 1，
    新產生度數 ≤ 1 的點再入隊。最後**沒被刪掉的點**就是 2-core；
    在基環樹裡，2-core 恰好就是那個唯一的環（樹的話會被全部刪光，環長 0）。
    剝葉的同時維護 dep[u] = u 下方（已被刪掉那一側）子樹的最大深度，
    以及 top1/top2 = u 各個分支深度的前兩大：
        「經過 u 的最長路徑」= top1[u] + top2[u]
    這是樹形 DP 求直徑的標準寫法。

    ### B2. 基環樹直徑 = max(每棵掛在環上的樹的直徑, 穿過環的最長路徑)
    設環上按順序為 c_0..c_{L-1}，環上點 c_i 掛著深度 val[i] = dep[c_i] 的樹。
    環上兩點 c_i、c_j 的距離是 min(|i-j|, L-|i-j|)（單位邊權），
    所以候選答案是 val[i] + val[j] + min(|i-j|, L-|i-j|)。
    把 val 複製一份接在後面（破環成鏈），枚舉 i、用**滑動窗口最大值**
    在 O(1) 內拿到「另一端的最大值」，整條式子 O(L) 就算完：
      - 若 d = j - i ≤ L/2：距離 = d，式子寫成 (val[j] + j) + (val[i] - i)
      - 若 d > L/2：距離 = L - d，式子寫成 (val[j] - j) + (val[i] + i + L)
    兩種情況分別對 (val[j]+j) 與 (val[j]-j) 做滑動窗口最大值即可。

應用場景：
    - 有向函數圖：置換群的循環分解、追蹤「不斷套用某個映射」的收斂行為、
      雜湊碰撞（生日攻擊）裡 rho 形的路徑、Pollard-Rho 因數分解的形狀本質上就是函數圖
    - 二進制跳躍：置換的 k 次冪、自動機跑 k 步的落點、連結串列/指標跳躍、LCA 倍增同理
    - 無向基環樹：n 點 n 邊的圖問題（城市的環狀道路 + 支線）、
      「Island / 基環樹森林直徑」這類經典題、網路拓撲裡單環骨架的直徑估計

複雜度：
    A 找環                O(n) 時間 / O(n) 空間
    A 二進制跳躍          O(n log K) 預處理，O(log K) 單次查詢
    B 剝葉 + 環上 DP      O(n + m + L) 時間 / O(n + m) 空間（L 為環長）
"""

import random
import re
import sys
from collections import deque
from typing import List, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 60        # 有向函數圖在 IO 模式下的最大點數
MAXN2 = 30       # 無向圖最大點數
MAXM2 = 200      # 無向圖最大邊數
KMAX = 10 ** 18  # 跳躍步數上限
K_BRUTE = 3000   # 超過這個步數就不跑暴力模擬（只比對倍增結果本身）
LOG = 62         # 2^62 > KMAX


# ---------------------------------------------------------------- A. 有向函數圖


def analyze_functional(nxt: List[int]) -> Tuple[List[int], List[int], List[int], List[List[int]]]:
    """回傳 (belong, dist_to_cycle, cycle_len, cycles)。cycles 已按最小點號排序。"""
    n = len(nxt)
    belong = [-1] * n
    dtc = [0] * n
    clen = [0] * n
    state = [0] * n          # 0 未訪問 / 1 在當前鏈上 / 2 已處理
    cycles: List[List[int]] = []

    for s in range(n):
        if state[s] != 0:
            continue
        path: List[int] = []
        pos = {}
        v = s
        while state[v] == 0:
            state[v] = 1
            pos[v] = len(path)
            path.append(v)
            v = nxt[v]

        if state[v] == 1:                     # 撞回當前鏈 → 找到新環
            cyc = path[pos[v]:]
            L = len(cyc)
            mn = min(cyc)
            k = cyc.index(mn)
            cycles.append(cyc[k:] + cyc[:k])  # 正規化：從最小點號開始
            for u in cyc:
                belong[u] = mn
                dtc[u] = 0
                clen[u] = L
                state[u] = 2
            start = pos[v] - 1
        else:                                 # 接到已處理好的環上
            start = len(path) - 1

        # 尾巴上的點倒序傳遞（兒子已經算好了）
        for j in range(start, -1, -1):
            w = nxt[path[j]]
            belong[path[j]] = belong[w]
            dtc[path[j]] = dtc[w] + 1
            clen[path[j]] = clen[w]
            state[path[j]] = 2

    cycles.sort(key=lambda c: c[0])
    return belong, dtc, clen, cycles


def build_lift(nxt: List[int]) -> List[List[int]]:
    """up[j][i] = 從 i 出發走 2^j 步的落點。"""
    n = len(nxt)
    up = [list(nxt)]
    for _ in range(1, LOG):
        prev = up[-1]
        up.append([prev[prev[i]] for i in range(n)])
    return up


def jump_all(nxt: List[int], k: int) -> List[int]:
    """從每個點出發走 k 步的落點。"""
    n = len(nxt)
    k = max(0, min(k, KMAX))
    if n == 0:
        return []
    up = build_lift(nxt)
    res = []
    for i in range(n):
        v = i
        x = k
        b = 0
        while x:
            if x & 1:
                v = up[b][v]
            x >>= 1
            b += 1
        res.append(v)
    return res


def jump_brute(nxt: List[int], k: int) -> List[int]:
    """暴力：老實一步一步走。"""
    res = []
    for i in range(len(nxt)):
        v = i
        for _ in range(k):
            v = nxt[v]
        res.append(v)
    return res


def functional_brute(nxt: List[int]) -> Tuple[List[int], List[int], List[int], List[Tuple[int, int]]]:
    """暴力：對每個點模擬走 n 步，記錄第一次撞見重複點的位置。"""
    n = len(nxt)
    belong, dtc, clen = [], [], []
    pairs = set()
    for i in range(n):
        seen = {}
        v = i
        step = 0
        while v not in seen:
            seen[v] = step
            v = nxt[v]
            step += 1
        first = seen[v]                                  # 進入環的步數
        cyc = [u for u, s in seen.items() if s >= first]
        L = step - first
        mn = min(cyc)
        belong.append(mn)
        dtc.append(first)
        clen.append(L)
        pairs.add((mn, L))
    return belong, dtc, clen, sorted(pairs)


# ---------------------------------------------------------------- B. 無向基環樹


def sliding_max(arr: List[int], size: int) -> List[int]:
    """out[k] = max(arr[k..k+size-1])，單調佇列 O(len(arr))。"""
    n = len(arr)
    if size <= 0 or size > n:
        return []
    out = []
    dq = deque()                       # 存下標，對應值遞減
    for i in range(n):
        while dq and arr[dq[-1]] <= arr[i]:
            dq.pop()
        dq.append(i)
        while dq[0] <= i - size:
            dq.popleft()
        if i >= size - 1:
            out.append(arr[dq[0]])
    return out


def ring_diameter_fast(vals: List[int]) -> int:
    """環上兩點 i<j，答案候選 vals[i] + vals[j] + min(d, L-d)。滑動窗口 O(L)。"""
    L = len(vals)
    if L <= 1:
        return 0
    v2 = vals + vals
    half = L // 2
    key_a = [v2[j] + j for j in range(2 * L)]
    win_a = sliding_max(key_a, half)                 # 窗口 [i+1, i+half]
    size_b = L - half - 1
    key_b = [v2[j] - j for j in range(2 * L)]
    win_b = sliding_max(key_b, size_b) if size_b > 0 else []
    best = 0
    for i in range(L):
        if half > 0:
            k = i + 1
            if k < len(win_a):
                best = max(best, vals[i] - i + win_a[k])
        if size_b > 0:
            k = i + half + 1
            if k < len(win_b):
                best = max(best, vals[i] + i + L + win_b[k])
    return best


def ring_diameter_brute(vals: List[int]) -> int:
    L = len(vals)
    best = 0
    for i in range(L):
        for j in range(i + 1, L):
            d = j - i
            best = max(best, vals[i] + vals[j] + min(d, L - d))
    return best


def ring_order(adj: List[List[int]], on: List[bool], nodes: List[int]) -> List[int]:
    """把環上的點排成一圈（從最小點號出發）。"""
    if not nodes:
        return []
    start = min(nodes)
    order = [start]
    prev = -1
    cur = start
    while len(order) < len(nodes):
        cand = [w for w in adj[cur] if on[w] and w != prev]
        if not cand:
            break
        nxt = min(cand)
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def analyze_pseudo_forest(n: int, edges: List[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    """回傳每個連通分量的 (最小點號, 環長, 直徑)，按最小點號排序。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    deg = [len(adj[i]) for i in range(n)]

    on = [True] * n
    dep = [0] * n            # 往「已被剝掉」那一側的最大深度
    top1 = [0] * n           # 各分支深度的最大值
    top2 = [0] * n           # 次大值
    q = deque(i for i in range(n) if deg[i] <= 1)
    while q:
        u = q.popleft()
        if not on[u]:
            continue
        on[u] = False
        for v in adj[u]:
            if not on[v]:
                continue
            branch = dep[u] + 1
            if branch > top1[v]:
                top2[v] = top1[v]
                top1[v] = branch
            elif branch > top2[v]:
                top2[v] = branch
            if branch > dep[v]:
                dep[v] = branch
            deg[v] -= 1
            if deg[v] == 1:
                q.append(v)

    comp = [-1] * n
    comps: List[List[int]] = []
    for s in range(n):
        if comp[s] != -1:
            continue
        cid = len(comps)
        comp[s] = cid
        stack = [s]
        nodes: List[int] = []
        while stack:
            x = stack.pop()
            nodes.append(x)
            for w in adj[x]:
                if comp[w] == -1:
                    comp[w] = cid
                    stack.append(w)
        comps.append(nodes)

    out = []
    for nodes in comps:
        cyc = [x for x in nodes if on[x]]
        inner = max(top1[x] + top2[x] for x in nodes)
        if len(cyc) <= 1:
            ring = 0
        else:
            order = ring_order(adj, on, cyc)
            ring = ring_diameter_fast([dep[x] for x in order])
        out.append((min(nodes), len(cyc), max(inner, ring)))
    out.sort(key=lambda t: t[0])
    return out


def pseudo_forest_brute(n: int, edges: List[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    """暴力：連通分量用 BFS；環長 = 落在某個環上的點數；直徑 = 全點對 BFS。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    comp = [-1] * n
    comps: List[List[int]] = []
    for s in range(n):
        if comp[s] != -1:
            continue
        cid = len(comps)
        comp[s] = cid
        nodes = [s]
        i = 0
        while i < len(nodes):
            x = nodes[i]
            i += 1
            for w in adj[x]:
                if comp[w] == -1:
                    comp[w] = cid
                    nodes.append(w)
        comps.append(nodes)

    def reachable(x: int, y: int, ban_u: int, ban_v: int) -> bool:
        """x 能否在不走 (ban_u, ban_v) 這條邊的情況下到達 y。"""
        seen = [False] * n
        seen[x] = True
        stack = [x]
        while stack:
            a = stack.pop()
            if a == y:
                return True
            for b in adj[a]:
                if (a == ban_u and b == ban_v) or (a == ban_v and b == ban_u):
                    continue
                if not seen[b]:
                    seen[b] = True
                    stack.append(b)
        return False

    out = []
    for nodes in comps:
        cyc_cnt = 0
        for v in nodes:
            on_cycle = False
            for u in set(adj[v]):
                if u == v:                       # 自環
                    on_cycle = True
                    break
                if reachable(u, v, u, v):        # 繞別條路還能回來 → v 在環上
                    on_cycle = True
                    break
            if on_cycle:
                cyc_cnt += 1
        dia = 0
        for s in nodes:
            dist = [-1] * n
            dist[s] = 0
            stack = [s]
            i = 0
            while i < len(stack):
                x = stack[i]
                i += 1
                for w in adj[x]:
                    if dist[w] == -1:
                        dist[w] = dist[x] + 1
                        stack.append(w)
            for t in nodes:
                if dist[t] > dia:
                    dia = dist[t]
        out.append((min(nodes), cyc_cnt, dia))
    out.sort(key=lambda t: t[0])
    return out


# ---------------------------------------------------------------- 輸出


def join(nums: List[int]) -> str:
    return " ".join(str(v) for v in nums)


def run(nxt: List[int], k: int, n2: int, edges: List[Tuple[int, int]]) -> str:
    belong, dtc, clen, cycles = analyze_functional(nxt)
    b_belong, b_dtc, b_clen, b_pairs = functional_brute(nxt)
    jumps = jump_all(nxt, k)

    ok = 1
    if belong != b_belong or dtc != b_dtc or clen != b_clen:
        ok = 0
    algo_pairs = sorted((c[0], len(c)) for c in cycles)
    if algo_pairs != b_pairs:
        ok = 0
    if k <= K_BRUTE and jumps != jump_brute(nxt, k):
        ok = 0
    for c in cycles:                                  # 環的形狀自洽性
        if len(set(c)) != len(c) or min(c) != c[0]:
            ok = 0
        for idx, u in enumerate(c):
            if nxt[u] != c[(idx + 1) % len(c)]:
                ok = 0

    info = analyze_pseudo_forest(n2, edges)
    binfo = pseudo_forest_brute(n2, edges)
    if info != binfo:
        ok = 0
    lines = [
        str(len(cycles)),
        " ".join("%d:%d" % (c[0], len(c)) for c in cycles),
        join(belong),
        join(dtc),
        join(clen),
        join(jumps),
        str(len(info)),
        " ".join("%d:%d:%d" % t for t in info),
        str(ok),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 內置測試


def gen_pseudo_forest(rng: random.Random, n: int) -> List[Tuple[int, int]]:
    """產生基環樹森林：隨機森林 + 每個分量最多加一條弦。"""
    edges: List[Tuple[int, int]] = []
    eset = set()
    parent = list(range(n))
    for i in range(1, n):
        if rng.random() < 0.85:
            p = rng.randrange(0, i)
            e = (p, i)
            edges.append(e)
            eset.add(e)
            eset.add((i, p))
    # 連通分量（並查集）
    uf = list(range(n))

    def find(x: int) -> int:
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x

    for u, v in edges:
        a, b = find(u), find(v)
        if a != b:
            uf[a] = b
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    for members in groups.values():
        if len(members) < 2 or rng.random() < 0.25:
            continue                                  # 保持是樹
        for _ in range(20):
            u = rng.choice(members)
            v = rng.choice(members)
            if u != v and (u, v) not in eset:
                edges.append((u, v))
                eset.add((u, v))
                eset.add((v, u))
                break
    return edges


def run_tests() -> None:
    # ---- 固定用例：有向函數圖 ----
    nxt = [1, 2, 0, 1]
    belong, dtc, clen, cycles = analyze_functional(nxt)
    assert cycles == [[0, 1, 2]]
    assert belong == [0, 0, 0, 0]
    assert dtc == [0, 0, 0, 1]
    assert clen == [3, 3, 3, 3]
    assert jump_all(nxt, 0) == [0, 1, 2, 3]
    assert jump_all(nxt, 1) == [1, 2, 0, 1]
    assert jump_all(nxt, 100) == jump_brute(nxt, 100)

    nxt2 = [1, 2, 3, 4, 2]           # 0→1→2→3→4→2，環 {2,3,4}
    belong2, dtc2, clen2, cy2 = analyze_functional(nxt2)
    assert cy2 == [[2, 3, 4]]
    assert dtc2 == [2, 1, 0, 0, 0]
    assert clen2 == [3, 3, 3, 3, 3]
    assert belong2 == [2] * 5

    assert analyze_functional([]) == ([], [], [], [])
    assert jump_all([], 5) == []
    assert analyze_functional([0])[3] == [[0]]
    assert jump_all([0], 10 ** 18) == [0]

    # 純大環：jump 應等價於模運算
    big = [(i + 1) % 50 for i in range(50)]
    assert jump_all(big, 10 ** 18) == [(i + 10 ** 18) % 50 for i in range(50)]

    # ---- 固定用例：無向基環樹 ----
    # 三角形 0-1-2-0，0 掛 3、3 掛 4
    e1 = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4)]
    info = analyze_pseudo_forest(5, e1)
    assert info == [(0, 3, 3)], info            # 4→3→0→2 距離 3
    assert info == pseudo_forest_brute(5, e1)

    # 純樹：鏈 0-1-2-3
    e2 = [(0, 1), (1, 2), (2, 3)]
    assert analyze_pseudo_forest(4, e2) == [(0, 0, 3)]
    assert analyze_pseudo_forest(4, e2) == pseudo_forest_brute(4, e2)

    # 孤立點
    assert analyze_pseudo_forest(1, []) == [(0, 0, 0)]
    assert analyze_pseudo_forest(0, []) == []

    # 四元環
    e3 = [(0, 1), (1, 2), (2, 3), (3, 0)]
    assert analyze_pseudo_forest(4, e3) == [(0, 4, 2)]
    assert analyze_pseudo_forest(4, e3) == pseudo_forest_brute(4, e3)

    # ---- 環上 DP：滑動窗口版 vs O(L²) 暴力 ----
    rng0 = random.Random(1)
    for _ in range(300):
        L = rng0.randint(1, 12)
        vals = [rng0.randint(0, 9) for _ in range(L)]
        assert ring_diameter_fast(vals) == ring_diameter_brute(vals)

    # ---- 隨機對拍 ----
    rng = random.Random(20261017)
    for _ in range(400):
        n = rng.randint(0, 20)
        nxt = [rng.randrange(0, n) for _ in range(n)] if n else []
        k = rng.choice([0, 1, 2, 3, rng.randint(0, 60), rng.randint(0, 2000)])
        n2 = rng.randint(0, 10)
        edges = gen_pseudo_forest(rng, n2) if n2 else []
        out = run(nxt, k, n2, edges)
        assert out.split("\n")[-1] == "1", (nxt, k, n2, edges, out)

    # ---- 長鏈只跑線性演算法 ----
    chain = [i + 1 for i in range(1999)] + [1999]
    belong3, dtc3, clen3, cy3 = analyze_functional(chain)
    assert cy3 == [[1999]]
    assert dtc3[0] == 1999
    assert clen3[0] == 1
    assert jump_all(chain, 10 ** 18)[0] == 1999


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
        return max(-(10 ** 19), min(10 ** 19, v))

    n = next_int(0)
    n = max(0, min(n, MAXN))
    nxt = []
    for _ in range(n):
        v = next_int(0)
        if v < 0 or v >= n:
            v = 0
        nxt.append(v)
    k = max(0, next_int(0))
    n2 = max(0, min(next_int(0), MAXN2))
    m2 = max(0, min(next_int(0), MAXM2))
    edges = []
    for _ in range(m2):
        u = next_int(0)
        v = next_int(0)
        if 0 <= u < n2 and 0 <= v < n2:
            edges.append((u, v))
    emit(run(nxt, k, n2, edges) + "\n")


def main() -> None:
    sys.setrecursionlimit(1000000)
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
