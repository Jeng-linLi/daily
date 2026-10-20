"""Dancing Links（DLX）與精確覆蓋：Knuth Algorithm X 的雙向環狀串列實作

題意：
    給一個 0/1 矩陣，**精確覆蓋（Exact Cover）**要選出若干列（row），
    使得每一欄（column）恰好被一個被選到的列覆蓋一次。

    本程式用 Knuth 的 **Algorithm X + Dancing Links（DLX）**求解，並提供兩種模式：

        mode 1  一般精確覆蓋：給定 n × m 的 0/1 矩陣，輸出找到的解數（上限 k）
                與字典序最小的解，並與 2^n 暴力枚舉對拍（n ≤ 14 時）
        mode 2  數獨求解：把 B²×B² 的數獨轉成精確覆蓋（4·N² 個約束欄），
                N = B²，輸出填完的盤面

思路：
    ### 1. Algorithm X：遞迴 + 回溯
    - 選一個「還沒被覆蓋」的欄 c（啟發式：選**剩餘候選列最少**的那一欄，
      這是所謂的 S-heuristic，能把分支因子壓到最小）。
    - 枚舉所有在 c 欄為 1 的列 r：選 r，把 r 覆蓋到的所有欄都從矩陣中刪掉，
      遞迴；回溯時再把它們放回去。

    ### 2. Dancing Links：O(1) 的刪除與還原
    樸素做法是每次複製一份矩陣，代價 O(nm)。Knuth 的巧思是把每個 1 做成一個節點，
    用**雙向環狀串列**橫向（同一列）與縱向（同一欄）串起來，加上一個「欄表頭」，
    再串成一條含 root 的環。刪除一個欄只要：

        R[L[c]] = R[c]; L[R[c]] = L[c]            // 從欄表頭環上摘掉

    而**還原**只需要反向做一次（因為被摘掉的節點指標根本沒被清掉，
    它還記得左右的鄰居）——這就是 "dancing links" 名字的由來：

        R[L[c]] = c; L[R[c]] = c

    這個技巧讓「覆蓋 / 取消覆蓋」是 O(該欄的 1 的個數)，而不是 O(nm)，
    實務上數獨這種規模（324 欄、729 列）能在毫秒級解完。

    ### 3. 確定性：兩語言必須長出一模一樣的結構
    - 節點插入順序固定：逐列（row id 遞增）逐欄加入，同一欄內的 D 串列
      因此就是 row id 遞增的順序。
    - 選欄用 S-heuristic，平手時取**欄編號最小**的那一個（從左往右掃、
      只有嚴格小於才更新）。
    - 因為 uncover 是 cover 的精確逆運算，欄表頭環的左右順序在回溯後完全一致，
      所以「從左往右掃」的結果可重現，Python 與 C++ 能逐位元組對拍。

    ### 4. 數獨怎麼轉成精確覆蓋
    N = B²，盤面 N×N。對每個 (r, c, d) 建一列（共 N³ 列），約束欄分四組各 N² 個：
        1. 格子 (r,c) 必須填一個數      2. 第 r 列必須有數字 d
        3. 第 c 行必須有數字 d          4. 第 b 個宮必須有數字 d
    已填的格子只保留對應的 d 那一列，其餘從候選中刪掉。
    精確覆蓋的一組解 ⟺ 一組合法的數獨填法。

    ### 5. 全 0 的列要事先剔除
    一個全 0 的列覆蓋不到任何欄，選它或不選它都不影響「每欄恰好一次」，
    會讓解的數量白白乘以 2。兩個版本都在建結構時**跳過全 0 的列**，
    暴力對拍也一樣略過，兩邊語義才一致。

應用場景：
    - 數獨 / KenKen / Star Battle / Hidato 等謎題求解器（DLX 是數獨的標準解法）
    - N 皇后、八皇后變體（轉成精確覆蓋後可瞬間數出全部解）
    - 多聯骨牌（Pentomino）鋪滿矩形棋盤
    - 排班問題：護士輪班、機組指派、會議室分配（本質是「集合劃分」）
    - VLSI 測試向量最小化、覆蓋陣列（covering array）生成
    - 把 NP 問題建模成精確覆蓋後，用 DLX 當作小型的組合搜尋引擎

複雜度：
    建結構                            O(n · m)（等於矩陣中 1 的個數）
    cover / uncover 一個欄            O(該欄目前的列數)
    Algorithm X 搜尋                  指數級（精確覆蓋是 NP-complete），
                                      但 S-heuristic 在實務上極快
    空間                              O(1 的個數 + m)
"""

import random
import re
import sys
from typing import List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 40           # 矩陣列數上限
MMAX = 40           # 矩陣欄數上限
KMAX = 10           # 解數上限
BF_NMAX = 14        # 超過這個列數就不做 2^n 暴力對拍
BMAX = 4            # 數獨宮的邊長上限（B=3 → 9×9，B=4 → 16×16）


# ---------------------------------------------------------------- DLX


class DLX:
    """Knuth 的 Dancing Links：環狀雙向串列 + Algorithm X。"""

    def __init__(self, ncol: int, rows: List[Tuple[int, List[int]]]) -> None:
        # 節點 0 是 root，1..ncol 是欄表頭，之後每個 1 一個節點
        self.ncol = ncol
        self.L = list(range(-1, ncol))          # L[j] = j-1（L[0] 下面修正）
        self.R = list(range(1, ncol + 2))       # R[j] = j+1
        self.U = list(range(ncol + 1))
        self.D = list(range(ncol + 1))
        self.C = list(range(ncol + 1))          # 節點所屬的欄（表頭的 C 是自己）
        self.ROW = [-1] * (ncol + 1)
        self.size = [0] * (ncol + 1)
        self.L[0] = ncol
        self.R[ncol] = 0

        for rid, cols in rows:
            if not cols:
                continue                        # 全 0 的列直接略過
            first = -1
            for c in cols:
                idx = len(self.L)
                self.L.append(0)
                self.R.append(0)
                self.U.append(0)
                self.D.append(0)
                self.C.append(c)
                self.ROW.append(rid)
                self.size.append(0)
                # 垂直：插到欄 c 的環狀串列尾部（表頭上方）
                h = c
                self.U[idx] = self.U[h]
                self.D[idx] = h
                self.D[self.U[h]] = idx
                self.U[h] = idx
                self.size[c] += 1
                # 水平：串成這一列的環
                if first < 0:
                    first = idx
                    self.L[idx] = idx
                    self.R[idx] = idx
                else:
                    self.L[idx] = self.L[first]
                    self.R[idx] = first
                    self.R[self.L[first]] = idx
                    self.L[first] = idx

    # -- 覆蓋 / 取消覆蓋 --

    def cover(self, c: int) -> None:
        L, R, U, D, C, size = self.L, self.R, self.U, self.D, self.C, self.size
        R[L[c]] = R[c]
        L[R[c]] = L[c]
        i = D[c]
        while i != c:
            j = R[i]
            while j != i:
                U[D[j]] = U[j]
                D[U[j]] = D[j]
                size[C[j]] -= 1
                j = R[j]
            i = D[i]

    def uncover(self, c: int) -> None:
        L, R, U, D, C, size = self.L, self.R, self.U, self.D, self.C, self.size
        i = U[c]
        while i != c:
            j = L[i]
            while j != i:
                size[C[j]] += 1
                U[D[j]] = j
                D[U[j]] = j
                j = L[j]
            i = U[i]
        R[L[c]] = c
        L[R[c]] = c

    # -- Algorithm X --

    def solve(self, limit: int) -> Tuple[int, Optional[List[int]]]:
        self.limit = limit
        self.cnt = 0
        self.best: Optional[List[int]] = None
        self._search([])
        return self.cnt, self.best

    def _search(self, cur: List[int]) -> None:
        if self.cnt >= self.limit:
            return
        if self.R[0] == 0:
            self.cnt += 1
            s = sorted(cur)
            if self.best is None or s < self.best:
                self.best = s
            return
        # S-heuristic：候選列最少的欄；平手取欄編號最小（從左往右掃，嚴格小於才更新）
        c = -1
        best_sz = -1
        j = self.R[0]
        while j != 0:
            if best_sz < 0 or self.size[j] < best_sz:
                best_sz = self.size[j]
                c = j
            j = self.R[j]
        if best_sz == 0:
            return                               # 有欄湊不出來，此路不通
        self.cover(c)
        i = self.D[c]
        while i != c:
            cur.append(self.ROW[i])
            j = self.R[i]
            while j != i:
                self.cover(self.C[j])
                j = self.R[j]
            self._search(cur)
            cur.pop()
            j = self.L[i]
            while j != i:
                self.uncover(self.C[j])
                j = self.L[j]
            i = self.D[i]
            if self.cnt >= self.limit:
                break
        self.uncover(c)


# ---------------------------------------------------------------- 暴力對拍


def brute(n: int, m: int, mat: List[List[int]]) -> Tuple[int, Optional[List[int]]]:
    """2^n 枚舉所有列的子集（n ≤ 14）。全 0 的列與 DLX 一樣略過。"""
    masks = []
    ids = []
    for i in range(n):
        mk = 0
        for j in range(m):
            if mat[i][j]:
                mk |= (1 << j)
        if mk:
            masks.append(mk)
            ids.append(i)
    full = (1 << m) - 1
    cnt = 0
    best: Optional[List[int]] = None
    k = len(masks)
    for mask in range(1 << k):
        acc = 0
        ok = True
        sel: List[int] = []
        for t in range(k):
            if (mask >> t) & 1:
                if acc & masks[t]:
                    ok = False
                    break
                acc |= masks[t]
                sel.append(ids[t])
        if ok and acc == full:
            cnt += 1
            if best is None or sel < best:
                best = sel
    if cnt == 0:
        return 0, None
    return cnt, best


def solve_exact_cover(n: int, m: int, mat: List[List[int]], k: int) -> Tuple[int, Optional[List[int]]]:
    rows = []
    for i in range(n):
        # DLX 的欄編號是 1..m，節點 0 被 root 佔用
        cols = [j + 1 for j in range(m) if mat[i][j]]
        if cols:
            rows.append((i, cols))
    dlx = DLX(m, rows)
    return dlx.solve(k)


# ---------------------------------------------------------------- 數獨


def sudoku_rows(B: int, grid: List[List[int]]) -> Tuple[int, List[Tuple[int, List[int]]], List[Tuple[int, int, int]]]:
    """把數獨轉成精確覆蓋。回傳 (欄數, 候選列, 每列對應的 (r, c, d))。"""
    N = B * B
    cells = N * N

    # 欄編號 1..4·N²（節點 0 是 root），所以每組都 +1
    def cell_col(r: int, c: int) -> int:
        return r * N + c + 1

    def row_col(r: int, d: int) -> int:
        return cells + r * N + d + 1

    def col_col(c: int, d: int) -> int:
        return 2 * cells + c * N + d + 1

    def box_col(r: int, c: int, d: int) -> int:
        b = (r // B) * B + (c // B)
        return 3 * cells + b * N + d + 1

    rows: List[Tuple[int, List[int]]] = []
    meta: List[Tuple[int, int, int]] = []
    rid = 0
    for r in range(N):
        for c in range(N):
            given = grid[r][c]
            for d in range(N):
                if given != 0 and given != d + 1:
                    continue
                rows.append((rid, [cell_col(r, c), row_col(r, d), col_col(c, d), box_col(r, c, d)]))
                meta.append((r, c, d))
                rid += 1
    return 4 * cells, rows, meta


def solve_sudoku(B: int, grid: List[List[int]]) -> Optional[List[List[int]]]:
    N = B * B
    ncol, rows, meta = sudoku_rows(B, grid)
    dlx = DLX(ncol, rows)
    cnt, best = dlx.solve(1)
    if cnt == 0 or best is None:
        return None
    out = [row[:] for row in grid]
    for rid in best:
        r, c, d = meta[rid]
        out[r][c] = d + 1
    # 合法性檢查（自洽檢查，兩語言都做）
    for i in range(N):
        seen = [0] * (N + 1)
        for j in range(N):
            seen[out[i][j]] += 1
        for d in range(1, N + 1):
            assert seen[d] == 1, ("row", i, d)
        seen = [0] * (N + 1)
        for j in range(N):
            seen[out[j][i]] += 1
        for d in range(1, N + 1):
            assert seen[d] == 1, ("col", i, d)
        seen = [0] * (N + 1)
        br = (i // B) * B
        bc = (i % B) * B
        for dr in range(B):
            for dc in range(B):
                seen[out[br + dr][bc + dc]] += 1
        for d in range(1, N + 1):
            assert seen[d] == 1, ("box", i, d)
    return out


# ---------------------------------------------------------------- 輸出組裝


def run_mode1(n: int, m: int, mat: List[List[int]], k: int) -> str:
    cnt, best = solve_exact_cover(n, m, mat, k)
    if best is None:
        line1 = str(cnt)
    else:
        line1 = str(cnt) + ("".join(" " + str(v) for v in best))
    if n > BF_NMAX:
        return line1 + "\n-1\n-1\n"
    bcnt, bbest = brute(n, m, mat)
    ok = 1 if cnt == min(bcnt, k) else 0
    if bcnt <= k and ok:
        ok = 1 if best == bbest else 0
    return line1 + "\n" + str(bcnt) + "\n" + str(ok) + "\n"


def run_mode2(B: int, grid: List[List[int]]) -> str:
    sol = solve_sudoku(B, grid)
    N = B * B
    if sol is None:
        return "no solution\n"
    out = []
    for r in range(N):
        out.append(" ".join(str(sol[r][c]) for c in range(N)))
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例 1：經典的小精確覆蓋（Knuth 論文裡的例子）----
    # 欄 A B C D E F G 共 7 欄；唯一解是 rows 0, 3, 5
    mat1 = [
        [0, 0, 1, 0, 1, 1, 0],
        [1, 0, 0, 1, 0, 0, 1],
        [0, 1, 1, 0, 0, 1, 0],
        [1, 0, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 1],
        [0, 0, 0, 1, 1, 0, 1],
    ]
    out = run_mode1(6, 7, mat1, 10)
    assert out.split("\n")[0] == "1 0 3 4", out
    assert out.split("\n")[1] == "1", out
    assert out.split("\n")[2] == "1", out

    # ---- 固定用例 2：無解（某一欄全 0）----
    mat2 = [[1, 0], [0, 0], [1, 0]]
    out = run_mode1(3, 2, mat2, 10)
    assert out.split("\n")[0] == "0", out
    assert out.split("\n")[1] == "0", out
    assert out.split("\n")[2] == "1", out

    # ---- 固定用例 3：欄數為 0（空集合也是一個解）----
    out = run_mode1(3, 0, [[], [], []], 10)
    assert out.split("\n")[0] == "1", out
    assert out.split("\n")[1] == "1", out

    # ---- 固定用例 4：三解，驗證解數上限 k ----
    mat4 = [[1, 0], [0, 1], [1, 0], [0, 1]]
    out = run_mode1(4, 2, mat4, 10)
    assert out.split("\n")[1] == "4", out        # 2 × 2 = 4 組
    out = run_mode1(4, 2, mat4, 2)
    assert out.split("\n")[0] == "2 0 1", out

    # ---- 數獨：經典 9×9 ----
    puzzle = [
        [5, 3, 0, 0, 7, 0, 0, 0, 0],
        [6, 0, 0, 1, 9, 5, 0, 0, 0],
        [0, 9, 8, 0, 0, 0, 0, 6, 0],
        [8, 0, 0, 0, 6, 0, 0, 0, 3],
        [4, 0, 0, 8, 0, 3, 0, 0, 1],
        [7, 0, 0, 0, 2, 0, 0, 0, 6],
        [0, 6, 0, 0, 0, 0, 2, 8, 0],
        [0, 0, 0, 4, 1, 9, 0, 0, 5],
        [0, 0, 0, 0, 8, 0, 0, 7, 9],
    ]
    txt = run_mode2(3, puzzle)
    lines = txt.rstrip("\n").split("\n")
    assert len(lines) == 9
    assert lines[0] == "5 3 4 6 7 8 9 1 2", lines[0]
    assert lines[8] == "3 4 5 2 8 6 1 7 9", lines[8]

    # ---- 數獨：4×4 ----
    p4 = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    txt4 = run_mode2(2, p4)
    assert txt4.rstrip("\n").split("\n")[0] == "1 2 3 4", txt4

    # ---- 數獨：無解 ----
    bad = [[1, 2, 3, 4], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    for c in range(4):
        bad[1][c] = 0
    bad[1][0] = 3                     # 第一宮有兩個 3 的需求衝突（其實仍可能無解）
    bad[1][1] = 1                     # (1,0)=3 與 (0,0)=1，(1,1)=1 同一列重複 1
    txtb = run_mode2(2, bad)
    assert txtb == "no solution\n", txtb

    # ---- 隨機對拍：DLX vs 2^n 暴力 ----
    rng = random.Random(20261020)
    for t in range(600):
        n = rng.randint(0, 10)
        m = rng.randint(0, 8)
        p = rng.choice([0.15, 0.3, 0.5, 0.7])
        mat = [[1 if rng.random() < p else 0 for _ in range(m)] for _ in range(n)]
        out = run_mode1(n, m, mat, 10)
        assert out.split("\n")[2] == "1", (t, n, m, p, out)

    # ---- 隨機對拍：稍微大一點的矩陣（n ≤ 14）----
    for t in range(200):
        n = rng.randint(8, 14)
        m = rng.randint(4, 10)
        p = rng.choice([0.2, 0.35, 0.5])
        mat = [[1 if rng.random() < p else 0 for _ in range(m)] for _ in range(n)]
        out = run_mode1(n, m, mat, 5)
        assert out.split("\n")[2] == "1", (t, n, m, p, out)

    # ---- 隨機數獨（挖空 9×9）----
    for t in range(6):
        rg = random.Random(3000 + t)
        base = solve_sudoku(3, [[0] * 9 for _ in range(9)])
        assert base is not None
        holes = rg.sample(range(81), 55)
        g = [row[:] for row in base]
        for h in holes:
            g[h // 9][h % 9] = 0
        sol = solve_sudoku(3, g)
        assert sol is not None, t
        # 挖空盤的解必須與原盤在已給定的格子上一致
        for r in range(9):
            for c in range(9):
                if g[r][c] != 0:
                    assert sol[r][c] == g[r][c], (t, r, c)


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

    mode = next_int(1)
    if mode == 2:
        B = next_int(3)
        if B < 1:
            B = 1
        if B > BMAX:
            B = BMAX
        N = B * B
        grid = [[0] * N for _ in range(N)]
        for r in range(N):
            for c in range(N):
                v = next_int(0)
                if v < 0:
                    v = 0
                if v > N:
                    v = N
                grid[r][c] = v
        emit(run_mode2(B, grid))
        return

    n = next_int(0)
    m = next_int(0)
    if n < 0:
        n = 0
    if m < 0:
        m = 0
    if n > NMAX:
        n = NMAX
    if m > MMAX:
        m = MMAX
    mat = []
    for _ in range(n):
        row = []
        for _ in range(m):
            v = next_int(0)
            row.append(1 if v != 0 else 0)
        mat.append(row)
    k = next_int(1)
    if k < 0:
        k = 0
    if k > KMAX:
        k = KMAX
    emit(run_mode1(n, m, mat, k))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
