"""高斯消元：GF(p) 上的線性方程組 + GF(2) 上的 XOR 線性基

題意：
    輸入一個 n×m 的整數矩陣 A、右端向量 b、一個質數模數 p，
    以及一組非負整數 X 與一個整數 k。輸出：
      1. A 在模 p 下的秩與行列式；
      2. 線性方程組 A·x ≡ b (mod p) 的解的個數分類與一組解；
      3. A 在模 p 下是否可逆；
      4. X 的 XOR 線性基大小，以及 X 能組出的最大異或值；
      5. X 能組出的值之中，第 k 小的是多少。

思路：
    ### 高斯消元到底在做什麼
    線性方程組有三種基本操作不改動解集合：
      ① 交換兩列（行列式變號）  ② 一列乘上非零常數  ③ 一列加上另一列的倍數
    高斯消元就是用這三種操作把矩陣化成**列階梯型**，再化成**簡化列階梯型（RREF）**：
    每個主元所在行只有該主元非零，於是可以「一眼讀出解」。

    ### 為什麼本題用模質數域 GF(p)
    在實數上做高斯消元會遇到浮點誤差，Python 與 C++ 的輸出不可能逐位元組一致。
    改成在**模質數 p** 的有限域上運算後，所有運算都是精確整數：
      - 除法 → 乘上乘法逆元；p 是質數時，費馬小定理給出 `a⁻¹ ≡ a^(p-2) (mod p)`；
      - 每一個非零元素都可逆，所以主元搜尋只要找「非零」即可，不存在「p 是合數時無逆元」的麻煩。
    這讓兩個語言的輸出**必然一致**，也讓行列式、秩、可逆性都有精確定義。

    ### 行列式怎麼從消元過程讀出來
    記錄每次列交換的次數 s 與每個主元 pᵢ（正規化之前的值）。
    三個操作對行列式的影響分別是 ×(-1)、×(1/pᵢ)、×1，而最後得到 RREF（可逆時為單位矩陣，行列式 1），
    所以 `det(A) = (-1)^s · Π pᵢ (mod p)`。若消元結束時秩 < n（方陣不滿秩），行列式為 0。

    ### 解的分類
    把增廣矩陣 [A | b] 一起消元，得到 rank(A) 與 rank([A|b])：
      - `rank(A) < rank([A|b])` → 出現「0 = 非零」的列 → **無解**；
      - `rank(A) = m`（未知數個數）→ **唯一解**；
      - 否則 → **無窮多解**，自由變數個數 `m - rank`，在 GF(p) 上共有 `p^(m-rank)` 組解。
    本題把自由變數一律補 0，保證輸出唯一（兩個語言才對得齊）。

    ### XOR 線性基：GF(2) 上的同一套高斯消元
    把每個數字看成 63 維的 GF(2) 向量，異或就是加法。
    「線性基」就是這組向量的一組基底，`basis[i]` 是最高位為 i 的基底向量。
    插入時從高位往低位消：`x` 的第 i 位是 1，若 `basis[i]` 為空就放進去，否則 `x ^= basis[i]` 繼續往下。
    這跟高斯消元的「拿主元去消掉其他列」完全同構。

    基底建好後：
      - **最大異或值**：從高位貪心，`res = max(res, res ^ basis[i])`；
      - **第 k 小**：先把基底化成「每個 basis[i] 的低位都不含其它基底的最高位」
        （`for i: for j<i: if basis[i] 的第 j 位為 1: basis[i] ^= basis[j]`），
        再把非零基底由小到大排成 d，則第 k 小（1-indexed）就是
        把 `k-1` 的二進位位元當成選擇向量：`ans ^= d[i]`（第 i 位為 1 時）。
        因為這恰好是前綴樹的前序走訪，順序就是數值大小順序。

應用場景：
    - 電路模擬（SPICE /  nodal analysis）與結構力學的線性方程組求解
    - 電腦圖學的座標變換與矩陣求逆、機器學習的最小平方迴歸正規方程
    - 密碼學：GF(2) 上的線性分析、MDS 矩陣、糾刪碼（Reed-Solomon 的編碼矩陣）
    - 競賽與遊戲：Nim 變體的必勝判定、最大異或路徑、子集異或統計
    - 資料去重與特徵壓縮：判斷一組特徵向量的張成空間、求極大線性無關組

複雜度：
    高斯消元（RREF）  O(n · m · min(n,m))，本題 n、m ≤ 40
    行列式 / 可逆性   隨消元一起得到，不額外花時間
    XOR 線性基插入    O(63) / 個數字
    最大異或值        O(63)
    第 k 小           O(63²) 的基底重整 + O(63)
    空間              O(n·m)

輸入格式（stdin，全部以空白分隔）：
    n m p
    a_11 ... a_1m        n 列，每列 m 個整數
    ...
    b_1 ... b_n          右端向量
    c                    參與線性基的數字個數
    x_1 ... x_c          這些數字（先截斷到 63 位元無號）
    k                    第 k 小異或值（1-indexed）

    非數字 token 一律視為「輸入到此為止」；讀不滿的那一列/那一組直接丟棄。
    數值截斷：n、m、c ≤ 40；p 會被提升到「不小於輸入值的最小質數」（p < 2 時取 2）。
輸出格式（stdout，共 6 行）：
    第 1 行  rank det              秩；行列式（n == m 時有意義，否則記 0）
    第 2 行  solkind free_cnt      0 = 無解、1 = 唯一解、2 = 無窮多解（後者給自由變數數）
    第 3 行  m 個整數               一組解（無解時為空行；自由變數補 0）
    第 4 行  invertible            A 在模 p 下是否可逆（1 / 0）
    第 5 行  basis_size max_xor    線性基大小、能組出的最大異或值
    第 6 行  kth_xor               第 k 小可組出的異或值（k 越界為 -1）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import itertools
import random
import re
import sys
from typing import List, Optional, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

MAXNM = 40
P_MAX = 2147483647          # 2^31 - 1，本身是質數
MAXB = 62                   # 線性基最高位（位元 0..62）
MASK = (1 << 63) - 1


# ---------------------------------------------------------------- 質數處理


def is_prime(v: int) -> bool:
    """試除法判定質數（兩個語言一致）。"""
    if v < 2:
        return False
    if v % 2 == 0:
        return v == 2
    d = 3
    while d * d <= v:
        if v % d == 0:
            return False
        d += 2
    return True


def to_prime(v: int) -> int:
    """把模數提升到「不小於 v 的最小質數」，確保 GF(p) 是個域。"""
    if v < 2:
        return 2
    while not is_prime(v):
        v += 1
    return v


def mod_inv(a: int, p: int) -> int:
    """費馬小定理求逆元：p 為質數且 a ≢ 0 (mod p)。"""
    return pow(a % p, p - 2, p)


# ---------------------------------------------------------------- 高斯消元


def rref(a: Sequence[Sequence[int]], b: Sequence[int], n: int, m: int,
         p: int) -> Tuple[int, int, List[List[int]]]:
    """對增廣矩陣 [A | b] 做簡化列階梯化。

    返回 (rank, det, aug)，其中 det 只有在 n == m 時有意義（否則回傳 0）。
    """
    aug = [[a[i][j] % p for j in range(m)] + [b[i] % p] for i in range(n)]
    rank = 0
    det = 1
    sign = 1
    col = 0
    while col < m and rank < n:
        piv = -1
        for i in range(rank, n):
            if aug[i][col] % p != 0:
                piv = i
                break
        if piv < 0:
            col += 1
            continue
        if piv != rank:
            aug[piv], aug[rank] = aug[rank], aug[piv]
            sign = -sign
        pivot = aug[rank][col] % p
        det = (det * pivot) % p
        inv = mod_inv(pivot, p)
        row = aug[rank]
        aug[rank] = [(v * inv) % p for v in row]
        row = aug[rank]
        for i in range(n):
            if i != rank and aug[i][col] % p != 0:
                f = aug[i][col] % p
                ai = aug[i]
                aug[i] = [(ai[j] - f * row[j]) % p for j in range(m + 1)]
        rank += 1
        col += 1

    if n == m:
        det = (sign * det) % p if rank == n else 0
    else:
        det = 0
    return rank, det, aug


def solve_linear(a: Sequence[Sequence[int]], b: Sequence[int], n: int, m: int,
                 p: int) -> Tuple[int, int, int, int, List[int]]:
    """解 A·x ≡ b (mod p)。

    返回 (rank, det, solkind, free_cnt, x)。
    solkind：0 = 無解、1 = 唯一解、2 = 無窮多解。
    """
    rank, det, aug = rref(a, b, n, m, p)
    for i in range(rank, n):
        if aug[i][m] % p != 0:
            return rank, det, 0, m - rank, []          # 0 = 非零 → 無解

    x = [0] * m
    which_col: List[int] = []
    col = 0
    r = 0
    while col < m and r < rank:                        # 重新定位主元欄
        while col < m and aug[r][col] % p == 0:
            col += 1
        if col >= m:
            break
        which_col.append(col)
        r += 1
        col += 1
    for i, c in enumerate(which_col):
        x[c] = aug[i][m] % p

    if rank == m:
        return rank, det, 1, 0, x
    return rank, det, 2, m - rank, x


def matrix_inverse(a: Sequence[Sequence[int]], n: int,
                   p: int) -> Optional[List[List[int]]]:
    """在模 p 下求逆矩陣；不可逆回傳 None。把 [A | I] 一起消元即可。"""
    aug = [[a[i][j] % p for j in range(n)] + [1 if i == j else 0 for j in range(n)]
           for i in range(n)]
    for col in range(n):
        piv = -1
        for i in range(col, n):
            if aug[i][col] % p != 0:
                piv = i
                break
        if piv < 0:
            return None
        aug[col], aug[piv] = aug[piv], aug[col]
        inv = mod_inv(aug[col][col] % p, p)
        aug[col] = [(v * inv) % p for v in aug[col]]
        for i in range(n):
            if i != col and aug[i][col] % p != 0:
                f = aug[i][col] % p
                aug[i] = [(aug[i][j] - f * aug[col][j]) % p for j in range(2 * n)]
    return [row[n:] for row in aug]


# ---------------------------------------------------------------- XOR 線性基


class LinearBasis:
    """GF(2) 上的線性基（63 位元）。"""

    __slots__ = ("basis",)

    def __init__(self) -> None:
        self.basis = [0] * (MAXB + 1)

    def insert(self, x: int) -> bool:
        """插入一個值；回傳是否真的擴大了基底（False 代表可被既有基底表示）。"""
        x &= MASK
        for i in range(MAXB, -1, -1):
            if (x >> i) & 1:
                if self.basis[i] == 0:
                    self.basis[i] = x
                    return True
                x ^= self.basis[i]
        return False

    def size(self) -> int:
        return sum(1 for v in self.basis if v != 0)

    def max_xor(self) -> int:
        res = 0
        for i in range(MAXB, -1, -1):
            if self.basis[i] and (res ^ self.basis[i]) > res:
                res ^= self.basis[i]
        return res

    def reduce_basis(self) -> List[int]:
        """重整基底，使 basis[i] 的低位不含其它基底的最高位；回傳由小到大的非零基底。"""
        b = list(self.basis)
        for i in range(MAXB + 1):
            if b[i] == 0:
                continue
            for j in range(i - 1, -1, -1):
                if b[j] and (b[i] >> j) & 1:
                    b[i] ^= b[j]
        return [v for v in b if v != 0]

    def kth_smallest(self, k: int) -> int:
        """第 k 小可組出的異或值（1-indexed，含 0）；k 越界回傳 -1。"""
        d = self.reduce_basis()
        total = 1 << len(d)
        if k < 1 or k > total:
            return -1
        kk = k - 1
        ans = 0
        for i, v in enumerate(d):
            if (kk >> i) & 1:
                ans ^= v
        return ans


def build_basis(xs: Sequence[int]) -> LinearBasis:
    lb = LinearBasis()
    for x in xs:
        lb.insert(x)
    return lb


# ---------------------------------------------------------------- 輸出


def run(a: List[List[int]], b: List[int], n: int, m: int, p: int,
        xs: List[int], k: int) -> str:
    rank, det, kind, free_cnt, x = solve_linear(a, b, n, m, p)
    invertible = 1 if (n == m and rank == n) else 0
    lb = build_basis(xs)
    lines = [
        f"{rank} {det}",
        f"{kind} {free_cnt}",
        " ".join(str(v) for v in x) if kind != 0 else "",
        str(invertible),
        f"{lb.size()} {lb.max_xor()}",
        str(lb.kth_smallest(k)),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 對拍用暴力


def bf_rank(a: Sequence[Sequence[int]], n: int, m: int, p: int) -> int:
    """另一種主元選取（由右往左）算秩，驗證秩與選取順序無關。"""
    mat = [[a[i][j] % p for j in range(m)] for i in range(n)]
    rank = 0
    col = m - 1
    row = 0
    while col >= 0 and row < n:
        piv = -1
        for i in range(row, n):
            if mat[i][col] % p != 0:
                piv = i
                break
        if piv < 0:
            col -= 1
            continue
        mat[row], mat[piv] = mat[piv], mat[row]
        inv = mod_inv(mat[row][col] % p, p)
        mat[row] = [(v * inv) % p for v in mat[row]]
        for i in range(n):
            if i != row and mat[i][col] % p != 0:
                f = mat[i][col] % p
                mat[i] = [(mat[i][j] - f * mat[row][j]) % p for j in range(m)]
        row += 1
        col -= 1
    return row


def bf_det(a: Sequence[Sequence[int]], n: int, p: int) -> int:
    """萊布尼茲公式（枚舉排列），僅供對拍（n ≤ 6）。"""
    total = 0
    for perm in itertools.permutations(range(n)):
        sign = 1
        # 計算排列奇偶
        seen = [False] * n
        for i in range(n):
            if not seen[i]:
                j = i
                length = 0
                while not seen[j]:
                    seen[j] = True
                    j = perm[j]
                    length += 1
                if length % 2 == 0:
                    sign = -sign
        prod = 1
        for i in range(n):
            prod = (prod * (a[i][perm[i]] % p)) % p
        total = (total + sign * prod) % p
    return total % p


def bf_solutions(a: Sequence[Sequence[int]], b: Sequence[int], n: int, m: int,
                 p: int) -> Tuple[int, Optional[List[int]]]:
    """枚舉所有 p^m 組 x，回傳 (解的個數, 字典序最小的解)；p^m 大時不要用。"""
    cnt = 0
    first: Optional[List[int]] = None
    for combo in itertools.product(range(p), repeat=m):
        ok = True
        for i in range(n):
            s = 0
            for j in range(m):
                s += (a[i][j] % p) * combo[j]
            if s % p != b[i] % p:
                ok = False
                break
        if ok:
            cnt += 1
            if first is None:
                first = list(combo)
    return cnt, first


def bf_xor_values(xs: Sequence[int]) -> List[int]:
    """枚舉所有子集異或值（含 0），排序後回傳；|xs| ≤ 16。"""
    vals = [0]
    for x in xs:
        vals = vals + [v ^ x for v in vals]
    return sorted(set(vals))


# ---------------------------------------------------------------- 測試


def check_system(a: List[List[int]], b: List[int], n: int, m: int, p: int) -> None:
    rank, det, kind, free_cnt, x = solve_linear(a, b, n, m, p)
    assert rank == bf_rank(a, n, m, p), (a, p)
    if n == m:
        assert det == bf_det(a, n, p) if n <= 6 else True, (a, p)
    if kind != 0:
        assert len(x) == m
        for i in range(n):
            s = sum((a[i][j] % p) * x[j] for j in range(m))
            assert s % p == b[i] % p, (a, b, p, x)
    if kind == 1:
        assert rank == m and free_cnt == 0
    if kind == 2:
        assert free_cnt == m - rank and free_cnt > 0
    if kind == 0:
        assert len(x) == 0

    if n == m:
        if rank == n:
            inv = matrix_inverse(a, n, p)
            assert inv is not None
            for i in range(n):
                for j in range(n):
                    s = sum((a[i][t] % p) * inv[t][j] for t in range(n)) % p
                    assert s == (1 if i == j else 0), (a, p)
        else:
            assert n == 0 or matrix_inverse(a, n, p) is None, (a, p)

    # 小到可以全枚舉時，直接數解的個數
    if p ** m <= 3000:
        cnt, first = bf_solutions(a, b, n, m, p)
        if kind == 0:
            assert cnt == 0, (a, b, p)
        elif kind == 1:
            assert cnt == 1, (a, b, p, cnt)
            assert x == first, (a, b, p, x, first)
        else:
            assert cnt == p ** (m - rank), (a, b, p, cnt, rank)
            assert first is not None


def check_basis(xs: List[int]) -> None:
    lb = build_basis(xs)
    vals = bf_xor_values(xs)
    assert len(vals) == (1 << lb.size()), (xs, vals, lb.size())
    assert lb.max_xor() == max(vals), (xs, vals)
    for k in range(1, len(vals) + 1):
        assert lb.kth_smallest(k) == vals[k - 1], (xs, k, vals)
    assert lb.kth_smallest(0) == -1
    assert lb.kth_smallest(len(vals) + 1) == -1


def run_tests() -> None:
    rng = random.Random(20261013)

    # ---- 固定用例（先跑過再寫死期望值）----
    assert to_prime(0) == 2 and to_prime(1) == 2 and to_prime(2) == 2
    assert to_prime(9) == 11 and to_prime(100) == 101
    assert to_prime(1000000007) == 1000000007

    assert run([[1, 2], [3, 4]], [5, 11], 2, 2, 7, [3, 5], 3) \
        == "2 5\n1 0\n1 2\n1\n2 6\n5"
    assert run([[1, 1], [2, 2]], [1, 2], 2, 2, 7, [7, 3], 2) \
        == "1 0\n2 1\n1 0\n0\n2 7\n3"
    assert run([[1, 1], [2, 2]], [1, 3], 2, 2, 7, [1], 1) \
        == "1 0\n0 1\n\n0\n1 1\n0"

    check_system([[1, 2], [3, 4]], [5, 11], 2, 2, 7)
    check_system([[1, 1], [2, 2]], [1, 2], 2, 2, 7)
    check_system([[1, 1], [2, 2]], [1, 3], 2, 2, 7)

    check_basis([3, 5])
    check_basis([7, 3])
    check_basis([1, 2, 3])
    check_basis([])

    # ---- 隨機對拍：小矩陣 + 全枚舉驗證 ----
    for _ in range(300):
        p = rng.choice([2, 3, 5, 7, 11])
        n = rng.randint(0, 4)
        m = rng.randint(0, 3)
        a = [[rng.randint(-6, 6) for _ in range(m)] for _ in range(n)]
        b = [rng.randint(-6, 6) for _ in range(n)]
        check_system(a, b, n, m, p)

    # ---- 隨機對拍：方陣行列式與萊布尼茲公式 ----
    for _ in range(150):
        p = rng.choice([7, 13, 101, 1000000007])
        n = rng.randint(1, 5)
        a = [[rng.randint(-9, 9) for _ in range(n)] for _ in range(n)]
        b = [rng.randint(-9, 9) for _ in range(n)]
        rank, det, _, _, _ = solve_linear(a, b, n, n, p)
        assert det == bf_det(a, n, p), (a, p)
        check_system(a, b, n, n, p)

    # ---- 隨機對拍：中矩陣（秩一致性 + 解的可行性 + 逆矩陣）----
    for _ in range(120):
        p = rng.choice([101, 1009, 1000000007])
        n = rng.randint(1, 12)
        m = rng.randint(1, 12)
        a = [[rng.randint(-20, 20) for _ in range(m)] for _ in range(n)]
        b = [rng.randint(-20, 20) for _ in range(n)]
        rank, det, kind, free_cnt, x = solve_linear(a, b, n, m, p)
        assert rank == bf_rank(a, n, m, p)
        if kind != 0:
            for i in range(n):
                s = sum((a[i][j] % p) * x[j] for j in range(m))
                assert s % p == b[i] % p
        if n == m and rank == n:
            inv = matrix_inverse(a, n, p)
            assert inv is not None
            for i in range(n):
                for j in range(n):
                    s = sum((a[i][t] % p) * inv[t][j] for t in range(n)) % p
                    assert s == (1 if i == j else 0)

    # ---- 隨機對拍：線性基 ----
    for _ in range(200):
        c = rng.randint(0, 12)
        xs = [rng.randrange(1 << rng.randint(1, 20)) for _ in range(c)]
        check_basis(xs)

    # ---- 邊界 ----
    assert run([], [], 0, 0, 7, [], 1) == "0 1\n1 0\n\n1\n0 0\n0"
    assert run([], [], 0, 3, 7, [], 1) == "0 0\n2 3\n0 0 0\n0\n0 0\n0"
    assert run([[0, 0]], [5], 1, 2, 7, [], 1) == "0 0\n0 2\n\n0\n0 0\n0"

    lb = LinearBasis()
    assert lb.size() == 0 and lb.max_xor() == 0 and lb.kth_smallest(1) == 0
    assert lb.kth_smallest(2) == -1
    lb2 = LinearBasis()
    lb2.insert(MASK)
    assert lb2.max_xor() == MASK
    assert lb2.insert(MASK) is False          # 重複插入不擴大基底

    # 63 位元極值：線性基大小上限
    lb3 = LinearBasis()
    for i in range(63):
        lb3.insert(1 << i)
    assert lb3.size() == 63
    assert lb3.max_xor() == MASK


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。"""
    tokens = raw.split()
    pos = 0

    def read_int() -> Optional[int]:
        nonlocal pos
        if pos < len(tokens) and NUM_RE.match(tokens[pos]):
            v = int(tokens[pos])
            pos += 1
            return v
        return None

    def clamp(v: Optional[int], lo: int, hi: int, default: int) -> int:
        if v is None:
            return default
        return max(lo, min(hi, v))

    n = clamp(read_int(), 0, MAXNM, 0)
    m = clamp(read_int(), 0, MAXNM, 0)
    p = to_prime(clamp(read_int(), 0, P_MAX, 2))

    a: List[List[int]] = []
    for _ in range(n):
        row: List[int] = []
        for _ in range(m):
            v = read_int()
            if v is None:
                break
            row.append(v % p)
        if len(row) < m:
            break                      # 這一列讀不滿 → 整列丟棄
        a.append(row)
    n = len(a)

    b: List[int] = []
    for _ in range(n):
        v = read_int()
        if v is None:
            break
        b.append(v % p)
    if len(b) < n:                     # 右端項不足 → 丟掉後面的列
        n = len(b)
        a = a[:n]

    c = clamp(read_int(), 0, MAXNM, 0)
    xs: List[int] = []
    for _ in range(c):
        v = read_int()
        if v is None:
            break
        xs.append(v & MASK)

    k = read_int()
    if k is None:
        k = 0

    emit(run(a, b, n, m, p, xs, k) + "\n")


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
