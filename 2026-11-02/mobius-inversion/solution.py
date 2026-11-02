"""莫比烏斯反演（Möbius Inversion）：求 gcd 恰好為 g 的數對數量

題意：
    給定 n 個正整數 a[1..n] 與 q 個查詢 g_1..g_q。
    對每個查詢 g，求有多少個**無序數對** (i, j)（1 <= i < j <= n）
    滿足 gcd(a[i], a[j]) == g。
    限制：1 <= n, q <= 100000，1 <= a[i] <= 100000，1 <= g <= 100000。

思路：
    ### 1. 從容斥講起
    設 cnt[d] = 滿足 d | a[i] 的元素個數。
    則「兩個數都是 d 的倍數」的數對數是 C(cnt[d], 2) = cnt[d](cnt[d]-1)/2。

    我們要的是 gcd **恰好等於** g，而不是「g 是它們的公因數」。
    容斥的做法是：
        答案(g) = Σ_{k>=1} μ(k) * C(cnt[g*k], 2)
    其中 μ 是莫比烏斯函數：
        μ(1) = 1
        μ(n) = 0              若 n 含平方因子
        μ(n) = (-1)^k         若 n 是 k 個相異質數的乘積
    這個式子的直觀意義是「對『g 是公因數』做帶權容斥」，
    權重就是 μ。

    ### 2. 莫比烏斯函數的定義式
    莫比烏斯反演定理：若 f(n) = Σ_{d|n} g(d)，則 g(n) = Σ_{d|n} μ(d) * f(n/d)。
    反過來，若 g(n) = Σ_{d|n} f(d)，則 f(n) = Σ_{d|n} μ(d) * g(n/d)。
    套到本題：
        f(g) = #{數對 : gcd 恰好為 g}
        g(d) = C(cnt[d], 2) = #{數對 : d 同時整除兩數}
    因為「gcd 是 d 的倍數」等價於「d 同時整除兩數」，所以 g(d) = C(cnt[d],2)，
    而 f(g) = Σ_{k>=1} μ(k) * g(g*k)，正是上面的容斥式。

    ### 3. 怎麼算 cnt[d]
    直接對每個 d 去數有幾個 a[i] 是 d 的倍數，O(n * maxA) 太慢。
    用「倒數篩」：對每個 a[i]，把它所有因數的 cnt 加一。
    因為 Σ_{d=1}^{M} M/d = O(M log M)，用一個類似埃氏篩的倍數列舉即可：
        對每個 d from 1 to M：cnt[d] = freq[d] + freq[2d] + freq[3d] + ...
    其中 freq[v] = 值 v 出現的次數。總共 O(M log M)。

    ### 4. 預先算好所有 g 的答案
    對固定 g，Σ_{k>=1} μ(k) * C(cnt[g*k], 2) 的 k 範圍是 1..M/g。
    樸素是 O(M/g)，q 個查詢總共 O(qM/g)，最壞 10^10，太慢。
    但注意：**隨著 k 增大，g*k 增大，cnt[g*k] 單調不增**，
    而 C(cnt, 2) 在 cnt <= 1 時恆為 0。
    所以一旦 cnt[g*k] <= 1，後面全是 0，可以直接砍掉。
    更重要的是：我們對**每個 g** 都預先算好答案（g 從 1 到 M），
    總共 Σ_{g=1}^{M} (M/g) = O(M log M) 次運算 —— 這已經夠快了。
    有了 ans[1..M] 之後，每個查詢 O(1) 回應。

    ### 5. 邊界
    cnt[d] 可能為 0 或 1，此時 C(cnt,2) = 0（沒有數對），
    計算時必須處理負數的乘積（C(1,2) = 1*0/2 = 0，C(0,2) = 0*(-1)/2 = 0，
    整數除法剛好也是 0，但保險起見先判 cnt >= 2）。

    莫比烏斯函數用**線性篩**求出，O(M)。

複雜度：
    線性篩求 μ            O(M)
    倍數列舉求 cnt[]     O(M log M)
    對所有 g 求 ans[]     O(M log M)
    每個查詢              O(1)
    總計                  時間 O(M log M + n + q)，空間 O(M + n)

應用場景（Applications）：
    - **數論題的批量求解**：一次反演算出所有 g 的答案，之後查詢 O(1)
    - **資訊安全**：RSA 的 φ(n) 計算、互質數對統計（金鑰生成前的檢查）
    - **生物統計**：族群中「親緣關係恰為某值」的配對數量
    - **資料清洗**：找出「共同特徵恰好為某類」的記錄對數
    - **競賽題型**：數位反演、歐拉函數求和、Divisor Sum 的各種變體
"""

import sys
from typing import List, Optional, Sequence

NMAX = 100000
QMAX = 100000
VAL_MAX = 100000


def mobius_sieve(m: int) -> List[int]:
    """線性篩求莫比烏斯函數 μ[1..m]。"""
    mu = [0] * (m + 1)
    if m >= 1:
        mu[1] = 1
    is_composite = bytearray(m + 1)
    primes: List[int] = []
    for i in range(2, m + 1):
        if not is_composite[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            x = i * p
            if x > m:
                break
            is_composite[x] = 1
            if i % p == 0:
                mu[x] = 0
                break
            mu[x] = -mu[i]
    return mu


def count_pairs_by_gcd(values: Sequence[int], max_query: int) -> List[int]:
    """計算 ans[g] = gcd 恰好為 g 的無序數對數量（g 從 0 到 max_query）。

    參數：
        values    —— 數列（正整數）
        max_query —— 需要計算的最大 g
    回傳：長度 max_query+1 的列表，ans[g] 為答案（ans[0] 恆為 0）。
    """
    n = len(values)
    if n < 2 or max_query <= 0:
        return [0] * (max_query + 1)

    m = max(max(values), max_query)
    # ---- 1. 線性篩求 μ ----
    mu = mobius_sieve(m)
    # ---- 2. 頻率表 ----
    freq = [0] * (m + 1)
    for v in values:
        freq[v] += 1
    # ---- 3. 倍數列舉求 cnt[d] = 有多少個元素是 d 的倍數 ----
    cnt = [0] * (m + 1)
    for d in range(1, m + 1):
        s = 0
        for v in range(d, m + 1, d):
            s += freq[v]
        cnt[d] = s
    # ---- 4. C(cnt[d], 2) ----
    pairs = [0] * (m + 1)
    for d in range(1, m + 1):
        c = cnt[d]
        if c >= 2:
            pairs[d] = c * (c - 1) // 2
    # ---- 5. 反演：ans[g] = Σ_k μ(k) * pairs[g*k] ----
    ans = [0] * (max_query + 1)
    for g in range(1, max_query + 1):
        total = 0
        k = 1
        limit = m // g
        while k <= limit:
            total += mu[k] * pairs[g * k]
            k += 1
        ans[g] = total
    return ans


def brute_force(values: Sequence[int], g: int) -> int:
    """參考實作：直接枚舉所有數對算 gcd（O(n^2 log V)，僅供對拍）。"""
    import math
    n = len(values)
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            if math.gcd(values[i], values[j]) == g:
                cnt += 1
    return cnt


def run_tests() -> None:
    # ---- 固定用例 ----
    # 1) 兩個數：gcd(6,10) = 2，所以 ans[2] = 1、ans[1] = 0
    ans = count_pairs_by_gcd([6, 10], 10)
    assert ans[1] == 0 and ans[2] == 1
    # 2) 兩個互質數：gcd(9,28) = 1
    assert count_pairs_by_gcd([9, 28], 10)[1] == 1
    # 3) 三個數，gcd 為 1 的有兩對：(6,10) gcd 2, (6,15) gcd 3, (10,15) gcd 5
    #    → gcd 1 = 0，gcd 2 = 1，gcd 3 = 1，gcd 5 = 1
    ans = count_pairs_by_gcd([6, 10, 15], 15)
    assert ans[1] == 0 and ans[2] == 1 and ans[3] == 1 and ans[5] == 1
    # 4) 全部相同
    assert count_pairs_by_gcd([7, 7, 7], 7)[7] == 3
    # 5) 只有一個元素 → 沒有數對
    assert count_pairs_by_gcd([5], 5)[1] == 0
    # 6) 空數列
    assert count_pairs_by_gcd([], 5)[1] == 0
    # 7) 數列 [2,4,8]：三對的 gcd 是 (2,4)=2, (2,8)=2, (4,8)=4
    ans = count_pairs_by_gcd([2, 4, 8], 8)
    assert ans[2] == 2 and ans[4] == 1 and ans[1] == 0
    # 8) 數列 [1,2,3]：三對的 gcd 分別是 1,1,1
    ans = count_pairs_by_gcd([1, 2, 3], 3)
    assert ans[1] == 3 and ans[2] == 0 and ans[3] == 0
    # 9) 數列 [12, 18, 24]：
    #    (12,18) gcd 6, (12,24) gcd 12, (18,24) gcd 6
    ans = count_pairs_by_gcd([12, 18, 24], 24)
    assert ans[6] == 2 and ans[12] == 1 and ans[1] == 0
    # 10) 數列 [4,4,4,4,4]，C(5,2) = 10 對，gcd 都是 4
    ans = count_pairs_by_gcd([4, 4, 4, 4, 4], 4)
    assert ans[4] == 10 and ans[1] == 0 and ans[2] == 0
    # 11) 大值域：gcd(100000, 99999) = 1
    ans = count_pairs_by_gcd([100000, 99999], 100000)
    assert ans[1] == 1 and ans[100000] == 0
    # 12) 相鄰整數：(8,10) 的 gcd 是 2，其餘五對都是 1
    ans = count_pairs_by_gcd([8, 9, 10, 11], 11)
    assert ans[1] == 5 and ans[2] == 1

    # ---- 與參考實作對拍 ----
    import random

    rng = random.Random(20261102)
    for _ in range(200):
        n = rng.randint(0, 12)
        values = [rng.randint(1, 24) for _ in range(n)]
        maxq = rng.randint(1, 24)
        ans = count_pairs_by_gcd(values, maxq)
        for g in range(1, maxq + 1):
            want = brute_force(values, g)
            assert ans[g] == want, (values, g, ans[g], want)

    # 值域較大
    rng = random.Random(99887766)
    for _ in range(80):
        n = rng.randint(0, 10)
        values = [rng.randint(1, 200) for _ in range(n)]
        maxq = rng.randint(1, 200)
        ans = count_pairs_by_gcd(values, maxq)
        for g in range(1, maxq + 1):
            want = brute_force(values, g)
            assert ans[g] == want, (values, g, ans[g], want)

    # 全相同的數（大量重複，cnt 很大）
    rng = random.Random(13572468)
    for _ in range(60):
        v = rng.randint(1, 30)
        n = rng.randint(2, 15)
        values = [v] * n
        maxq = max(v, rng.randint(1, 30))
        ans = count_pairs_by_gcd(values, maxq)
        want_total = n * (n - 1) // 2
        assert ans[v] == want_total, (values, v, ans[v], want_total)
        for g in range(1, maxq + 1):
            if g != v:
                assert ans[g] == 0, (values, g, ans[g])

    # 所有數都是同一個大數的倍數
    rng = random.Random(24681357)
    for _ in range(40):
        base = rng.randint(2, 12)
        n = rng.randint(2, 12)
        values = [base * rng.randint(1, 8) for _ in range(n)]
        maxq = max(values) + 5
        ans = count_pairs_by_gcd(values, maxq)
        total = 0
        for g in range(1, maxq + 1):
            want = brute_force(values, g)
            assert ans[g] == want, (values, g, ans[g], want)
            total += want
        assert total == n * (n - 1) // 2, (values, total)

    # 檢查莫比烏斯函數本身
    mu = mobius_sieve(100)
    assert mu[1] == 1
    assert mu[2] == -1 and mu[3] == -1 and mu[4] == 0 and mu[5] == -1
    assert mu[6] == 1 and mu[7] == -1 and mu[8] == 0 and mu[9] == 0
    assert mu[10] == 1 and mu[12] == 0 and mu[30] == -1
    # Σ_{d|n} μ(d) = [n == 1]
    for n in range(1, 101):
        s = 0
        d = 1
        while d * d <= n:
            if n % d == 0:
                s += mu[d]
                if d * d != n:
                    s += mu[n // d]
            d += 1
        assert s == (1 if n == 1 else 0), (n, s)


# ---------------------------------------------------------------- IO 模式

def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；非純整數或超出 ±10^18 回 None。"""
    if not tok:
        return None
    i = 0
    if tok[0] in "+-":
        i = 1
    if i == len(tok):
        return None
    for ch in tok[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(tok)
    if val > 10 ** 18 or val < -(10 ** 18):
        return None
    return val


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0
    out: List[str] = []

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    n_opt = _parse_int(nxt() or "")
    if n_opt is None or not (1 <= n_opt <= NMAX):
        bad()
        return
    q_opt = _parse_int(nxt() or "")
    if q_opt is None or not (0 <= q_opt <= QMAX):
        bad()
        return
    n = n_opt
    q = q_opt

    values: List[int] = []
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None or not (1 <= v <= VAL_MAX):
            bad()
            return
        values.append(v)

    queries: List[int] = []
    for _ in range(q):
        g = _parse_int(nxt() or "")
        if g is None or not (1 <= g <= VAL_MAX):
            bad()
            return
        queries.append(g)

    maxq = max(queries) if queries else 1
    ans = count_pairs_by_gcd(values, maxq)
    for g in queries:
        out.append(str(ans[g]))
    for line in out:
        emit(line + "\n")


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
