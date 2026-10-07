// 狀態壓縮 DP（TSP 旅行商 / 最短哈密頓路徑 / 集合劃分 / SOS 子集和 DP）
//
// 題意：
//     給定 n 個點（n ≤ 14）與一張**有向帶權圖**的距離矩陣 `dist`，以及一個長度 2^n 的成本陣列 `cost`，
//     用狀態壓縮動態規劃解決四個問題：
//       1. **TSP 最短迴路**：從 0 號點出發，訪問每個點恰好一次，最後回到 0 號點的最小總代價與具體路徑；
//       2. **最短哈密頓路徑**：從 0 號點出發訪問全部點，但**不要求**回到 0 號點；
//       3. **集合劃分**：把 {0..n-1} 劃分成若干非空子集，總成本 = Σ cost[子集]，求最小總成本與具體劃分；
//       4. **SOS DP（Sum over Subsets）**：對每個 mask 求 f[mask] = Σ_{sub ⊆ mask} cost[sub]。
//
// 思路：
//     ### 狀態壓縮的核心
//     當 n 不大（≤ 20 左右）時，把「一個子集」編碼成一個整數 mask（第 i 位為 1 表示 i 在集合內），
//     於是以「子集」為狀態的 DP 就變成一維陣列下標，位運算代替集合運算：
//       - 加入元素 i：`mask | (1 << i)`
//       - 移除元素 i：`mask & ~(1 << i)`
//       - 枚舉 mask 的所有子集：`for (sub = mask; ; sub = (sub - 1) & mask)`（含空集，O(2^popcount)）
//
//     ### TSP：dp[mask][u]
//     `dp[mask][u]` = 從 0 出發、已訪問集合恰為 mask、當前停在 u 的最小代價。
//     轉移：枚舉下一個未訪問的點 v，
//         dp[mask | 1<<v][v] = min(dp[mask | 1<<v][v], dp[mask][u] + dist[u][v])
//     共 2^n · n 個狀態、每狀態 n 次轉移 → **O(2^n · n^2)** 時間、**O(2^n · n)** 空間。
//     相比枚舉所有排列的 O(n!)，n = 14 時從 8.7e10 降到約 3.2e6，這就是狀態壓縮的威力。
//     答案 = `min_v dp[ALL][v] + dist[v][0]`（補上回到起點那一跳）。
//     路徑用 `parent[mask][v]` 記錄前驅，回溯即可得到。
//
//     ### 最短哈密頓路徑
//     同一張 dp 表，只是答案不補最後一跳：`min_v dp[ALL][v]`。
//
//     ### 集合劃分：dp[mask]
//     `dp[mask]` = 把 mask 這個集合劃分成若干非空子集的最小總成本。
//     枚舉包含 mask 最低位的那個子集 s（這樣每個劃分只被算一次，避免重複計數）：
//         dp[mask] = min_{s ⊆ mask, s 含 lowestbit(mask)} ( dp[mask ^ s] + cost[s] )
//     枚舉量是 3^n 的子集和級別 → **O(3^n)**，比枚舉所有劃分（Bell 數，B(14) ≈ 1.9e8）小得多。
//
//     ### SOS DP
//     `f[mask] = Σ_{sub ⊆ mask} a[sub]` 若暴力枚舉子集是 O(3^n)；
//     SOS DP 按位處理：對每一位 i，把「去掉第 i 位」的結果累加進來：
//         for i in 0..n-1: for mask: if mask 的第 i 位為 1: f[mask] += f[mask ^ (1<<i)]
//     → **O(n · 2^n)**。常用於「帶位掩碼的數位 DP / 子集計數 / 相容性統計」。
//
//     ### 平手規則（保證 Python 與 C++ 輸出逐字節一致）
//       - 狀態枚舉順序、v 的枚舉順序皆為升序；
//       - 轉移只在**嚴格更優**時更新（先來先佔）；
//       - 集合劃分用「降序枚舉子集 + 嚴格更優才更新」。
//
// 應用場景：
//     物流配送路線規劃、晶片鑽孔 / 焊接路徑優化（TSP 的直接應用）、
//     任務排程與分批（集合劃分）、競賽中的位掩碼計數（SOS DP）、基因組組裝的 overlap 圖。
//
// 複雜度：
//     TSP / 哈密頓路徑   O(2^n · n^2) 時間、O(2^n · n) 空間
//     集合劃分           O(3^n) 時間、O(2^n) 空間
//     SOS DP             O(n · 2^n) 時間、O(2^n) 空間
//
// 輸入格式（stdin，全部以空白分隔）：
//     n
//     n × n 個整數：dist[i][j]（i 行 j 列；可為負，但不要求對稱）
//     2^n 個整數：cost[0] .. cost[2^n - 1]
//     讀不到那麼多時，缺的部分補 0；遇到非整數 token 視為輸入結束。
//     為避免狀態數爆炸，n 會被截斷到 14。
// 輸出格式（stdout）：
//     第 1 行：TSP 最短迴路代價
//     第 2 行：TSP 路徑（0-indexed，空格分隔，首尾皆為 0）
//     第 3 行：最短哈密頓路徑代價
//     第 4 行：哈密頓路徑（起點為 0）
//     第 5 行：集合劃分最小成本
//     第 6 行：劃分組數
//     第 7 行起：每組一行（組內元素升序，組間按組內最小元素升序），共「組數」行
//     最後一行：SOS DP 的 f[mask]（0 ≤ mask < 2^n，空格分隔）
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static const long long INF = 1000000000000000000LL;   // 與 Python 版的 INF 完全一致
static const int MAX_N = 14;                          // 狀態數上限

// 嚴格整數規則：只接受 [+-]?digits，與 Python 的 re 版完全一致
static bool tryLL(const string &tok, long long &out) {
    if (tok.empty()) return false;
    size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i >= tok.size()) return false;
    for (size_t k = i; k < tok.size(); ++k) {
        if (!isdigit((unsigned char)tok[k])) return false;
    }
    bool neg = (tok[0] == '-');
    string digits = tok.substr(i);
    while (digits.size() > 1 && digits[0] == '0') digits.erase(digits.begin());
    if (digits.size() > 18) digits = digits.substr(digits.size() - 18);   // 保底截斷，避免 UB
    long long v = 0;
    for (char c : digits) v = v * 10 + (c - '0');
    out = neg ? -v : v;
    return true;
}

static int popcount_int(int x) { return __builtin_popcount((unsigned)x); }

static vector<int> bits_of(int mask) {
    vector<int> out;
    int i = 0;
    while (mask) {
        if (mask & 1) out.push_back(i);
        mask >>= 1;
        ++i;
    }
    return out;
}

// ---------------------------------------------------------------- TSP

struct TspResult {
    long long cost;
    vector<int> path;
};

static TspResult tsp(const vector<vector<long long>> &dist, int n) {
    TspResult res;
    res.cost = 0;
    if (n == 0) return res;
    int size = 1 << n;
    vector<long long> dp((size_t)size * n, INF);
    vector<int> par((size_t)size * n, -1);
    dp[(size_t)1 * n + 0] = 0;
    for (int mask = 0; mask < size; ++mask) {
        if (!(mask & 1)) continue;
        size_t base = (size_t)mask * n;
        for (int u = 0; u < n; ++u) {
            long long cur = dp[base + u];
            if (cur == INF) continue;
            for (int v = 0; v < n; ++v) {
                if (mask >> v & 1) continue;
                int nmask = mask | (1 << v);
                long long cand = cur + dist[u][v];
                size_t idx = (size_t)nmask * n + v;
                if (cand < dp[idx]) {
                    dp[idx] = cand;
                    par[idx] = u;
                }
            }
        }
    }
    int full = size - 1;
    int best_v = -1;
    long long best_c = INF;
    for (int v = 0; v < n; ++v) {
        long long cur = dp[(size_t)full * n + v];
        if (cur == INF) continue;
        long long c = cur + dist[v][0];
        if (c < best_c) { best_c = c; best_v = v; }
    }
    if (best_v < 0) return res;
    int mask = full, u = best_v;
    while (true) {
        res.path.push_back(u);
        int p = par[(size_t)mask * n + u];
        if (p == -1) break;
        mask ^= (1 << u);
        u = p;
    }
    reverse(res.path.begin(), res.path.end());
    res.path.push_back(0);                       // 補上回到起點
    res.cost = best_c;
    return res;
}

static long long tsp_bruteforce(const vector<vector<long long>> &dist, int n) {
    if (n == 0) return 0;
    if (n == 1) return dist[0][0];
    vector<int> perm;
    for (int i = 1; i < n; ++i) perm.push_back(i);
    long long best = INF;
    do {
        long long c = dist[0][perm[0]];
        for (int i = 0; i + 1 < (int)perm.size(); ++i) c += dist[perm[i]][perm[i + 1]];
        c += dist[perm.back()][0];
        if (c < best) best = c;
    } while (next_permutation(perm.begin(), perm.end()));
    return best;
}

static TspResult hamiltonian_path(const vector<vector<long long>> &dist, int n) {
    TspResult res;
    res.cost = 0;
    if (n == 0) return res;
    int size = 1 << n;
    vector<long long> dp((size_t)size * n, INF);
    vector<int> par((size_t)size * n, -1);
    dp[(size_t)1 * n + 0] = 0;
    for (int mask = 0; mask < size; ++mask) {
        if (!(mask & 1)) continue;
        size_t base = (size_t)mask * n;
        for (int u = 0; u < n; ++u) {
            long long cur = dp[base + u];
            if (cur == INF) continue;
            for (int v = 0; v < n; ++v) {
                if (mask >> v & 1) continue;
                int nmask = mask | (1 << v);
                long long cand = cur + dist[u][v];
                size_t idx = (size_t)nmask * n + v;
                if (cand < dp[idx]) {
                    dp[idx] = cand;
                    par[idx] = u;
                }
            }
        }
    }
    int full = size - 1;
    int best_v = -1;
    long long best_c = INF;
    for (int v = 0; v < n; ++v) {
        long long cur = dp[(size_t)full * n + v];
        if (cur == INF) continue;
        if (cur < best_c) { best_c = cur; best_v = v; }
    }
    if (best_v < 0) return res;
    int mask = full, u = best_v;
    while (true) {
        res.path.push_back(u);
        int p = par[(size_t)mask * n + u];
        if (p == -1) break;
        mask ^= (1 << u);
        u = p;
    }
    reverse(res.path.begin(), res.path.end());
    res.cost = best_c;
    return res;
}

static long long hamiltonian_bruteforce(const vector<vector<long long>> &dist, int n) {
    if (n == 0) return 0;
    if (n == 1) return 0;
    vector<int> perm;
    for (int i = 1; i < n; ++i) perm.push_back(i);
    long long best = INF;
    do {
        long long c = dist[0][perm[0]];
        for (int i = 0; i + 1 < (int)perm.size(); ++i) c += dist[perm[i]][perm[i + 1]];
        if (c < best) best = c;
    } while (next_permutation(perm.begin(), perm.end()));
    return best;
}

// ---------------------------------------------------------------- 集合劃分

static pair<long long, vector<vector<int>>> min_partition_cost(const vector<long long> &cost, int n) {
    int size = 1 << n;
    vector<long long> dp(size, INF);
    vector<int> choice(size, 0);
    dp[0] = 0;
    for (int mask = 1; mask < size; ++mask) {
        int low = mask & -mask;
        int rest = mask ^ low;
        int sub = rest;
        while (true) {
            int s = sub | low;
            long long prev = dp[mask ^ s];
            if (prev != INF) {
                long long cand = prev + cost[s];
                if (cand < dp[mask]) {           // 嚴格更優才更新 → 確定性
                    dp[mask] = cand;
                    choice[mask] = s;
                }
            }
            if (sub == 0) break;
            sub = (sub - 1) & rest;
        }
    }
    vector<vector<int>> groups;
    int mask = size - 1;
    while (mask) {
        int s = choice[mask];
        groups.push_back(bits_of(s));
        mask ^= s;
    }
    sort(groups.begin(), groups.end(), [](const vector<int> &a, const vector<int> &b) {
        return a[0] < b[0];
    });
    return make_pair(dp[size - 1], groups);
}

static void partition_rec(int i, int n, const vector<long long> &cost,
                          vector<int> &groups, long long &best) {
    if (i == n) {
        long long total = 0;
        for (int g : groups) total += cost[g];
        if (total < best) best = total;
        return;
    }
    for (int k = 0; k < (int)groups.size(); ++k) {
        groups[k] |= (1 << i);
        partition_rec(i + 1, n, cost, groups, best);
        groups[k] ^= (1 << i);
    }
    groups.push_back(1 << i);
    partition_rec(i + 1, n, cost, groups, best);
    groups.pop_back();
}

static long long partition_bruteforce(const vector<long long> &cost, int n) {
    if (n == 0) return 0;
    long long best = INF;
    vector<int> groups;
    partition_rec(0, n, cost, groups, best);
    return best;
}

// ---------------------------------------------------------------- SOS DP

static vector<long long> sos_dp(const vector<long long> &a, int n) {
    int size = 1 << n;
    vector<long long> f(size, 0);
    for (int i = 0; i < size && i < (int)a.size(); ++i) f[i] = a[i];
    for (int i = 0; i < n; ++i) {
        int bit = 1 << i;
        for (int mask = 0; mask < size; ++mask) {
            if (mask & bit) f[mask] += f[mask ^ bit];
        }
    }
    return f;
}

static vector<long long> sos_bruteforce(const vector<long long> &a, int n) {
    int size = 1 << n;
    vector<long long> out(size, 0);
    for (int mask = 0; mask < size; ++mask) {
        long long total = 0;
        int sub = mask;
        while (true) {
            total += (sub < (int)a.size() ? a[sub] : 0);
            if (sub == 0) break;
            sub = (sub - 1) & mask;
        }
        out[mask] = total;
    }
    return out;
}

// ---------------------------------------------------------------- IO 與測試

static string join_ints(const vector<int> &v) {
    string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += " ";
        out += to_string(v[i]);
    }
    return out;
}

static string join_ll(const vector<long long> &v) {
    string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += " ";
        out += to_string(v[i]);
    }
    return out;
}

static void run_io(const string &raw) {
    vector<string> toks;
    {
        istringstream iss(raw);
        string t;
        while (iss >> t) toks.push_back(t);
    }
    size_t pos = 0;
    auto nxt = [&]() -> long long {
        long long v = 0;
        if (pos < toks.size()) {
            long long parsed = 0;
            if (tryLL(toks[pos], parsed)) v = parsed;
        }
        ++pos;
        return v;
    };

    long long nll = nxt();
    int n = (int)nll;
    if (n < 0) n = 0;
    if (n > MAX_N) n = MAX_N;
    vector<vector<long long>> dist(n, vector<long long>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) dist[i][j] = nxt();
    int size = 1 << n;
    vector<long long> cost(size, 0);
    for (int i = 0; i < size; ++i) cost[i] = nxt();

    TspResult tr = tsp(dist, n);
    TspResult hr = hamiltonian_path(dist, n);
    pair<long long, vector<vector<int>>> pr = min_partition_cost(cost, n);
    vector<long long> f = sos_dp(cost, n);

    vector<string> out;
    out.push_back(to_string(tr.cost));
    out.push_back(join_ints(tr.path));
    out.push_back(to_string(hr.cost));
    out.push_back(join_ints(hr.path));
    out.push_back(to_string(pr.first));
    out.push_back(to_string((int)pr.second.size()));
    for (const vector<int> &g : pr.second) out.push_back(join_ints(g));
    out.push_back(join_ll(f));
    for (size_t i = 0; i < out.size(); ++i) cout << out[i] << "\n";
}

static void run_tests() {
    assert(popcount_int(0) == 0);
    assert(popcount_int(0b101101) == 4);
    assert(bits_of(0b10110) == vector<int>({1, 2, 4}));

    // ---- SOS DP 對拍 ----
    vector<long long> a{1, 2, 3, 4};
    assert(sos_dp(a, 2) == vector<long long>({1, 3, 4, 10}));
    assert(sos_dp(a, 2) == sos_bruteforce(a, 2));

    // ---- 固定用例：對稱 TSP ----
    vector<vector<long long>> d4 = {
        {0, 10, 15, 20},
        {10, 0, 35, 25},
        {15, 35, 0, 30},
        {20, 25, 30, 0},
    };
    assert(tsp(d4, 4).cost == 80);                  // 0→1→3→2→0 = 10+25+30+15
    assert(tsp(d4, 4).cost == tsp_bruteforce(d4, 4));
    assert(hamiltonian_path(d4, 4).cost == 65);     // 0→1→3→2 = 10+25+30
    assert(hamiltonian_path(d4, 4).cost == hamiltonian_bruteforce(d4, 4));

    vector<int> tp = tsp(d4, 4).path;
    assert(tp.front() == 0 && tp.back() == 0);
    {
        vector<int> mid(tp.begin() + 1, tp.end() - 1);
        sort(mid.begin(), mid.end());
        assert(mid == vector<int>({1, 2, 3}));
    }
    vector<int> hp = hamiltonian_path(d4, 4).path;
    assert(hp.front() == 0);
    {
        vector<int> s = hp;
        sort(s.begin(), s.end());
        assert(s == vector<int>({0, 1, 2, 3}));
    }

    // ---- 邊界 ----
    assert(tsp(vector<vector<long long>>(), 0).cost == 0);
    assert(hamiltonian_path(vector<vector<long long>>(), 0).cost == 0);
    assert(min_partition_cost(vector<long long>({0}), 0).first == 0);
    assert(sos_dp(vector<long long>({7}), 0) == vector<long long>({7}));
    {
        vector<vector<long long>> d1{{5}};
        assert(tsp(d1, 1).cost == 5);               // 只有一個點：0 → 0
        assert(tsp(d1, 1).path == vector<int>({0, 0}));
        assert(hamiltonian_path(d1, 1).cost == 0);
        assert(hamiltonian_path(d1, 1).path == vector<int>({0}));
    }

    mt19937 rng(20261007);

    // ---- 隨機對拍：TSP / 哈密頓路徑 vs 排列枚舉（n ≤ 7）----
    for (int t = 0; t < 120; ++t) {
        int n = 1 + (int)(rng() % 7);
        vector<vector<long long>> dist(n, vector<long long>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) dist[i][j] = (long long)(rng() % 31);
        TspResult tr = tsp(dist, n);
        assert(tr.cost == tsp_bruteforce(dist, n));
        assert(tr.path.front() == 0 && tr.path.back() == 0);
        {
            vector<int> mid(tr.path.begin() + 1, tr.path.end() - 1);
            sort(mid.begin(), mid.end());
            vector<int> want;
            for (int i = 1; i < n; ++i) want.push_back(i);
            assert(mid == want);
        }
        TspResult hr = hamiltonian_path(dist, n);
        assert(hr.cost == hamiltonian_bruteforce(dist, n));
        assert(hr.path.front() == 0);
        {
            vector<int> s = hr.path;
            sort(s.begin(), s.end());
            vector<int> want;
            for (int i = 0; i < n; ++i) want.push_back(i);
            assert(s == want);
        }
    }

    // ---- 隨機對拍：含負權邊 ----
    for (int t = 0; t < 120; ++t) {
        int n = 1 + (int)(rng() % 7);
        vector<vector<long long>> dist(n, vector<long long>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) dist[i][j] = (long long)(rng() % 41) - 20;
        assert(tsp(dist, n).cost == tsp_bruteforce(dist, n));
        assert(hamiltonian_path(dist, n).cost == hamiltonian_bruteforce(dist, n));
    }

    // ---- 隨機對拍：集合劃分 vs 枚舉所有劃分（n ≤ 6）----
    for (int t = 0; t < 120; ++t) {
        int n = 1 + (int)(rng() % 6);
        int size = 1 << n;
        vector<long long> cost(size);
        for (int i = 0; i < size; ++i) cost[i] = (long long)(rng() % 41) - 10;
        pair<long long, vector<vector<int>>> pr = min_partition_cost(cost, n);
        assert(pr.first == partition_bruteforce(cost, n));
        vector<int> covered;
        for (const vector<int> &g : pr.second)
            for (int x : g) covered.push_back(x);
        sort(covered.begin(), covered.end());
        {
            vector<int> want;
            for (int i = 0; i < n; ++i) want.push_back(i);
            assert(covered == want);              // 劃分必須恰好覆蓋每個元素一次
        }
        long long recomputed = 0;
        for (const vector<int> &g : pr.second) {
            int m = 0;
            for (int x : g) m |= (1 << x);
            recomputed += cost[m];
        }
        assert(recomputed == pr.first);
    }

    // ---- 隨機對拍：SOS DP vs 暴力（n ≤ 8）----
    for (int t = 0; t < 120; ++t) {
        int n = (int)(rng() % 9);
        int size = 1 << n;
        vector<long long> arr(size);
        for (int i = 0; i < size; ++i) arr[i] = (long long)(rng() % 21) - 10;
        assert(sos_dp(arr, n) == sos_bruteforce(arr, n));
    }

    cout << "all tests passed" << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ostringstream oss;
    oss << cin.rdbuf();
    string raw = oss.str();

    bool has_input = raw.find_first_not_of(" \t\r\n") != string::npos;
    if (has_input) run_io(raw);
    else run_tests();
    return 0;
}
