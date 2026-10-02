// 差分陣列與前綴和：一維區間加 / 二維子矩陣加 / 子矩陣求和 / 和為 k 的子數組個數
//
// 題意：
//     給定長度 n 的整數陣列，做 m 次區間加；再給定 R×C 矩陣，做 M 次子矩陣加；
//     最後回答一次子矩陣求和查詢，並統計「和恰為 k 的連續子數組」個數。
//     所有修改 O(1)、所有查詢 O(1)。
//
// 思路：
//     - 差分 d[i] = a[i] − a[i−1]：區間 [l, r] += x 等價於 d[l] += x、d[r+1] −= x，
//       做完所有操作後對 d 求一次前綴和即還原陣列。二維推廣為四個角的加減
//       （容斥：右下角被減兩次，要加回來）。
//     - 前綴和 S[i] = a[0]+…+a[i−1]：區間和 = S[r+1] − S[l]；
//       二維子矩陣和 = S[r2+1][c2+1] − S[r1][c2+1] − S[r2+1][c1] + S[r1][c1]。
//     - 和為 k 的子數組個數：區間和 = S[j+1] − S[i]，轉成「數有多少對前綴和相差 k」；
//       邊掃邊用哈希表記已出現的前綴和次數，初始要放 {0: 1} 代表空前綴。
//       含負數時不能用滑動窗口（窗口不具單調性）。
//     記憶點：差分用來「改」，前綴和用來「查」，兩者都把區間操作拆成端點單點操作。
//
// 輸入格式（stdin，數字以空白分隔即可）：
//     n m
//     a1 a2 … an
//     l1 r1 x1        （m 行，1-indexed 閉區間）
//     …
//     k
//     R C M
//     b(1,1) … b(1,C) （R 行，每行 C 個整數）
//     …
//     r1 c1 r2 c2 v   （M 行，1-indexed 閉區間）
//     …
//     qr1 qc1 qr2 qc2 （最後一行：子矩陣求和查詢，1-indexed 閉區間）
// 輸出格式（stdout）：
//     第 1 行：一維區間加後的陣列（n 個數；n = 0 時輸出空行）
//     第 2 行：該陣列的前綴和
//     第 3 行：和恰為 k 的連續子數組個數
//     第 4 … R+3 行：二維子矩陣加後的矩陣（每行 C 個數；R = 0 時不輸出）
//     最後一行：查詢子矩陣的元素和
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using ll = long long;

// ---------------------------------------------------------------- 輸入掃描
struct Scanner {
    vector<ll> tok;
    size_t pos = 0;
    explicit Scanner(const string& data) {
        istringstream in(data);
        ll x;
        while (in >> x) tok.push_back(x);   // 非數字 token 會讓讀取停在此處，等同截斷
    }
    ll nxt() {                              // 輸入截斷時補 0，避免中途崩掉
        if (pos < tok.size()) return tok[pos++];
        return 0;
    }
};

// ---------------------------------------------------------------- 一維差分
// 由原陣列構造長度 n+1 的差分陣列（第 n 位是哨兵，吸收 r = n−1 的減項）
static vector<ll> buildDiff1(int n, const vector<ll>& a) {
    vector<ll> d(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        d[i] += a[i];
        d[i + 1] -= a[i];
    }
    return d;
}

// 一次區間加 [l, r] += x（0-indexed 閉區間）：只動兩個端點，O(1)
static void rangeAdd1(vector<ll>& d, int n, ll l, ll r, ll x) {
    if (n <= 0) return;
    l = max<ll>(0, l);
    r = min<ll>(n - 1, r);
    if (l > r) return;
    d[(size_t)l] += x;
    d[(size_t)(r + 1)] -= x;                // r+1 <= n，哨兵位一定存在
}

// 對差分陣列求前綴和，還原出最終陣列
static vector<ll> restore1(const vector<ll>& d, int n) {
    vector<ll> res(n, 0);
    ll cur = 0;
    for (int i = 0; i < n; ++i) {
        cur += d[i];
        res[i] = cur;
    }
    return res;
}

// 暴力版區間加，只在測試裡當基準
static vector<ll> rangeAddNaive(int n, const vector<ll>& a, const vector<vector<ll>>& ops) {
    vector<ll> res = a;
    for (auto& op : ops) {
        ll l = max<ll>(0, op[0]), r = min<ll>(n - 1, op[1]), x = op[2];
        for (ll i = l; i <= r; ++i) res[(size_t)i] += x;
    }
    return res;
}

// 前綴和：out[i] = a[0] + … + a[i]
static vector<ll> prefixSum1(const vector<ll>& a) {
    vector<ll> out;
    ll s = 0;
    for (ll x : a) {
        s += x;
        out.push_back(s);
    }
    return out;
}

// 和恰為 k 的連續子數組個數：前綴和 + 哈希表，O(n)
static ll countSubarraySumK(const vector<ll>& a, ll k) {
    unordered_map<ll, ll> freq;
    freq[0] = 1;                            // 空前綴
    ll s = 0, ans = 0;
    for (ll x : a) {
        s += x;
        auto it = freq.find(s - k);
        if (it != freq.end()) ans += it->second;
        freq[s] += 1;
    }
    return ans;
}

// 暴力枚舉所有子數組，O(n^2)，只在測試裡當基準
static ll countSubarraySumKNaive(const vector<ll>& a, ll k) {
    int n = (int)a.size();
    ll ans = 0;
    for (int i = 0; i < n; ++i) {
        ll s = 0;
        for (int j = i; j < n; ++j) {
            s += a[j];
            if (s == k) ans++;
        }
    }
    return ans;
}

// ---------------------------------------------------------------- 二維差分
using Mat = vector<vector<ll>>;

// 由原矩陣構造 (R+1)×(C+1) 的二維差分陣列（多出來的一行一列是哨兵）
static Mat buildDiff2(int R, int C, const Mat& b) {
    Mat d(R + 1, vector<ll>(C + 1, 0));
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < C; ++j) {
            ll v = b[i][j];
            d[i][j] += v;
            d[i + 1][j] -= v;
            d[i][j + 1] -= v;
            d[i + 1][j + 1] += v;
        }
    return d;
}

// 一次子矩陣加 [r1..r2]×[c1..c2] += v：動四個角，O(1)
static void submatrixAdd2(Mat& d, int R, int C, ll r1, ll c1, ll r2, ll c2, ll v) {
    if (R <= 0 || C <= 0) return;
    r1 = max<ll>(0, r1);
    c1 = max<ll>(0, c1);
    r2 = min<ll>(R - 1, r2);
    c2 = min<ll>(C - 1, c2);
    if (r1 > r2 || c1 > c2) return;
    d[r1][c1] += v;
    d[r2 + 1][c1] -= v;
    d[r1][c2 + 1] -= v;
    d[r2 + 1][c2 + 1] += v;                 // 被減兩次，加回來
}

// 二維前綴和還原：g[i][j] = d[i][j] + 上 + 左 − 左上
static Mat restore2(const Mat& d, int R, int C) {
    Mat g(R, vector<ll>(C, 0));
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < C; ++j) {
            ll up = i ? g[i - 1][j] : 0;
            ll left = j ? g[i][j - 1] : 0;
            ll diag = (i && j) ? g[i - 1][j - 1] : 0;
            g[i][j] = d[i][j] + up + left - diag;
        }
    return g;
}

// 二維前綴和表 S，S[i][j] = 左上角 (0,0) 到 (i−1,j−1) 的矩形和
static Mat prefix2(const Mat& g, int R, int C) {
    Mat S(R + 1, vector<ll>(C + 1, 0));
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < C; ++j)
            S[i + 1][j + 1] = g[i][j] + S[i][j + 1] + S[i + 1][j] - S[i][j];
    return S;
}

// 子矩陣求和（容斥），O(1)
static ll submatrixSum(const Mat& S, int R, int C, ll r1, ll c1, ll r2, ll c2) {
    if (R <= 0 || C <= 0) return 0;
    r1 = max<ll>(0, r1);
    c1 = max<ll>(0, c1);
    r2 = min<ll>(R - 1, r2);
    c2 = min<ll>(C - 1, c2);
    if (r1 > r2 || c1 > c2) return 0;
    return S[r2 + 1][c2 + 1] - S[r1][c2 + 1] - S[r2 + 1][c1] + S[r1][c1];
}

// 暴力版子矩陣加，只在測試裡當基準
static Mat submatrixAddNaive(int R, int C, const Mat& b, const vector<vector<ll>>& ops) {
    Mat g = b;
    for (auto& op : ops) {
        ll r1 = max<ll>(0, op[0]), c1 = max<ll>(0, op[1]);
        ll r2 = min<ll>(R - 1, op[2]), c2 = min<ll>(C - 1, op[3]), v = op[4];
        for (ll i = r1; i <= r2; ++i)
            for (ll j = c1; j <= c2; ++j) g[(size_t)i][(size_t)j] += v;
    }
    return g;
}

// 暴力版子矩陣求和，只在測試裡當基準
static ll submatrixSumNaive(const Mat& g, int R, int C, ll r1, ll c1, ll r2, ll c2) {
    if (R <= 0 || C <= 0) return 0;
    r1 = max<ll>(0, r1);
    c1 = max<ll>(0, c1);
    r2 = min<ll>(R - 1, r2);
    c2 = min<ll>(C - 1, c2);
    if (r1 > r2 || c1 > c2) return 0;
    ll s = 0;
    for (ll i = r1; i <= r2; ++i)
        for (ll j = c1; j <= c2; ++j) s += g[(size_t)i][(size_t)j];
    return s;
}

// ---------------------------------------------------------------- 小工具
static void printRow(const vector<ll>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) cout << ' ';
        cout << v[i];
    }
    cout << '\n';
}

static mt19937 rngEngine;

static int rndInt(int lo, int hi) {          // [lo, hi]
    return lo + (int)(rngEngine() % (unsigned)(hi - lo + 1));
}

// ---------------------------------------------------------------- IO 模式
static void runIo(const string& data) {
    Scanner sc(data);
    int n = (int)sc.nxt();
    int m = (int)sc.nxt();
    vector<ll> a;
    for (int i = 0; i < max(0, n); ++i) a.push_back(sc.nxt());
    vector<vector<ll>> ops;
    for (int i = 0; i < max(0, m); ++i) {
        ll l = sc.nxt() - 1;                // 輸入是 1-indexed，轉 0-indexed
        ll r = sc.nxt() - 1;
        ll x = sc.nxt();
        ops.push_back({l, r, x});
    }

    n = max(0, n);
    vector<ll> d1 = buildDiff1(n, a);
    for (auto& op : ops) rangeAdd1(d1, n, op[0], op[1], op[2]);
    vector<ll> arr = restore1(d1, n);
    printRow(arr);
    printRow(prefixSum1(arr));

    ll k = sc.nxt();
    cout << countSubarraySumK(arr, k) << '\n';

    int R = (int)sc.nxt();
    int C = (int)sc.nxt();
    int M = (int)sc.nxt();
    Mat b(max(0, R), vector<ll>(max(0, C), 0));
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < C; ++j) b[i][j] = sc.nxt();
    vector<vector<ll>> ops2;
    for (int i = 0; i < max(0, M); ++i) {
        ll r1 = sc.nxt() - 1, c1 = sc.nxt() - 1;
        ll r2 = sc.nxt() - 1, c2 = sc.nxt() - 1;
        ll v = sc.nxt();
        ops2.push_back({r1, c1, r2, c2, v});
    }
    ll qr1 = sc.nxt() - 1, qc1 = sc.nxt() - 1;
    ll qr2 = sc.nxt() - 1, qc2 = sc.nxt() - 1;

    if (R > 0 && C > 0) {
        Mat d2 = buildDiff2(R, C, b);
        for (auto& op : ops2) submatrixAdd2(d2, R, C, op[0], op[1], op[2], op[3], op[4]);
        Mat g = restore2(d2, R, C);
        for (auto& row : g) printRow(row);
        cout << submatrixSum(prefix2(g, R, C), R, C, qr1, qc1, qr2, qc2) << '\n';
    } else {
        cout << 0 << '\n';
    }
}

// ---------------------------------------------------------------- 測試
static void runTests() {
    // 固定用例：README 的示例
    vector<ll> a = {1, 2, 3, 4, 5};
    vector<vector<ll>> ops = {{1, 3, 10}, {3, 4, 1}};     // [1,3] += 10，[3,4] += 1
    vector<ll> d = buildDiff1(5, a);
    for (auto& op : ops) rangeAdd1(d, 5, op[0], op[1], op[2]);
    vector<ll> got = restore1(d, 5);
    assert(got == rangeAddNaive(5, a, ops));
    assert((got == vector<ll>{1, 12, 13, 15, 6}));
    assert((prefixSum1(got) == vector<ll>{1, 13, 26, 41, 47}));
    assert(countSubarraySumK(got, 13) == 2);              // [1, 12] 與 [13]
    assert(countSubarraySumK(got, 13) == countSubarraySumKNaive(got, 13));

    // 空陣列
    assert(restore1(buildDiff1(0, {}), 0).empty());
    assert(prefixSum1({}).empty());
    assert(countSubarraySumK({}, 0) == 0);
    assert(countSubarraySumK({}, 7) == 0);

    // 前綴和為 0 的經典陷阱：全 0 陣列 + k = 0 → 每個子數組都算
    assert(countSubarraySumK({0, 0, 0}, 0) == 6);         // 3*4/2
    assert(countSubarraySumK({3, 4, -7, 1}, 0) == 1);     // 只有 [3, 4, -7]
    assert(countSubarraySumK({3, 4, -7, 1}, 0) == countSubarraySumKNaive({3, 4, -7, 1}, 0));

    // 二維固定用例：3×3 全 0，兩次子矩陣加
    Mat b(3, vector<ll>(3, 0));
    vector<vector<ll>> ops2 = {{0, 0, 1, 1, 5}, {1, 1, 2, 2, 3}};
    Mat d2 = buildDiff2(3, 3, b);
    for (auto& op : ops2) submatrixAdd2(d2, 3, 3, op[0], op[1], op[2], op[3], op[4]);
    Mat g = restore2(d2, 3, 3);
    assert(g == submatrixAddNaive(3, 3, b, ops2));
    assert((g == Mat{{5, 5, 0}, {5, 8, 3}, {0, 3, 3}}));
    Mat S = prefix2(g, 3, 3);
    assert(submatrixSum(S, 3, 3, 0, 0, 2, 2) == 32);
    assert(submatrixSum(S, 3, 3, 0, 0, 2, 2) == submatrixSumNaive(g, 3, 3, 0, 0, 2, 2));
    assert(submatrixSum(S, 3, 3, 1, 1, 1, 1) == 8);
    assert(submatrixSum(S, 3, 3, 0, 0, 0, 2) == 10);

    // 隨機對拍
    rngEngine.seed(20261002);
    for (int t = 0; t < 500; ++t) {
        int n = rndInt(0, 8), m = rndInt(0, 5);
        vector<ll> a;
        for (int i = 0; i < n; ++i) a.push_back(rndInt(-9, 9));
        vector<vector<ll>> ops;
        for (int i = 0; i < m; ++i) {
            if (n == 0) {
                ops.push_back({0, -1, (ll)rndInt(-5, 5)});       // 空區間，應被忽略
            } else {
                int l = rndInt(0, n - 1), r = rndInt(l, n - 1);
                ops.push_back({l, r, (ll)rndInt(-9, 9)});
            }
        }
        if (n > 0 && rndInt(0, 2) == 0) ops.push_back({-3, n + 3, 7});   // 越界區間

        vector<ll> d = buildDiff1(n, a);
        for (auto& op : ops) rangeAdd1(d, n, op[0], op[1], op[2]);
        vector<ll> arr = restore1(d, n);
        assert(arr == rangeAddNaive(n, a, ops));

        ll k = rndInt(-12, 12);
        assert(countSubarraySumK(arr, k) == countSubarraySumKNaive(arr, k));

        int R = rndInt(0, 4), C = rndInt(0, 4);
        Mat b(R, vector<ll>(C, 0));
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < C; ++j) b[i][j] = rndInt(-9, 9);
        vector<vector<ll>> ops2;
        int M = rndInt(0, 4);
        for (int i = 0; i < M; ++i) {
            if (R == 0 || C == 0) {
                ops2.push_back({0, 0, -1, -1, (ll)rndInt(-5, 5)});
            } else {
                int r1 = rndInt(0, R - 1), r2 = rndInt(r1, R - 1);
                int c1 = rndInt(0, C - 1), c2 = rndInt(c1, C - 1);
                ops2.push_back({r1, c1, r2, c2, (ll)rndInt(-9, 9)});
            }
        }
        Mat d2 = buildDiff2(R, C, b);
        for (auto& op : ops2) submatrixAdd2(d2, R, C, op[0], op[1], op[2], op[3], op[4]);
        Mat g = restore2(d2, R, C);
        if (R > 0 && C > 0) {
            assert(g == submatrixAddNaive(R, C, b, ops2));
            Mat S = prefix2(g, R, C);
            for (int q = 0; q < 3; ++q) {
                int r1 = rndInt(0, R - 1), r2 = rndInt(r1, R - 1);
                int c1 = rndInt(0, C - 1), c2 = rndInt(c1, C - 1);
                assert(submatrixSum(S, R, C, r1, c1, r2, c2) ==
                       submatrixSumNaive(g, R, C, r1, c1, r2, c2));
            }
        }
    }

    cout << "all tests passed" << '\n';
}

int main() {
    string data, line;
    bool hasInput = false;
    while (getline(cin, line)) {
        data += line;
        data += '\n';
        if (!line.empty()) hasInput = true;
    }
    if (hasInput) runIo(data);
    else runTests();
    return 0;
}
