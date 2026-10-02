// 背包變體：多重背包（二進制拆分）/ 完全背包 / 湊齊容量的最少件數
//
// 題意：
//     給定 n 種物品（重量 w、價值 v、數量上限 c）與容量 C，依次回答三問：
//       1. 多重背包：每種最多取 c_i 個，容量 C 內的最大價值 + 一組達到最優的選取個數；
//       2. 完全背包：每種無限取，容量 C 內的最大價值 + 一組達到最優的選取個數；
//       3. 每種無限取，湊出恰好總重 C 的最少件數（無解輸出 −1）。
//
// 思路：
//     0-1 背包的 1D 寫法是容量「倒序」遍歷，讓 dp[cap − w] 讀到上一輪（還沒放過本物品）的值。
//     把方向反過來就得到完全背包：正序時 dp[cap − w] 是本輪已更新過的值，
//     相當於本物品可被反覆加入，一次掃描等價於枚舉取 0,1,2,… 個。O(n·C)。
//     多重背包既不能倒序（會退化成 0-1）也不能正序（會變成無限），標準做法是二進制拆分：
//     把 c 拆成 1, 2, 4, …, 剩餘 的若干「捆」，每捆當成一個獨立 0-1 物品（k·w, k·v）。
//     1..c 的任意整數都能由這組 2 的冪唯一表示，因此捆的 0-1 組合恰好覆蓋「取 0..c 個」，
//     複雜度由 O(n·C·c) 降到 O(n·C·log c)。
//     最少件數則是把 max 換成 min、初值換 INF，dp[0] = 0，容量正序。
//
//     還原方案：
//       - 多重背包用二維 DP（dp[i][cap]），從 i = n 往回，枚舉每種取了幾個 t，
//         找第一個滿足 dp[i][cap] == dp[i−1][cap − t·w] + t·v 的 t（取最小的 t，輸出唯一）。
//       - 完全背包用 1D DP + par[cap] 記錄「該容量最後被哪種物品更新」，從 cap = C 往回減。
//         par 回溯對 0-1 / 多重不成立（會重複用同一捆而超出數量上限），只對「可重複使用」成立。
//
//     約定：重量必須 ≥ 1；w = 0 的物品直接忽略（否則價值無界）。
//
// 輸入格式（stdin）：
//     n C
//     w1 v1 c1
//     w2 v2 c2
//     …… （共 n 行）
// 輸出格式（stdout）：
//     第 1 行：多重背包的最大價值
//     第 2 行：多重背包的選取個數（n 個非負整數；n = 0 時輸出空行）
//     第 3 行：完全背包的最大價值
//     第 4 行：完全背包的選取個數（n 個非負整數）
//     第 5 行：湊出恰好總重 C 的最少件數（無解輸出 -1）
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static const int INF = 1000000000;

// ---------------------------------------------------------------- 多重背包
// 二維 DP（O(n·C·c)）：回傳 {最大價值, 每種物品取幾個}
// dp[i][cap] = 只考慮前 i 種、總重不超過 cap 的最大價值
static pair<int, vector<int>> multipleKnapsack2d(int nf, int C, const vector<int>& w,
                                                 const vector<int>& v, const vector<int>& c) {
    vector<vector<int>> dp(nf + 1, vector<int>(C + 1, 0));
    for (int i = 1; i <= nf; ++i) {
        int wi = w[i - 1], vi = v[i - 1], ci = c[i - 1];
        for (int cap = 0; cap <= C; ++cap) {
            int best = dp[i - 1][cap];                  // 取 0 個
            for (int t = 1; t <= ci && t * wi <= cap; ++t)
                best = max(best, dp[i - 1][cap - t * wi] + t * vi);
            dp[i][cap] = best;
        }
    }
    // 回溯：從第 nf 種往回，取滿足轉移式的最小 t
    vector<int> take(nf, 0);
    int cap = C;
    for (int i = nf; i >= 1; --i) {
        int wi = w[i - 1], vi = v[i - 1], ci = c[i - 1];
        int chosen = 0;
        for (int t = 0; t <= min(ci, cap / wi); ++t) {
            if (dp[i][cap] == dp[i - 1][cap - t * wi] + t * vi) {
                chosen = t;
                break;
            }
        }
        take[i - 1] = chosen;
        cap -= chosen * wi;
    }
    return {dp[nf][C], take};
}

// 二進制拆分 + 0-1 背包（O(n·C·log c)），只求最大價值（測試裡與 2D 版對拍）
static int multipleKnapsackSplit(int nf, int C, const vector<int>& w,
                                 const vector<int>& v, const vector<int>& c) {
    vector<int> dp(C + 1, 0);
    for (int i = 0; i < nf; ++i) {
        int wi = w[i], vi = v[i], rest = c[i], k = 1;
        while (rest > 0) {
            int take = min(k, rest);                    // 拆成 1, 2, 4, …, 剩餘
            int bw = take * wi, bv = take * vi;
            for (int cap = C; cap >= bw; --cap)         // 0-1 背包：容量倒序
                dp[cap] = max(dp[cap], dp[cap - bw] + bv);
            rest -= take;
            k <<= 1;
        }
    }
    return dp[C];
}

// ---------------------------------------------------------------- 完全背包
// 1D DP 容量正序（O(n·C)）：回傳 {最大價值, 每種物品取幾個}
static pair<int, vector<int>> completeKnapsack(int nf, int C, const vector<int>& w,
                                               const vector<int>& v) {
    vector<int> dp(C + 1, 0);
    vector<int> par(C + 1, -1);                         // par[cap] = 該容量最後被哪種物品更新
    for (int i = 0; i < nf; ++i) {
        int wi = w[i], vi = v[i];
        for (int cap = wi; cap <= C; ++cap) {           // 正序：允許同一種物品被重複取用
            int cand = dp[cap - wi] + vi;
            if (cand > dp[cap]) {
                dp[cap] = cand;
                par[cap] = i;
            }
        }
    }
    vector<int> take(nf, 0);
    int cap = C;
    while (cap > 0) {
        int i = par[cap];
        if (i < 0) break;
        take[i] += 1;
        cap -= w[i];
    }
    return {dp[C], take};
}

// 每種物品無限取，湊出恰好總重 C 的最少件數；無解回傳 -1
static int minItemsExact(int nf, int C, const vector<int>& w) {
    vector<int> dp(C + 1, INF);
    dp[0] = 0;
    for (int i = 0; i < nf; ++i) {
        int wi = w[i];
        for (int cap = wi; cap <= C; ++cap)
            if (dp[cap - wi] + 1 < dp[cap]) dp[cap] = dp[cap - wi] + 1;
    }
    return dp[C] >= INF ? -1 : dp[C];
}

// ---------------------------------------------------------------- 暴力基準（只用於測試）
// 枚舉所有 (t1, …, tn) 組合求最大價值
static int bruteBest(int nf, int C, const vector<int>& w, const vector<int>& v,
                     const vector<int>& c) {
    int best = 0;
    vector<int> t(nf, 0);
    auto eval = [&]() {
        int weight = 0, value = 0;
        for (int i = 0; i < nf; ++i) {
            weight += t[i] * w[i];
            value += t[i] * v[i];
        }
        if (weight <= C) best = max(best, value);
    };
    // 混合進位計數：第 i 位的進位上限是 c[i]
    while (true) {
        eval();
        int i = 0;
        while (i < nf) {
            t[i]++;
            if (t[i] <= c[i]) break;
            t[i] = 0;
            i++;
        }
        if (i == nf) break;
    }
    return best;
}

// Bellman-Ford 式鬆弛求最少件數
static int bruteMinItems(int nf, int C, const vector<int>& w) {
    vector<int> dp(C + 1, INF);
    dp[0] = 0;
    for (int round = 0; round < C; ++round) {
        bool changed = false;
        for (int cap = 1; cap <= C; ++cap)
            for (int i = 0; i < nf; ++i)
                if (cap - w[i] >= 0 && dp[cap - w[i]] + 1 < dp[cap]) {
                    dp[cap] = dp[cap - w[i]] + 1;
                    changed = true;
                }
        if (!changed) break;
    }
    return dp[C] >= INF ? -1 : dp[C];
}

// ---------------------------------------------------------------- 小工具
static void printVec(const vector<int>& a) {
    for (size_t i = 0; i < a.size(); ++i) {
        if (i) cout << ' ';
        cout << a[i];
    }
    cout << '\n';
}

static mt19937 rngEngine;
static int rndInt(int lo, int hi) { return lo + (int)(rngEngine() % (unsigned)(hi - lo + 1)); }

// ---------------------------------------------------------------- IO 模式
static void runIo(const string& data) {
    istringstream in(data);
    int n, C;
    if (!(in >> n >> C)) return;
    n = max(0, n);
    C = max(0, C);
    vector<int> wAll(n, 0), vAll(n, 0), cAll(n, 0);
    for (int i = 0; i < n; ++i) {
        // 逐個讀取：輸入截斷時缺的那個欄位補 0（與 Python 版行為一致），
        // 不能寫成 if (in >> a >> b >> d)，否則整行的前兩個數字也會被丟掉
        int a = 0, b = 0, d = 0;
        in >> a;
        in >> b;
        in >> d;
        wAll[i] = a;
        vAll[i] = b;
        cAll[i] = d;
    }

    vector<int> idx;
    for (int i = 0; i < n; ++i)
        if (wAll[i] >= 1) idx.push_back(i);             // 忽略 w = 0 的物品（價值無界）
    int nf = (int)idx.size();
    vector<int> w(nf), v(nf), c(nf);
    for (int j = 0; j < nf; ++j) {
        w[j] = wAll[idx[j]];
        v[j] = vAll[idx[j]];
        c[j] = max(0, cAll[idx[j]]);
    }

    pair<int, vector<int>> multi = multipleKnapsack2d(nf, C, w, v, c);
    vector<int> cntMulti(n, 0);
    for (int j = 0; j < nf; ++j) cntMulti[idx[j]] = multi.second[j];
    cout << multi.first << '\n';
    printVec(cntMulti);

    pair<int, vector<int>> full = completeKnapsack(nf, C, w, v);
    vector<int> cntFull(n, 0);
    for (int j = 0; j < nf; ++j) cntFull[idx[j]] = full.second[j];
    cout << full.first << '\n';
    printVec(cntFull);

    cout << minItemsExact(nf, C, w) << '\n';
}

// ---------------------------------------------------------------- 測試
static void runTests() {
    // README 示例：(2,3,2) (3,4,1) (4,5,3)，C = 10
    {
        vector<int> w = {2, 3, 4}, v = {3, 4, 5}, c = {2, 1, 3};
        int C = 10;
        pair<int, vector<int>> r = multipleKnapsack2d(3, C, w, v, c);
        assert(r.first == 13);                          // 1 個 w=2 + 2 個 w=4 → 重 10、值 13
        assert((r.second == vector<int>{1, 0, 2}));
        assert(r.first == multipleKnapsackSplit(3, C, w, v, c));
        assert(r.first == bruteBest(3, C, w, v, c));
        pair<int, vector<int>> f = completeKnapsack(3, C, w, v);
        assert(f.first == 15);                          // 5 個 w=2 → 重 10、值 15
        int wsum = 0, vsum = 0;
        for (int i = 0; i < 3; ++i) { wsum += f.second[i] * w[i]; vsum += f.second[i] * v[i]; }
        assert(wsum <= C && vsum == f.first);
        assert(minItemsExact(3, C, w) == 3);            // 4 + 4 + 2 = 10，三件
    }

    // 空輸入
    assert(multipleKnapsack2d(0, 0, {}, {}, {}).first == 0);
    assert(multipleKnapsackSplit(0, 0, {}, {}, {}) == 0);
    assert(completeKnapsack(0, 0, {}, {}).first == 0);
    assert(minItemsExact(0, 0, {}) == 0);
    assert(minItemsExact(1, 5, {3}) == -1);             // 3 湊不出 5

    // 容量 0：價值 0、件數 0
    {
        pair<int, vector<int>> r = multipleKnapsack2d(2, 0, {1, 2}, {5, 9}, {3, 3});
        assert(r.first == 0 && (r.second == vector<int>{0, 0}));
        pair<int, vector<int>> f = completeKnapsack(2, 0, {1, 2}, {5, 9});
        assert(f.first == 0 && (f.second == vector<int>{0, 0}));
        assert(minItemsExact(2, 0, {1, 2}) == 0);
    }

    // 數量上限為 0 時該物品完全不能用
    assert(multipleKnapsack2d(2, 5, {1, 2}, {5, 9}, {0, 3}).first == 18);   // 只能用第二種：2 個 w=2

    // 裝不下任何東西
    assert(multipleKnapsack2d(1, 3, {5}, {100}, {1}) == make_pair(0, vector<int>{0}));
    assert(completeKnapsack(1, 3, {5}, {100}) == make_pair(0, vector<int>{0}));

    // 隨機對拍：2D 版 vs 二進制拆分版 vs 暴力枚舉
    rngEngine.seed(20261002);
    for (int t = 0; t < 400; ++t) {
        int nf = rndInt(0, 4), C = rndInt(0, 12);
        vector<int> w(nf), v(nf), c(nf);
        for (int i = 0; i < nf; ++i) {
            w[i] = rndInt(1, 5);
            v[i] = rndInt(0, 9);
            c[i] = rndInt(0, 3);
        }

        pair<int, vector<int>> r = multipleKnapsack2d(nf, C, w, v, c);
        assert(r.first == multipleKnapsackSplit(nf, C, w, v, c));
        assert(r.first == bruteBest(nf, C, w, v, c));
        int wsum = 0, vsum = 0;
        for (int i = 0; i < nf; ++i) {
            assert(r.second[i] >= 0 && r.second[i] <= c[i]);
            wsum += r.second[i] * w[i];
            vsum += r.second[i] * v[i];
        }
        assert(wsum <= C);
        assert(vsum == r.first);

        // 完全背包 = 把數量上限設成「最多能裝幾個」之後的多重背包
        vector<int> capC(nf);
        for (int i = 0; i < nf; ++i) capC[i] = C / w[i];
        pair<int, vector<int>> f = completeKnapsack(nf, C, w, v);
        assert(f.first == bruteBest(nf, C, w, v, capC));
        wsum = 0; vsum = 0;
        for (int i = 0; i < nf; ++i) {
            wsum += f.second[i] * w[i];
            vsum += f.second[i] * v[i];
        }
        assert(wsum <= C);
        assert(vsum == f.first);

        assert(minItemsExact(nf, C, w) == bruteMinItems(nf, C, w));
    }

    cout << "all tests passed" << '\n';
}

int main() {
    string data, line;
    while (getline(cin, line)) {
        data += line;
        data += '\n';
    }
    // 只有含非空白內容才走 IO 模式（與 Python 版的 raw.strip() 對齊）
    if (data.find_first_not_of(" \t\r\n") != string::npos) runIo(data);
    else runTests();
    return 0;
}
