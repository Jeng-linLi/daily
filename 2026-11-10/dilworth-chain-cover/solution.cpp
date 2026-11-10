// Dilworth 定理：二維偏序集的最小鏈覆蓋 = 最長反鏈（C++17）
//
// 命令（IO 模式，從 stdin 讀入）：
//   chain n x1 y1 x2 y2 ... xn yn
//       -> 輸出「最小鏈覆蓋 最長反鏈」兩個數（依 Dilworth 定理恆相等）
//
// 偏序定義：(a, b) ≼ (c, d) 當且僅當 a <= c 且 b <= d（支配序）。
//   鏈 = 兩兩可比的點集；反鏈 = 兩兩互不支配的點集。
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * 去重後按 a 升序排序，則反鏈 = b 的最長嚴格遞減子序列（LDS，O(n^2) DP）。
//   * 最小鏈覆蓋 = DAG 最小路徑覆蓋 = n − 最大匹配（匈牙利算法）。
//   * 兩個完全獨立的算法（LDS vs 匹配）結果恆相等——Dilworth 定理的構造性驗證。

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
typedef long long ll;
typedef pair<ll, ll> pll;

vector<pll> dedupSort(vector<pll> pts) {
    set<pll> s(pts.begin(), pts.end());
    return vector<pll>(s.begin(), s.end());
}

int longestAntichain(const vector<pll>& pts) {
    int n = (int)pts.size();
    if (n == 0) return 0;
    vector<int> dp(n, 1);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            if (pts[j].second > pts[i].second)  // 嚴格遞減
                dp[i] = max(dp[i], dp[j] + 1);
    return *max_element(dp.begin(), dp.end());
}

// 匈牙利算法：左側 n 點，adj[i] = 右側鄰居
int hungarian(int n, const vector<vector<int>>& adj) {
    vector<int> matchR(n, -1);

    struct TryK {
        const vector<vector<int>>& adj;
        vector<int>& matchR;
        vector<char> vis;
        bool go(int i) {
            for (int j : adj[i]) {
                if (vis[j]) continue;
                vis[j] = 1;
                if (matchR[j] == -1 || go(matchR[j])) {
                    matchR[j] = i;
                    return true;
                }
            }
            return false;
        }
    };

    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        TryK tk{adj, matchR, vector<char>(n, 0)};
        if (tk.go(i)) ++cnt;
    }
    return cnt;
}

int minChainCover(const vector<pll>& pts) {
    int n = (int)pts.size();
    if (n == 0) return 0;
    vector<vector<int>> adj(n);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (pts[i].second <= pts[j].second)  // a 已升序，b_i <= b_j 即支配
                adj[i].push_back(j);
    return n - hungarian(n, adj);
}

pair<int, int> solve(vector<pll> pts) {
    vector<pll> v = dedupSort(pts);
    return {minChainCover(v), longestAntichain(v)};
}

// 自製 LCG
static ll g_state = 20261110;
ll rnd(ll m) {
    g_state = (g_state * 1103515245 + 12345) & 0x7FFFFFFFLL;
    return g_state % m;
}

void runTests() {
    // 1) 全鏈 -> 覆蓋 1、反鏈 1
    assert(solve({{1, 1}, {2, 2}, {3, 3}}) == make_pair(1, 1));
    // 2) 全反鏈 -> 覆蓋 = 反鏈 = n
    assert(solve({{1, 5}, {2, 4}, {3, 3}, {4, 2}, {5, 1}}) == make_pair(5, 5));
    // 3) 手驗混合例
    assert(solve({{1, 1}, {1, 3}, {2, 2}, {3, 4}}) == make_pair(2, 2));
    // 4) 重複點去重
    assert(solve({{1, 1}, {1, 1}, {2, 2}, {1, 1}}) == make_pair(1, 1));
    // 5) 空集
    assert(solve({}) == make_pair(0, 0));

    // 6) Dilworth 定理驗證：LDS 與匹配兩獨立算法恆相等
    for (int t = 0; t < 200; ++t) {
        int n = (int)rnd(20) + 1;
        vector<pll> pts;
        for (int i = 0; i < n; ++i)
            pts.push_back({rnd(15) + 1, rnd(15) + 1});
        auto res = solve(pts);
        assert(res.first == res.second);
    }

    // 7) LDS 與暴力子集枚舉對拍
    auto bruteAnti = [](vector<pll> pts) {
        int n = (int)pts.size(), best = 0;
        for (int mask = 0; mask < (1 << n); ++mask) {
            vector<pll> sel;
            for (int i = 0; i < n; ++i)
                if (mask >> i & 1) sel.push_back(pts[i]);
            bool ok = true;
            for (size_t a = 0; a < sel.size() && ok; ++a)
                for (size_t b = a + 1; b < sel.size() && ok; ++b) {
                    const pll& p = sel[a];
                    const pll& q = sel[b];
                    if ((p.first <= q.first && p.second <= q.second) ||
                        (q.first <= p.first && q.second <= p.second))
                        ok = false;
                }
            if (ok) best = max(best, (int)sel.size());
        }
        return best;
    };
    for (int t = 0; t < 50; ++t) {
        int n = (int)rnd(11) + 2;
        vector<pll> pts;
        for (int i = 0; i < n; ++i)
            pts.push_back({rnd(6), rnd(6)});
        vector<pll> v = dedupSort(pts);
        assert(longestAntichain(v) == bruteAnti(v));
    }

    // 8) 匈牙利與暴力最大匹配對拍
    auto bruteMatch = [](int n, const vector<vector<int>>& adj) {
        vector<int> matchR(n, -1);
        struct Go {
            int n, best;
            const vector<vector<int>>& adj;
            vector<int>& matchR;
            void go(int i) {
                if (i == n) {
                    int cnt = 0;
                    for (int x : matchR) if (x != -1) ++cnt;
                    best = max(best, cnt);
                    return;
                }
                go(i + 1);
                for (int j : adj[i]) {
                    if (matchR[j] == -1) {
                        matchR[j] = i;
                        go(i + 1);
                        matchR[j] = -1;
                    }
                }
            }
        };
        Go g{n, 0, adj, matchR};
        g.go(0);
        return g.best;
    };
    for (int t = 0; t < 50; ++t) {
        int n = (int)rnd(9) + 1;
        vector<pll> pts;
        for (int i = 0; i < n; ++i)
            pts.push_back({rnd(8), rnd(8)});
        vector<pll> v = dedupSort(pts);
        int nn = (int)v.size();
        vector<vector<int>> adj(nn);
        for (int i = 0; i < nn; ++i)
            for (int j = i + 1; j < nn; ++j)
                if (v[i].second <= v[j].second) adj[i].push_back(j);
        assert(hungarian(nn, adj) == bruteMatch(nn, adj));
    }
}

int main() {
    string cmd;
    if (cin >> cmd) {
        bool any = true;
        while (any) {
            if (cmd == "chain") {
                int n;
                if (!(cin >> n)) break;
                vector<pll> pts;
                bool ok = true;
                for (int i = 0; i < n; ++i) {
                    ll x, y;
                    if (!(cin >> x >> y)) { ok = false; break; }
                    pts.push_back({x, y});
                }
                if (!ok) break;
                auto res = solve(pts);
                cout << res.first << ' ' << res.second << '\n';
            } else {
                break;  // 未知命令
            }
            if (!(cin >> cmd)) break;
        }
    } else {
        runTests();
        printf("all tests passed\n");
    }
    return 0;
}
