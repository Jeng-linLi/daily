// 蓄水池抽樣 Reservoir Sampling（Algorithm R）—— 流式均勻採樣（C++17）
//
// 命令（IO 模式，從 stdin 讀入）：
//   sample n k seed
//       -> 對 1..n 的流做蓄水池抽樣（隨機決策由自製 LCG 產生），
//          輸出被選中的 k 個編號（升序，空格分隔）。
//
// 應用場景：A/B 測試分組、日誌流抽樣、大數據集單遍均勻子採樣——
// 只需 O(k) 內存、單遍掃描、事先不知道流長 n 也能保證每個元素
// 被選中概率恰為 k/n。
//
// 算法（Algorithm R）：前 k 個元素直接入池；第 i 個元素（i > k）以
// k/i 的概率替換池中隨機一個位置。歸納可證池中恆為均勻樣本。
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <set>

using namespace std;
typedef long long ll;

struct LCG {
    // 自製線性同餘產生器（與倉庫其他題目一致，兩語言行為完全相同）
    ll state;
    explicit LCG(ll seed) : state(seed & 0x7FFFFFFFLL) {}
    ll nextRand(ll m) {
        state = (state * 1103515245 + 12345) & 0x7FFFFFFFLL;
        return state % m;
    }
};

vector<int> reservoirSample(int n, int k, ll seed) {
    LCG rng(seed);
    vector<int> res;
    for (int i = 1; i <= n; ++i) {
        if (i <= k) {
            res.push_back(i);
        } else if (rng.nextRand(i) < k) {
            res[(size_t)rng.nextRand(k)] = i;
        }
    }
    return res;
}

vector<int> solve(int n, int k, ll seed) {
    if (n <= 0 || k <= 0) return {};
    if (k >= n) {
        vector<int> all(n);
        for (int i = 0; i < n; ++i) all[i] = i + 1;
        return all;
    }
    vector<int> r = reservoirSample(n, k, seed);
    sort(r.begin(), r.end());
    return r;
}

// 自製 LCG（測試用）
static ll g_state = 20261110;
ll rnd(ll m) {
    g_state = (g_state * 1103515245 + 12345) & 0x7FFFFFFFLL;
    return g_state % m;
}

void runTests() {
    // 1) 結構性質：恰好 k 個、無重複、升序
    for (int t = 0; t < 200; ++t) {
        int n = (int)rnd(50) + 1;
        int k = (int)rnd(n + 1);  // k 可以為 0 或 n
        ll seed = rnd(1000000);
        vector<int> r = solve(n, k, seed);
        assert((int)r.size() == min(k, n));
        assert(is_sorted(r.begin(), r.end()));
        assert(adjacent_find(r.begin(), r.end()) == r.end());  // 無重複
        for (int x : r) assert(1 <= x && x <= n);
    }

    // 2) 邊界：k=0 空、k>=n 全選、n=1
    assert(solve(10, 0, 1).empty());
    {
        vector<int> full = solve(10, 10, 42);
        for (int i = 0; i < 10; ++i) assert(full[i] == i + 1);
        full = solve(10, 99, 42);
        for (int i = 0; i < 10; ++i) assert(full[i] == i + 1);
    }
    assert(solve(1, 1, 7) == vector<int>{1});
    assert(solve(0, 0, 7).empty());

    // 3) 確定性：同一 seed 輸出一致
    assert(solve(20, 5, 12345) == solve(20, 5, 12345));

    // 4) 均勻性 Monte Carlo：n=10, k=3，每個元素期望 30% 入選
    //    20000 次實驗 -> 每個編號期望 6000 次，二項 std 約 65，容差 ±10% = 9σ
    const int trials = 20000;
    vector<ll> counts(11, 0);  // 1..10
    for (int t = 0; t < trials; ++t)
        for (int x : reservoirSample(10, 3, t + 1)) ++counts[x];
    for (int x = 1; x <= 10; ++x)
        assert(5400 <= counts[x] && counts[x] <= 6600);

    // 5) 組合多樣性：120 種 3-子集應全部出現過（期望每種 166.7 次）
    set<vector<int>> seen;
    for (int t = 0; t < trials; ++t) {
        vector<int> r = reservoirSample(10, 3, t + 1000000);
        sort(r.begin(), r.end());
        seen.insert(r);
    }
    assert((int)seen.size() == 120);
}

int main() {
    string cmd;
    if (cin >> cmd) {
        bool any = true;
        while (any) {
            if (cmd == "sample") {
                ll n, k, seed;
                if (!(cin >> n >> k >> seed)) break;
                vector<int> r = solve((int)n, (int)k, seed);
                for (size_t i = 0; i < r.size(); ++i) {
                    if (i) cout << ' ';
                    cout << r[i];
                }
                cout << '\n';
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
