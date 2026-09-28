// 0-1 背包（Knapsack 0-1，動態規劃）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：定義 dp[i][c] = 只考慮前 i 個物品、容量爲 c 時的最大價值。對第 i 個物品：
//     不選 dp[i][c] = dp[i-1][c]；選 dp[i][c] = dp[i-1][c-w[i]] + v[i]（需 c >= w[i]），取最大。
//   物品只能選一次，轉移只依賴 i-1 層，不會同層自我疊加（區別於完全背包）。
//
//   空間可壓到一維 dp[c]，但**容量 c 必須倒序遍歷**（W -> w[i]）：
//   倒序保證 dp[c-w[i]] 讀到的仍是 i-1 層舊值；正序則會把同一物品重複放入，退化成完全背包。
//   要還原方案必須保留二維表：dp[i][c] != dp[i-1][c] 說明第 i 個物品被選了，反推即可。
//
// 輸入：第一行 n W；接下來 n 行，每行 w_i v_i
// 輸出：第一行最大總價值；第二行被選中物品下標（0-based，升序，空格分隔）
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

// 只求最大價值：一維滾動數組，時間 O(n*W)，空間 O(W)
int knapsackMax(const vector<int>& w, const vector<int>& v, int capacity) {
    vector<int> dp(capacity + 1, 0);
    for (size_t i = 0; i < w.size(); ++i) {
        // 容量倒序；正序會讓同一物品被重複選取，退化爲完全背包
        for (int c = capacity; c >= w[i]; --c) {
            int cand = dp[c - w[i]] + v[i];
            if (cand > dp[c]) dp[c] = cand;
        }
    }
    return dp[capacity];
}

// 求最大價值並還原一種選取方案：保留二維表，空間 O(n*W)
pair<int, vector<int>> knapsackWithItems(const vector<int>& w, const vector<int>& v, int capacity) {
    int n = static_cast<int>(w.size());
    vector<vector<int>> dp(n + 1, vector<int>(capacity + 1, 0));

    for (int i = 1; i <= n; ++i) {
        for (int c = 0; c <= capacity; ++c) {
            int best = dp[i - 1][c];                       // 不選第 i-1 個物品
            if (c >= w[i - 1]) {
                int take = dp[i - 1][c - w[i - 1]] + v[i - 1];  // 選它
                if (take > best) best = take;
            }
            dp[i][c] = best;
        }
    }

    // 反向還原：dp[i][c] 比 dp[i-1][c] 大，說明第 i-1 個物品被選中了
    vector<int> chosen;
    int c = capacity;
    for (int i = n; i >= 1; --i) {
        if (dp[i][c] != dp[i - 1][c]) {
            chosen.push_back(i - 1);
            c -= w[i - 1];
        }
    }
    reverse(chosen.begin(), chosen.end());
    return {dp[n][capacity], chosen};
}

// 對照用的指數級枚舉，僅用於小規模測試驗證
int bruteForce(const vector<int>& w, const vector<int>& v, int capacity) {
    int n = static_cast<int>(w.size());
    int best = 0;
    for (int mask = 0; mask < (1 << n); ++mask) {
        int tw = 0, tv = 0;
        for (int i = 0; i < n; ++i) {
            if (mask >> i & 1) { tw += w[i]; tv += v[i]; }
        }
        if (tw <= capacity && tv > best) best = tv;
    }
    return best;
}

// 與 Python 版完全一致的固定隨機序列，保證兩版跑同一批用例
struct LCG {
    unsigned long long s;
    LCG(unsigned long long seed) : s(seed) {}
    int next(int lo, int hi) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<int>((s >> 33) % static_cast<unsigned long long>(hi - lo + 1));
    }
};

int main() {
    int n, W;
    if (cin >> n >> W) {  // IO 模式
        vector<int> w(n), v(n);
        for (int i = 0; i < n; ++i) cin >> w[i] >> v[i];
        auto res = knapsackWithItems(w, v, W);
        vector<int>& chosen = res.second;
        int total = 0;
        for (int i : chosen) total += v[i];
        cout << total << "\n";
        for (size_t i = 0; i < chosen.size(); ++i) {
            if (i) cout << " ";
            cout << chosen[i];
        }
        cout << "\n";
        return 0;
    }

    // 經典用例：容量 10 的最優值爲 12，但存在多解
    //   (a) 2 號 + 3 號：w=4+6=10, v=5+7=12
    //   (b) 0 號 + 1 號 + 2 號：w=2+3+4=9, v=3+4+5=12
    // 反推時「dp[i][c] == dp[i-1][c] 視爲未選」，會優先得到 (b)，故只斷言值最優、方案合法
    {
        vector<int> w = {2, 3, 4, 6};
        vector<int> v = {3, 4, 5, 7};
        auto res = knapsackWithItems(w, v, 10);
        assert(res.first == 12);
        assert(knapsackMax(w, v, 10) == 12);
        int tw = 0, tv = 0;
        for (int i : res.second) { tw += w[i]; tv += v[i]; }
        assert(tw <= 10);
        assert(tv == 12);
    }

    // 容量爲 0 / 物品爲空
    assert(knapsackMax({}, {}, 10) == 0);
    assert(knapsackMax({5}, {9}, 0) == 0);
    assert(knapsackWithItems({5}, {9}, 0).first == 0);
    assert(knapsackWithItems({5}, {9}, 0).second.empty());

    // 單件裝不下
    assert(knapsackMax({5}, {9}, 4) == 0);

    // 全部都能裝下
    assert(knapsackMax({1, 2, 3}, {1, 2, 3}, 10) == 6);

    // 按價值密度貪心會選錯：密度最高的是 0 號(2.0)，但最優解是 1+2 號 w=10 v=18
    assert(knapsackMax({5, 4, 6}, {10, 7, 11}, 10) == 18);
    assert(knapsackMax({5, 4, 6}, {10, 7, 11}, 10) == bruteForce({5, 4, 6}, {10, 7, 11}, 10));

    // 零重量物品：價值白拿，且不會造成死循環
    assert(knapsackMax({0, 3}, {5, 4}, 3) == 9);

    // 與暴力枚舉隨機對照：同時校驗最優值一致、方案可行且達到最優值
    LCG rng(20260922ULL);
    for (int t = 0; t < 300; ++t) {
        int len = rng.next(1, 10);
        int cap = rng.next(0, 20);
        vector<int> ws(len), vs(len);
        for (int i = 0; i < len; ++i) {
            ws[i] = rng.next(0, 8);
            vs[i] = rng.next(0, 20);
        }
        auto res = knapsackWithItems(ws, vs, cap);
        assert(res.first == bruteForce(ws, vs, cap));   // 最優值與暴力一致
        assert(knapsackMax(ws, vs, cap) == res.first);  // 一維版與二維版一致
        int tw = 0, tv = 0;
        for (int i : res.second) { tw += ws[i]; tv += vs[i]; }
        assert(tw <= cap);                              // 方案不超容量
        assert(tv == res.first);                        // 方案確實達到最優值
    }

    cout << "all tests passed" << endl;
    return 0;
}
