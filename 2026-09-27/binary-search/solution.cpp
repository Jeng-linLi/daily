// 二分查找與二分答案（lower_bound / upper_bound / 精確查找 / 最小化最大值）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：二分查找的本質是**在一個單調的判定函數上找分界點**。用左閉右開區間
//   [lo, hi) 維護「還沒確定」的部分，每次取中點把區間砍一半：
//     - lower_bound：判定 P(i) = (a[i] >= x)，找第一個使 P 爲真的 i；
//     - upper_bound：判定 P(i) = (a[i] >  x)，找第一個使 P 爲真的 i。
//   兩者只差一個比較符號（< 與 <=），這是最容易寫錯的地方。有了它們之後，
//   精確查找 = lower_bound 取出位置後判等；元素個數 = upper - lower。
//
//   二分答案同理，判定函數換成「給定上限 limit，能否切成 <= k 段且每段和 <= limit」。
//   可行性關於 limit 單調，故可二分最小可行值。下界 max(a)，上界 sum(a)。
//
// 輸入（空白分隔）：n / a1..an（升序，n=0 時省略）/ q / x1..xq
// 輸出：每個查詢一行 `<lower> <upper>`
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

// 第一個 >= x 的下標；不存在返回 n。時間 O(log n)，空間 O(1)
int lowerBound(const vector<int>& a, int x) {
    int lo = 0, hi = static_cast<int>(a.size());
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] < x) lo = mid + 1;   // a[mid] 太小，答案在右半邊
        else hi = mid;                  // a[mid] >= x，mid 本身可能是答案
    }
    return lo;
}

// 第一個 > x 的下標；不存在返回 n。時間 O(log n)，空間 O(1)
int upperBound(const vector<int>& a, int x) {
    int lo = 0, hi = static_cast<int>(a.size());
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] <= x) lo = mid + 1;  // 與 lowerBound 唯一的差別：相等也算「太小」
        else hi = mid;
    }
    return lo;
}

// 精確查找：返回第一個等於 x 的下標，不存在返回 -1
int binarySearch(const vector<int>& a, int x) {
    int i = lowerBound(a, x);
    return (i < static_cast<int>(a.size()) && a[i] == x) ? i : -1;
}

// 等於 x 的元素個數
int countEqual(const vector<int>& a, int x) {
    return upperBound(a, x) - lowerBound(a, x);
}

// 貪心判定：能否切成不超過 k 段且每段和 <= limit。時間 O(n)
bool canSplit(const vector<int>& a, int k, int limit) {
    int parts = 1, cur = 0;
    for (int v : a) {
        if (v > limit) return false;  // 單個元素就超過 limit，這一段無論如何裝不下
        if (cur + v <= limit) {
            cur += v;               // 還能裝下，繼續往當前段裏塞
        } else {
            ++parts;                // 裝不下了，在這裡切一刀
            cur = v;
            if (parts > k) return false;
        }
    }
    return true;
}

// 切成 k 個非空連續段，最小化最大段的和。時間 O(n log S)，S = sum(a) - max(a)
int minMaxSplit(const vector<int>& a, int k) {
    if (a.empty()) return 0;
    if (k > static_cast<int>(a.size())) k = static_cast<int>(a.size());  // 段數超過元素個數沒意義
    int lo = a[0], hi = 0;
    for (int v : a) {
        if (v > lo) lo = v;
        hi += v;
    }
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (canSplit(a, k, mid)) hi = mid;   // mid 可行，答案 <= mid
        else lo = mid + 1;                   // mid 不可行，答案 > mid
    }
    return lo;
}

// 對照用的 O(n^2 * k) 動態規劃，僅用於小規模測試驗證
int minMaxSplitBrute(const vector<int>& a, int k) {
    int n = static_cast<int>(a.size());
    if (n == 0) return 0;
    if (k > n) k = n;
    vector<int> pre(n + 1, 0);
    for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + a[i];
    const int INF = 1e9;
    vector<vector<int>> dp(n + 1, vector<int>(k + 1, INF));
    dp[0][0] = 0;
    for (int i = 1; i <= n; ++i) {
        for (int p = 1; p <= k && p <= i; ++p) {
            int best = INF;
            for (int j = p - 1; j < i; ++j) {         // 最後一段是 (j, i]
                if (dp[j][p - 1] == INF) continue;
                int val = max(dp[j][p - 1], pre[i] - pre[j]);
                if (val < best) best = val;
            }
            dp[i][p] = best;
        }
    }
    return dp[n][k];
}

// 線性掃描版，用於對拍
int lowerBoundBrute(const vector<int>& a, int x) {
    for (int i = 0; i < static_cast<int>(a.size()); ++i)
        if (a[i] >= x) return i;
    return static_cast<int>(a.size());
}

int upperBoundBrute(const vector<int>& a, int x) {
    for (int i = 0; i < static_cast<int>(a.size()); ++i)
        if (a[i] > x) return i;
    return static_cast<int>(a.size());
}

// 與 Python 版同規模的固定隨機序列（LCG），兩版各自獨立與暴力解對拍
struct LCG {
    unsigned long long s;
    LCG(unsigned long long seed) : s(seed) {}
    int next(int lo, int hi) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<int>((s >> 33) % static_cast<unsigned long long>(hi - lo + 1));
    }
};

int main() {
    int n;
    if (cin >> n) {  // IO 模式
        vector<int> a(n);
        for (int i = 0; i < n; ++i) {
            if (!(cin >> a[i])) a[i] = 0;      // 輸入被截斷時用 0 兜底
        }
        int q = 0;
        if (!(cin >> q)) q = 0;
        for (int i = 0; i < q; ++i) {
            int x = 0;
            if (!(cin >> x)) x = 0;
            cout << lowerBound(a, x) << " " << upperBound(a, x) << "\n";
        }
        return 0;
    }

    // README 示例：a = [1, 2, 2, 2, 4, 7]
    {
        vector<int> a = {1, 2, 2, 2, 4, 7};
        assert(lowerBound(a, 2) == 1);         // 第一個 >= 2 的位置
        assert(upperBound(a, 2) == 4);         // 第一個 > 2 的位置
        assert(countEqual(a, 2) == 3);         // 一共 3 個 2
        assert(binarySearch(a, 2) == 1);
        assert(binarySearch(a, 3) == -1);      // 3 不存在
        assert(countEqual(a, 3) == 0);

        // 邊界：比所有元素都小 / 都大
        assert(lowerBound(a, 0) == 0);
        assert(upperBound(a, 0) == 0);
        assert(lowerBound(a, 9) == 6);
        assert(upperBound(a, 9) == 6);
        assert(countEqual(a, 9) == 0);

        // 落在兩個元素之間的空隙裏
        assert(lowerBound(a, 3) == 4);
        assert(upperBound(a, 3) == 4);
        assert(lowerBound(a, 5) == 5);
        assert(upperBound(a, 5) == 5);

        // 命中唯一元素 / 命中最大元素
        assert(lowerBound(a, 1) == 0);
        assert(upperBound(a, 1) == 1);
        assert(lowerBound(a, 7) == 5);
        assert(upperBound(a, 7) == 6);
    }

    // 空數組：任何查詢都返回 0
    {
        vector<int> e;
        assert(lowerBound(e, 5) == 0);
        assert(upperBound(e, 5) == 0);
        assert(binarySearch(e, 5) == -1);
        assert(countEqual(e, 5) == 0);
    }

    // 單元素數組
    {
        vector<int> one = {5};
        assert(lowerBound(one, 5) == 0);
        assert(upperBound(one, 5) == 1);
        assert(binarySearch(one, 4) == -1);
    }

    // 二分答案：LeetCode 410 的經典用例
    assert(minMaxSplit({7, 2, 5, 10, 8}, 2) == 18);   // [7,2,5] | [10,8]
    assert(minMaxSplit({1, 2, 3, 4, 5}, 2) == 9);     // [1,2,3] | [4,5]
    assert(minMaxSplit({1, 4, 4}, 3) == 4);           // 每段一個元素
    assert(minMaxSplit({7, 2, 5, 10, 8}, 5) == 10);   // 段數 >= n 時就是 max(a)
    assert(minMaxSplit({7, 2, 5, 10, 8}, 9) == 10);
    assert(minMaxSplit({}, 3) == 0);                  // 空數組
    assert(minMaxSplit({5}, 1) == 5);

    LCG rng(20260927ULL);

    // 隨機對拍一：lower_bound / upper_bound 與線性掃描逐一比對
    for (int t = 0; t < 300; ++t) {
        int n = rng.next(0, 12);
        vector<int> a;
        for (int i = 0; i < n; ++i) a.push_back(rng.next(0, 10));
        sort(a.begin(), a.end());                     // 升序，故意製造重複元素
        for (int x = -2; x <= 12; ++x) {
            assert(lowerBound(a, x) == lowerBoundBrute(a, x));
            assert(upperBound(a, x) == upperBoundBrute(a, x));
            assert(0 <= lowerBound(a, x));
            assert(lowerBound(a, x) <= upperBound(a, x));
            assert(upperBound(a, x) <= n);
            // 精確查找與計數要和 lower/upper 自洽
            if (binarySearch(a, x) == -1) {
                assert(countEqual(a, x) == 0);
                assert(lowerBound(a, x) == upperBound(a, x));
            } else {
                assert(a[binarySearch(a, x)] == x);
                assert(binarySearch(a, x) == lowerBound(a, x));
                assert(countEqual(a, x) == upperBound(a, x) - lowerBound(a, x));
            }
        }
    }

    // 隨機對拍二：二分答案與 O(n^2*k) 的 DP 暴力解比對
    for (int t = 0; t < 200; ++t) {
        int n = rng.next(0, 8);
        vector<int> a;
        for (int i = 0; i < n; ++i) a.push_back(rng.next(0, 20));  // 不要求有序，非負即可
        for (int k = 1; k <= n; ++k) {
            int got = minMaxSplit(a, k);
            int exp = minMaxSplitBrute(a, k);
            assert(got == exp);
            assert(canSplit(a, k, got));               // 最優值確實可行
            if (got > 0) assert(!canSplit(a, k, got - 1));  // 再小一點就不行了
        }
        if (n == 0) assert(minMaxSplit(a, 1) == 0);
    }

    cout << "all tests passed" << endl;
    return 0;
}
