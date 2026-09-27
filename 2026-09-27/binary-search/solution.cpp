// 二分查找与二分答案（lower_bound / upper_bound / 精确查找 / 最小化最大值）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：二分查找的本质是**在一个单调的判定函数上找分界点**。用左闭右开区间
//   [lo, hi) 维护「还没确定」的部分，每次取中点把区间砍一半：
//     - lower_bound：判定 P(i) = (a[i] >= x)，找第一个使 P 为真的 i；
//     - upper_bound：判定 P(i) = (a[i] >  x)，找第一个使 P 为真的 i。
//   两者只差一个比较符号（< 与 <=），这是最容易写错的地方。有了它们之后，
//   精确查找 = lower_bound 取出位置后判等；元素个数 = upper - lower。
//
//   二分答案同理，判定函数换成「给定上限 limit，能否切成 <= k 段且每段和 <= limit」。
//   可行性关于 limit 单调，故可二分最小可行值。下界 max(a)，上界 sum(a)。
//
// 输入（空白分隔）：n / a1..an（升序，n=0 时省略）/ q / x1..xq
// 输出：每个查询一行 `<lower> <upper>`
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

// 第一个 >= x 的下标；不存在返回 n。时间 O(log n)，空间 O(1)
int lowerBound(const vector<int>& a, int x) {
    int lo = 0, hi = static_cast<int>(a.size());
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] < x) lo = mid + 1;   // a[mid] 太小，答案在右半边
        else hi = mid;                  // a[mid] >= x，mid 本身可能是答案
    }
    return lo;
}

// 第一个 > x 的下标；不存在返回 n。时间 O(log n)，空间 O(1)
int upperBound(const vector<int>& a, int x) {
    int lo = 0, hi = static_cast<int>(a.size());
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] <= x) lo = mid + 1;  // 与 lowerBound 唯一的差别：相等也算「太小」
        else hi = mid;
    }
    return lo;
}

// 精确查找：返回第一个等于 x 的下标，不存在返回 -1
int binarySearch(const vector<int>& a, int x) {
    int i = lowerBound(a, x);
    return (i < static_cast<int>(a.size()) && a[i] == x) ? i : -1;
}

// 等于 x 的元素个数
int countEqual(const vector<int>& a, int x) {
    return upperBound(a, x) - lowerBound(a, x);
}

// 贪心判定：能否切成不超过 k 段且每段和 <= limit。时间 O(n)
bool canSplit(const vector<int>& a, int k, int limit) {
    int parts = 1, cur = 0;
    for (int v : a) {
        if (v > limit) return false;  // 单个元素就超过 limit，这一段无论如何装不下
        if (cur + v <= limit) {
            cur += v;               // 还能装下，继续往当前段里塞
        } else {
            ++parts;                // 装不下了，在这里切一刀
            cur = v;
            if (parts > k) return false;
        }
    }
    return true;
}

// 切成 k 个非空连续段，最小化最大段的和。时间 O(n log S)，S = sum(a) - max(a)
int minMaxSplit(const vector<int>& a, int k) {
    if (a.empty()) return 0;
    if (k > static_cast<int>(a.size())) k = static_cast<int>(a.size());  // 段数超过元素个数没意义
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

// 对照用的 O(n^2 * k) 动态规划，仅用于小规模测试验证
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
            for (int j = p - 1; j < i; ++j) {         // 最后一段是 (j, i]
                if (dp[j][p - 1] == INF) continue;
                int val = max(dp[j][p - 1], pre[i] - pre[j]);
                if (val < best) best = val;
            }
            dp[i][p] = best;
        }
    }
    return dp[n][k];
}

// 线性扫描版，用于对拍
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

// 与 Python 版同规模的固定随机序列（LCG），两版各自独立与暴力解对拍
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
            if (!(cin >> a[i])) a[i] = 0;      // 输入被截断时用 0 兜底
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
        assert(lowerBound(a, 2) == 1);         // 第一个 >= 2 的位置
        assert(upperBound(a, 2) == 4);         // 第一个 > 2 的位置
        assert(countEqual(a, 2) == 3);         // 一共 3 个 2
        assert(binarySearch(a, 2) == 1);
        assert(binarySearch(a, 3) == -1);      // 3 不存在
        assert(countEqual(a, 3) == 0);

        // 边界：比所有元素都小 / 都大
        assert(lowerBound(a, 0) == 0);
        assert(upperBound(a, 0) == 0);
        assert(lowerBound(a, 9) == 6);
        assert(upperBound(a, 9) == 6);
        assert(countEqual(a, 9) == 0);

        // 落在两个元素之间的空隙里
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

    // 空数组：任何查询都返回 0
    {
        vector<int> e;
        assert(lowerBound(e, 5) == 0);
        assert(upperBound(e, 5) == 0);
        assert(binarySearch(e, 5) == -1);
        assert(countEqual(e, 5) == 0);
    }

    // 单元素数组
    {
        vector<int> one = {5};
        assert(lowerBound(one, 5) == 0);
        assert(upperBound(one, 5) == 1);
        assert(binarySearch(one, 4) == -1);
    }

    // 二分答案：LeetCode 410 的经典用例
    assert(minMaxSplit({7, 2, 5, 10, 8}, 2) == 18);   // [7,2,5] | [10,8]
    assert(minMaxSplit({1, 2, 3, 4, 5}, 2) == 9);     // [1,2,3] | [4,5]
    assert(minMaxSplit({1, 4, 4}, 3) == 4);           // 每段一个元素
    assert(minMaxSplit({7, 2, 5, 10, 8}, 5) == 10);   // 段数 >= n 时就是 max(a)
    assert(minMaxSplit({7, 2, 5, 10, 8}, 9) == 10);
    assert(minMaxSplit({}, 3) == 0);                  // 空数组
    assert(minMaxSplit({5}, 1) == 5);

    LCG rng(20260927ULL);

    // 随机对拍一：lower_bound / upper_bound 与线性扫描逐一比对
    for (int t = 0; t < 300; ++t) {
        int n = rng.next(0, 12);
        vector<int> a;
        for (int i = 0; i < n; ++i) a.push_back(rng.next(0, 10));
        sort(a.begin(), a.end());                     // 升序，故意制造重复元素
        for (int x = -2; x <= 12; ++x) {
            assert(lowerBound(a, x) == lowerBoundBrute(a, x));
            assert(upperBound(a, x) == upperBoundBrute(a, x));
            assert(0 <= lowerBound(a, x));
            assert(lowerBound(a, x) <= upperBound(a, x));
            assert(upperBound(a, x) <= n);
            // 精确查找与计数要和 lower/upper 自洽
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

    // 随机对拍二：二分答案与 O(n^2*k) 的 DP 暴力解比对
    for (int t = 0; t < 200; ++t) {
        int n = rng.next(0, 8);
        vector<int> a;
        for (int i = 0; i < n; ++i) a.push_back(rng.next(0, 20));  // 不要求有序，非负即可
        for (int k = 1; k <= n; ++k) {
            int got = minMaxSplit(a, k);
            int exp = minMaxSplitBrute(a, k);
            assert(got == exp);
            assert(canSplit(a, k, got));               // 最优值确实可行
            if (got > 0) assert(!canSplit(a, k, got - 1));  // 再小一点就不行了
        }
        if (n == 0) assert(minMaxSplit(a, 1) == 0);
    }

    cout << "all tests passed" << endl;
    return 0;
}
