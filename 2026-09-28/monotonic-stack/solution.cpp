// 單調棧（下一個更大元素 / 每日溫度 / 柱狀圖最大矩形）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：單調棧解決「爲每個元素找左/右邊第一個滿足某種大小關係的元素」。
//   棧裏始終保持單調序列，新元素入棧前先把被它破壞單調性的元素彈出去，
//   而那些被彈出的元素，答案恰好就是當前這個新元素。每個下標入棧出棧各一次 → O(n)。
//   1) 下一個更大元素：棧中值單調遞減，遇到更大的就把小的彈掉，答案記爲當前下標；
//   2) 每日溫度：同源，寫進答案的是距離（下標差）；
//   3) 柱狀圖最大矩形：單調遞增棧 + 末尾高度 0 的哨兵，彈棧時用
//      「寬 = i - 新棧頂 - 1」結算以該柱爲高的最大矩形。
//
// 輸入（空白分隔）：n / a1..an（n=0 時省略）
// 輸出：第 1 行下一個更大元素下標；第 2 行每日溫度；第 3 行最大矩形面積
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

// 每個位置右邊第一個嚴格大於它的元素下標；沒有則 -1。時間 O(n)，空間 O(n)
vector<int> nextGreaterIndex(const vector<int>& a) {
    int n = static_cast<int>(a.size());
    vector<int> res(n, -1);
    vector<int> st;                       // 存下標，對應值單調遞減
    for (int i = 0; i < n; ++i) {
        while (!st.empty() && a[st.back()] < a[i]) {
            res[st.back()] = i;           // a[i] 就是它們右邊第一個更大的
            st.pop_back();
        }
        st.push_back(i);
    }
    return res;
}

// 與 nextGreaterIndex 同源，但輸出距離而非下標；沒有則 0。時間 O(n)，空間 O(n)
vector<int> dailyTemperatures(const vector<int>& a) {
    int n = static_cast<int>(a.size());
    vector<int> res(n, 0);
    vector<int> st;
    for (int i = 0; i < n; ++i) {
        while (!st.empty() && a[st.back()] < a[i]) {
            int j = st.back();
            st.pop_back();
            res[j] = i - j;               // 隔了多少天才等到更暖的一天
        }
        st.push_back(i);
    }
    return res;
}

// 柱狀圖最大矩形面積（LeetCode 84）。單調遞增棧 + 高度 0 的哨兵。時間 O(n)，空間 O(n)
long long largestRectangle(const vector<int>& h) {
    int n = static_cast<int>(h.size());
    vector<int> st;                       // 存下標，對應值單調遞增
    long long best = 0;
    for (int i = 0; i <= n; ++i) {
        int cur = (i < n) ? h[i] : 0;     // 哨兵：高度 0 會彈出所有柱子
        while (!st.empty() && h[st.back()] > cur) {  // 嚴格大於才彈：相等高度留在棧裏，避免漏解
            int top = st.back();
            st.pop_back();
            int left = st.empty() ? -1 : st.back();
            long long width = i - left - 1;          // 向右延伸到 i-1，向左到 left+1
            long long area = static_cast<long long>(h[top]) * width;
            if (area > best) best = area;
        }
        st.push_back(i);
    }
    return best;
}

// ---------------- 對照用的 O(n^2) 暴力實現 ----------------

vector<int> nextGreaterIndexBrute(const vector<int>& a) {
    int n = static_cast<int>(a.size());
    vector<int> res(n, -1);
    for (int i = 0; i < n; ++i)
        for (int k = i + 1; k < n; ++k)
            if (a[k] > a[i]) { res[i] = k; break; }
    return res;
}

vector<int> dailyTemperaturesBrute(const vector<int>& a) {
    int n = static_cast<int>(a.size());
    vector<int> res(n, 0);
    for (int i = 0; i < n; ++i)
        for (int k = i + 1; k < n; ++k)
            if (a[k] > a[i]) { res[i] = k - i; break; }
    return res;
}

long long largestRectangleBrute(const vector<int>& h) {
    int n = static_cast<int>(h.size());
    long long best = 0;
    for (int i = 0; i < n; ++i) {
        int left = i, right = i;
        while (left - 1 >= 0 && h[left - 1] >= h[i]) --left;
        while (right + 1 < n && h[right + 1] >= h[i]) ++right;
        long long area = static_cast<long long>(h[i]) * (right - left + 1);
        if (area > best) best = area;
    }
    return best;
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

static void printVec(const vector<int>& v) {
    for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        if (i) cout << " ";
        cout << v[i];
    }
    cout << "\n";
}

int main() {
    int n;
    if (cin >> n) {  // IO 模式
        vector<int> a(n);
        for (int i = 0; i < n; ++i) {
            if (!(cin >> a[i])) a[i] = 0;      // 輸入被截斷時用 0 兜底
        }
        printVec(nextGreaterIndex(a));
        printVec(dailyTemperatures(a));
        cout << largestRectangle(a) << "\n";
        return 0;
    }

    // README 示例：a = {2, 1, 2, 4, 3}
    {
        vector<int> a = {2, 1, 2, 4, 3};
        assert(nextGreaterIndex(a) == vector<int>({3, 2, 3, -1, -1}));
        assert(dailyTemperatures(a) == vector<int>({3, 1, 1, 0, 0}));
        assert(largestRectangle(a) == 6);      // 高 2 寬 3
    }

    // 經典用例
    {
        vector<int> t = {73, 74, 75, 71, 69, 72, 76, 73};
        assert(nextGreaterIndex(t) == vector<int>({1, 2, 6, 5, 5, 6, -1, -1}));
        assert(dailyTemperatures(t) == vector<int>({1, 1, 4, 2, 1, 1, 0, 0}));
    }
    assert(largestRectangle({2, 1, 5, 6, 2, 3}) == 10);   // LeetCode 84 官方用例
    assert(largestRectangle({2, 4}) == 4);
    assert(largestRectangle({1, 1, 1, 1}) == 4);
    assert(largestRectangle({5}) == 5);
    assert(largestRectangle({0}) == 0);
    assert(largestRectangle({0, 0, 0}) == 0);

    // 遞減 / 遞增 / 全相同：三種極端形態
    assert(nextGreaterIndex({5, 4, 3, 2, 1}) == vector<int>({-1, -1, -1, -1, -1}));
    assert(dailyTemperatures({5, 4, 3, 2, 1}) == vector<int>({0, 0, 0, 0, 0}));
    assert(largestRectangle({5, 4, 3, 2, 1}) == 9);       // 高 3 寬 3
    assert(nextGreaterIndex({1, 2, 3, 4, 5}) == vector<int>({1, 2, 3, 4, -1}));
    assert(dailyTemperatures({1, 2, 3, 4, 5}) == vector<int>({1, 1, 1, 1, 0}));
    assert(largestRectangle({1, 2, 3, 4, 5}) == 9);       // 高 3 寬 3
    assert(nextGreaterIndex({3, 3, 3}) == vector<int>({-1, -1, -1}));  // 嚴格大於，相等不算
    assert(dailyTemperatures({3, 3, 3}) == vector<int>({0, 0, 0}));
    assert(largestRectangle({3, 3, 3}) == 9);

    // 空數組
    assert(nextGreaterIndex({}).empty());
    assert(dailyTemperatures({}).empty());
    assert(largestRectangle({}) == 0);

    // 答案自洽：nextGreater 與 dailyTemperatures 必須指向同一個位置
    {
        vector<vector<int>> cases = {{2, 1, 2, 4, 3}, {1}, {4, 2, 9, 1, 7}, {0, 0, 5, 0}};
        for (const auto& arr : cases) {
            vector<int> ng = nextGreaterIndex(arr);
            vector<int> dt = dailyTemperatures(arr);
            for (int i = 0; i < static_cast<int>(arr.size()); ++i) {
                if (ng[i] == -1) {
                    assert(dt[i] == 0);
                } else {
                    assert(dt[i] == ng[i] - i);
                    assert(arr[ng[i]] > arr[i]);
                    for (int k = i + 1; k < ng[i]; ++k) assert(arr[k] <= arr[i]);  // 中間沒有更大的
                }
            }
        }
    }

    LCG rng(20260928ULL);

    // 隨機對拍：三個函數全部與 O(n^2) 暴力解比對
    for (int t = 0; t < 500; ++t) {
        int n2 = rng.next(0, 40);
        vector<int> arr;
        for (int i = 0; i < n2; ++i) arr.push_back(rng.next(0, 12));   // 值域小 → 大量重複
        assert(nextGreaterIndex(arr) == nextGreaterIndexBrute(arr));
        assert(dailyTemperatures(arr) == dailyTemperaturesBrute(arr));
        assert(largestRectangle(arr) == largestRectangleBrute(arr));
        if (n2 > 0) {                                                  // 面積的上下界
            int mx = *max_element(arr.begin(), arr.end());
            assert(mx <= largestRectangle(arr));
            assert(largestRectangle(arr) <= static_cast<long long>(mx) * n2);
        }
    }

    cout << "all tests passed" << endl;
    return 0;
}
