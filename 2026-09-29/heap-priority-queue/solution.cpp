// 二元堆與優先佇列（堆排序 / Top-K / 合併 K 個有序陣列 / 雙堆維護中位數）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：二元堆是用**陣列**存的完全二元樹：下標 i 的父節點是 (i-1)/2，
//   左右孩子是 2i+1、2i+2。完全二元樹高度 floor(log2 n)，所以上浮 / 下沉最多走 log n 層，
//   push / pop 是 O(log n)，top 是 O(1)。堆只維護**偏序**（父節點優於兩個孩子），
//   不維護全序 —— 只取最大/最小值時沒必要把整個集合排好，這是它比有序陣列便宜的原因。
//
//   自底向上建堆（從最後一個非葉節點 (n-2)/2 往前 sift_down）是 O(n) 而非 O(n log n)：
//   大部分節點深度很淺，級數求和後收斂到 2n。
//
//   三個應用：
//     1. 堆排序：建堆後反覆彈出堆頂（原地版把堆頂換到末尾再下沉，額外空間 O(1)，不穩定）；
//     2. Top-K 最小：維護大小不超過 k 的**最大堆**，堆頂就是當前候選裏最大的，
//        新元素比它小就替換 —— O(n log k)，n 大而 k 小時遠快於排序；
//     3. 合併 K 個有序陣列：把每個陣列的當前頭元素放進最小堆，彈出最小的再補上下一個，
//        O(N log K)。多路歸併、定時器與任務調度都是同一個套路。
//   彩蛋 MedianMaintainer：最大堆存較小的一半、最小堆存較大的一半，插入 O(log n)、查詢 O(1)。
//
// 輸入（空白分隔）：n k / a1..an / m / 接著 m 段「長度 L + L 個元素」
// 輸出：堆排序結果、最小的 k 個元素、合併後的結果（均爲升序）
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

using ll = long long;

// Compare(a, b) 爲真表示 a 的優先級更高（更應該靠近堆頂）
template <class Compare>
class BinaryHeap {
public:
    explicit BinaryHeap(Compare cmp = Compare()) : cmp(cmp) {}

    size_t size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    const vector<ll>& raw() const { return data; }
    bool better(ll a, ll b) const { return cmp(a, b); }

    ll top() const { return data.front(); }

    void push(ll x) {                              // O(log n)
        data.push_back(x);
        siftUp((int)data.size() - 1);
    }

    ll pop() {                                     // O(log n)
        ll topVal = data.front();
        ll last = data.back();
        data.pop_back();
        if (!data.empty()) {
            data[0] = last;
            siftDown(0);
        }
        return topVal;
    }

    ll replaceTop(ll x) {                          // 彈出堆頂並插入 x，只下沉一次
        ll topVal = data.front();
        data[0] = x;
        siftDown(0);
        return topVal;
    }

    void heapify(const vector<ll>& values) {       // 自底向上建堆，O(n)
        data = values;
        for (int i = ((int)data.size() - 2) / 2; i >= 0; --i) siftDown(i);
    }

private:
    vector<ll> data;
    Compare cmp;

    void siftUp(int i) {                           // 與父節點比較，更優則上移
        ll x = data[i];
        while (i > 0) {
            int parent = (i - 1) / 2;
            if (!cmp(x, data[parent])) break;
            data[i] = data[parent];
            i = parent;
        }
        data[i] = x;
    }

    void siftDown(int i) {                         // 與更優的那個孩子交換
        int n = (int)data.size();
        ll x = data[i];
        while (true) {
            int left = 2 * i + 1;
            if (left >= n) break;
            int right = left + 1;
            int child = left;
            if (right < n && cmp(data[right], data[left])) child = right;
            if (!cmp(data[child], x)) break;
            data[i] = data[child];
            i = child;
        }
        data[i] = x;
    }
};

using MinHeap = BinaryHeap<less<ll>>;      // 堆頂最小
using MaxHeap = BinaryHeap<greater<ll>>;   // 堆頂最大

template <class Compare>
bool isValidHeap(const BinaryHeap<Compare>& h) {
    const vector<ll>& d = h.raw();
    for (size_t i = 1; i < d.size(); ++i)
        if (h.better(d[i], d[(i - 1) / 2])) return false;
    return true;
}

// ---------------- 應用 1：堆排序 ----------------

vector<ll> heapSort(const vector<ll>& values) {
    MinHeap h;
    h.heapify(values);
    vector<ll> out;
    while (!h.empty()) out.push_back(h.pop());
    return out;
}

// 原地堆排序：建最大堆，把堆頂換到末尾再下沉。額外空間 O(1)，不穩定
vector<ll> heapSortInPlace(const vector<ll>& values) {
    vector<ll> a = values;
    int n = (int)a.size();
    auto siftDown = [&](int i, int size) {
        ll x = a[i];
        while (true) {
            int left = 2 * i + 1;
            if (left >= size) break;
            int right = left + 1;
            int child = left;
            if (right < size && a[right] > a[left]) child = right;
            if (a[child] <= x) break;
            a[i] = a[child];
            i = child;
        }
        a[i] = x;
    };
    for (int i = (n - 2) / 2; i >= 0; --i) siftDown(i, n);   // 1) 建最大堆
    for (int end = n - 1; end > 0; --end) {                   // 2) 堆頂換到末尾
        swap(a[0], a[end]);
        siftDown(0, end);
    }
    return a;
}

// ---------------- 應用 2：Top-K 最小 ----------------

vector<ll> topKSmallest(const vector<ll>& values, int k) {
    if (k <= 0) return {};
    MaxHeap h;                                    // 堆頂是當前 k 個候選裏最大的
    for (ll v : values) {
        if ((int)h.size() < k) h.push(v);
        else if (v < h.top()) h.replaceTop(v);
    }
    vector<ll> out = h.raw();
    sort(out.begin(), out.end());                 // 輸出升序
    return out;
}

// ---------------- 應用 3：合併 K 個有序陣列 ----------------

vector<ll> mergeKSorted(const vector<vector<ll>>& arrays) {
    using Item = tuple<ll, int, int>;             // (值, 陣列編號, 元素下標)
    priority_queue<Item, vector<Item>, greater<Item>> pq;
    for (int i = 0; i < (int)arrays.size(); ++i)
        if (!arrays[i].empty()) pq.emplace(arrays[i][0], i, 0);

    vector<ll> out;
    while (!pq.empty()) {
        auto [val, i, j] = pq.top();
        pq.pop();
        out.push_back(val);
        if (j + 1 < (int)arrays[i].size()) pq.emplace(arrays[i][j + 1], i, j + 1);
    }
    return out;
}

// ---------------- 彩蛋：雙堆維護數據流中位數 ----------------

class MedianMaintainer {
public:
    void add(ll x) {
        if (lo.empty() || x <= lo.top()) lo.push(x);
        else hi.push(x);
        if (lo.size() > hi.size() + 1) {
            hi.push(lo.top());
            lo.pop();
        } else if (hi.size() > lo.size()) {
            lo.push(hi.top());
            hi.pop();
        }
    }

    // 空時返回 NaN；否則奇數個返回中間值，偶數個返回中間兩個的平均
    double median() const {
        if (lo.empty() && hi.empty()) return numeric_limits<double>::quiet_NaN();
        if (lo.size() == hi.size()) return ((double)lo.top() + (double)hi.top()) / 2.0;
        return (double)lo.top();
    }

private:
    priority_queue<ll> lo;                                          // 最大堆：較小的一半
    priority_queue<ll, vector<ll>, greater<ll>> hi;                 // 最小堆：較大的一半
};

// ---------------- IO ----------------

template <typename T>
static string joinInts(const vector<T>& v) {
    ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}

static void runIo(const string& data) {
    istringstream iss(data);
    auto nextInt = [&]() -> ll {
        ll x = 0;
        if (!(iss >> x)) x = 0;                   // 輸入被截斷時用 0 兜底
        return x;
    };
    ll n = nextInt();
    ll k = nextInt();
    vector<ll> a;
    for (ll i = 0; i < n; ++i) a.push_back(nextInt());
    ll m = nextInt();
    vector<vector<ll>> arrays;
    for (ll i = 0; i < m; ++i) {
        ll len = nextInt();
        vector<ll> arr;
        for (ll j = 0; j < len; ++j) arr.push_back(nextInt());
        arrays.push_back(arr);
    }

    cout << joinInts(heapSort(a)) << '\n';                     // n = 0 時輸出空行
    cout << joinInts(topKSmallest(a, (int)k)) << '\n';
    cout << joinInts(mergeKSorted(arrays)) << '\n';
}

static bool nanOrClose(double got, double want) {
    if (isnan(got)) return false;
    return fabs(got - want) < 1e-9;
}

static void runTests() {
    // README 示例：a = [5, 3, 8, 1, 4]，k = 3，兩個有序陣列 [1,4,7] 與 [2,3,9]
    vector<ll> a = {5, 3, 8, 1, 4};
    assert(heapSort(a) == vector<ll>({1, 3, 4, 5, 8}));
    assert(heapSortInPlace(a) == vector<ll>({1, 3, 4, 5, 8}));
    assert(topKSmallest(a, 3) == vector<ll>({1, 3, 4}));
    assert(mergeKSorted({{1, 4, 7}, {2, 3, 9}}) == vector<ll>({1, 2, 3, 4, 7, 9}));

    // 堆的基本性質
    MinHeap h;
    h.heapify({3, 1, 6, 5, 2, 4});
    assert(isValidHeap(h));
    assert(h.top() == 1);
    vector<ll> drained;
    while (!h.empty()) drained.push_back(h.pop());
    assert(drained == vector<ll>({1, 2, 3, 4, 5, 6}));

    MaxHeap hm;
    hm.heapify({3, 1, 6, 5, 2, 4});
    assert(isValidHeap(hm));
    assert(hm.top() == 6);
    vector<ll> drainedMax;
    while (!hm.empty()) drainedMax.push_back(hm.pop());
    assert(drainedMax == vector<ll>({6, 5, 4, 3, 2, 1}));

    // 邊界：空堆 / 單元素 / 重複值 / 負數
    MinHeap empty;
    assert(empty.size() == 0);
    assert(heapSort({}).empty());
    assert(heapSortInPlace({}).empty());
    assert(heapSort({7}) == vector<ll>({7}));
    assert(heapSort({2, 2, 2}) == vector<ll>({2, 2, 2}));
    assert(heapSort({-3, 0, -1, 5}) == vector<ll>({-3, -1, 0, 5}));
    assert(heapSort({1, 2, 3, 4, 5}) == vector<ll>({1, 2, 3, 4, 5}));
    assert(heapSort({5, 4, 3, 2, 1}) == vector<ll>({1, 2, 3, 4, 5}));

    // Top-K 邊界
    assert(topKSmallest(a, 0).empty());
    assert(topKSmallest(a, -2).empty());                       // k 爲負：空
    assert(topKSmallest(a, 100) == vector<ll>({1, 3, 4, 5, 8}));
    assert(topKSmallest({}, 3).empty());
    assert(topKSmallest({4, 4, 1, 4}, 2) == vector<ll>({1, 4}));

    // 合併：空陣列 / 全空 / 單個陣列 / 值全相等
    assert(mergeKSorted({}).empty());
    assert(mergeKSorted({{}, {}}).empty());
    assert(mergeKSorted({{1, 2, 3}}) == vector<ll>({1, 2, 3}));
    assert(mergeKSorted({{}, {1, 2}, {}}) == vector<ll>({1, 2}));
    assert(mergeKSorted({{1, 1}, {1, 1}}) == vector<ll>({1, 1, 1, 1}));

    // 中位數維護器
    MedianMaintainer mm;
    assert(isnan(mm.median()));
    for (ll x : {5LL, 2LL, 9LL, 1LL, 7LL}) mm.add(x);
    assert(nanOrClose(mm.median(), 5.0));                      // [1,2,5,7,9] → 5
    MedianMaintainer mm2;
    for (ll x : {1LL, 2LL, 3LL, 4LL}) mm2.add(x);
    assert(nanOrClose(mm2.median(), 2.5));                     // 偶數個取中間兩個的平均

    mt19937 rng(20260929);

    // 隨機對拍：全部與 sort() 的結果比對
    for (int iter = 0; iter < 400; ++iter) {
        int n = (int)(rng() % 41);
        vector<ll> vals;
        for (int i = 0; i < n; ++i) vals.push_back((ll)(rng() % 21) - 10);
        vector<ll> sortedVals = vals;
        sort(sortedVals.begin(), sortedVals.end());
        assert(heapSort(vals) == sortedVals);
        assert(heapSortInPlace(vals) == sortedVals);

        int k = -1 + (int)(rng() % (n + 4));
        int take = min(max(k, 0), n);                          // k 可能爲負或超過 n
        vector<ll> wantK(sortedVals.begin(), sortedVals.begin() + take);
        assert(topKSmallest(vals, k) == wantK);

        MinHeap hr;
        hr.heapify(vals);
        assert(isValidHeap(hr));
        vector<ll> drainedR;
        while (!hr.empty()) drainedR.push_back(hr.pop());
        assert(drainedR == sortedVals);

        // 合併隨機個有序陣列
        int m = (int)(rng() % 6);
        vector<vector<ll>> arrays;
        vector<ll> flat;
        for (int i = 0; i < m; ++i) {
            int size = (int)(rng() % 7);
            vector<ll> arr;
            for (int j = 0; j < size; ++j) arr.push_back((ll)(rng() % 21) - 10);
            sort(arr.begin(), arr.end());
            arrays.push_back(arr);
            for (ll v : arr) flat.push_back(v);
        }
        sort(flat.begin(), flat.end());
        assert(mergeKSorted(arrays) == flat);

        // 中位數維護器與「排序後取中位」逐前綴比對
        MedianMaintainer mmr;
        vector<ll> seen;
        for (ll x : vals) {
            mmr.add(x);
            seen.push_back(x);
            sort(seen.begin(), seen.end());
            int sz = (int)seen.size();
            double want = (sz % 2 == 1)
                              ? (double)seen[sz / 2]
                              : ((double)seen[sz / 2 - 1] + (double)seen[sz / 2]) / 2.0;
            assert(nanOrClose(mmr.median(), want));
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
