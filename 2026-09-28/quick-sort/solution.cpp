// 快速排序與快速選擇（隨機化 pivot + 三路分區 + Lomuto 分區）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：快排是「分治 + 原地分區」。每輪挑一個 pivot 把區間切成
//   < pivot / == pivot / > pivot 三段，只對左右兩段遞歸。
//   1) Lomuto 分區：寫法最直觀，但重複元素多時兩段嚴重不平衡，全相同退化 O(n^2)；
//   2) 三路分區（荷蘭國旗）：相等元素一次整段歸位，全相同元素時是 O(n)；
//   3) 隨機化 pivot：避免已排序輸入退化，期望 O(n log n)；
//   4) 先遞歸短半邊、循環處理長半邊 → 棧深度 O(log n)；
//   5) Quickselect：只遞歸包含 k 的那一段，平均 O(n) 拿到第 k 小。
//
// 輸入（空白分隔）：n / a1..an（n=0 時省略）/ k
// 輸出：第 1 行升序排序結果（n=0 時空行）；第 2 行第 k 小（1-based，越界輸出 -1）
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

const int SMALL = 16;  // 小區間閾值：低於這個長度改用插入排序

// 全局隨機數引擎，固定種子保證可復現
static mt19937 rng_engine(20260928);

static int randInt(int lo, int hi) {
    uniform_int_distribution<int> dist(lo, hi);
    return dist(rng_engine);
}

// 對 a[lo..hi] 做插入排序。小區間裡常數比快排更小
void insertionSort(vector<int>& a, int lo, int hi) {
    for (int i = lo + 1; i <= hi; ++i) {
        int x = a[i];
        int j = i - 1;
        while (j >= lo && a[j] > x) {
            a[j + 1] = a[j];
            --j;
        }
        a[j + 1] = x;
    }
}

// Lomuto 分區：以 a[hi] 爲 pivot，把 < pivot 的換到左邊，返回 pivot 最終下標
int lomutoPartition(vector<int>& a, int lo, int hi) {
    int pivot = a[hi];
    int i = lo;
    for (int j = lo; j < hi; ++j) {
        if (a[j] < pivot) {
            swap(a[i], a[j]);
            ++i;
        }
    }
    swap(a[i], a[hi]);
    return i;
}

// 三路分區：a[lo..lt-1] < pivot，a[lt..gt] == pivot，a[gt+1..hi] > pivot
pair<int, int> partition3(vector<int>& a, int lo, int hi) {
    int pivot = a[randInt(lo, hi)];
    int lt = lo, i = lo, gt = hi;
    while (i <= gt) {
        if (a[i] < pivot) {
            swap(a[lt], a[i]);
            ++lt;
            ++i;
        } else if (a[i] > pivot) {
            swap(a[i], a[gt]);
            --gt;              // 換過來的元素還沒看過，i 不動
        } else {
            ++i;
        }
    }
    return {lt, gt};
}

// 三路快排主過程：先遞歸短半邊，長半邊用循環，棧深度 O(log n)
void quickSort3(vector<int>& a, int lo, int hi) {
    while (lo < hi) {
        if (hi - lo + 1 <= SMALL) {
            insertionSort(a, lo, hi);
            return;
        }
        auto [lt, gt] = partition3(a, lo, hi);
        if (lt - lo < hi - gt) {          // 左段更短
            quickSort3(a, lo, lt - 1);
            lo = gt + 1;
        } else {                          // 右段更短
            quickSort3(a, gt + 1, hi);
            hi = lt - 1;
        }
    }
}

// 三路快排（原地）。平均 O(n log n)，最壞 O(n^2)，棧空間 O(log n)
void quickSort(vector<int>& a) { quickSort3(a, 0, static_cast<int>(a.size()) - 1); }

// Lomuto 版快排（原地），作爲對照實現保留
void quickSortLomuto(vector<int>& a, int lo, int hi) {
    while (lo < hi) {
        if (hi - lo + 1 <= SMALL) {
            insertionSort(a, lo, hi);
            return;
        }
        int p = lomutoPartition(a, lo, hi);
        if (p - lo < hi - p) {
            quickSortLomuto(a, lo, p - 1);
            lo = p + 1;
        } else {
            quickSortLomuto(a, p + 1, hi);
            hi = p - 1;
        }
    }
}

void quickSortLomuto(vector<int>& a) {
    if (a.empty()) return;
    quickSortLomuto(a, 0, static_cast<int>(a.size()) - 1);
}

// 第 k 小元素（k 爲 0-based），原地修改 a。平均 O(n)
int quickSelect(vector<int>& a, int k) {
    int lo = 0, hi = static_cast<int>(a.size()) - 1;
    while (lo <= hi) {
        if (hi - lo + 1 <= SMALL) {
            insertionSort(a, lo, hi);
            return a[k];
        }
        auto [lt, gt] = partition3(a, lo, hi);
        if (k < lt) hi = lt - 1;          // k 落在「小於」段
        else if (k > gt) lo = gt + 1;     // k 落在「大於」段
        else return a[k];                 // k 落在「等於」段，pivot 就是答案
    }
    return -1;                            // k 越界，理論上不會走到
}

// 第 k 小元素（k 爲 1-based）；k 越界返回 -1。內部拷貝以免破壞原數組
int kthSmallest(const vector<int>& a, int k) {
    if (k < 1 || k > static_cast<int>(a.size())) return -1;
    vector<int> b = a;
    return quickSelect(b, k - 1);
}

int main() {
    int n;
    if (cin >> n) {  // IO 模式
        vector<int> a(n);
        for (int i = 0; i < n; ++i) {
            if (!(cin >> a[i])) a[i] = 0;      // 輸入被截斷時用 0 兜底
        }
        int k = 0;
        if (!(cin >> k)) k = 0;

        quickSort(a);
        for (int i = 0; i < n; ++i) {          // n = 0 時輸出空行
            if (i) cout << " ";
            cout << a[i];
        }
        cout << "\n";
        cout << (a.empty() ? -1 : kthSmallest(a, k)) << "\n";
        return 0;
    }

    // README 示例
    {
        vector<int> a = {5, 3, 8, 3, 1, 9, 3};
        vector<int> b = a, c = a;
        quickSort(b);
        quickSortLomuto(c);
        vector<int> exp = {1, 3, 3, 3, 5, 8, 9};
        assert(b == exp);
        assert(c == exp);
        assert(kthSmallest(a, 1) == 1);
        assert(kthSmallest(a, 3) == 3);
        assert(kthSmallest(a, 7) == 9);
        assert(kthSmallest(a, 0) == -1);       // k 越界
        assert(kthSmallest(a, 8) == -1);
    }

    // 空數組 / 單元素
    {
        vector<int> e;
        quickSort(e);
        quickSortLomuto(e);
        assert(e.empty());
        assert(kthSmallest(e, 1) == -1);
        vector<int> one = {42};
        quickSort(one);
        assert(one[0] == 42);
        assert(kthSmallest(one, 1) == 42);
    }

    // 已排序 / 逆序 / 全相同 —— 三種最容易觸發退化的輸入
    {
        vector<vector<int>> cases = {
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
            {10, 9, 8, 7, 6, 5, 4, 3, 2, 1},
            vector<int>(20, 7),
            {-3, -1, -2, -5, -4},
            {0, 0, -1, 1, 0},
        };
        for (const auto& arr : cases) {
            vector<int> exp = arr;
            sort(exp.begin(), exp.end());
            vector<int> b = arr, c = arr;
            quickSort(b);
            quickSortLomuto(c);
            assert(b == exp);
            assert(c == exp);
            for (int k = 1; k <= static_cast<int>(arr.size()); ++k) {
                assert(kthSmallest(arr, k) == exp[k - 1]);
            }
        }
    }

    // Lomuto 分區的返回值必須自洽
    {
        vector<int> b = {4, 2, 7, 2, 9, 1};
        int p = lomutoPartition(b, 0, static_cast<int>(b.size()) - 1);
        for (int i = 0; i < p; ++i) assert(b[i] < b[p]);
        for (int i = p + 1; i < static_cast<int>(b.size()); ++i) assert(b[i] >= b[p]);
    }

    // 三路分區：切出來的三段必須真的有序
    {
        vector<int> c = {5, 1, 5, 3, 5, 2, 5};
        auto [lt, gt] = partition3(c, 0, static_cast<int>(c.size()) - 1);
        int pivotVal = c[lt];
        for (int i = 0; i < lt; ++i) assert(c[i] < pivotVal);
        for (int i = lt; i <= gt; ++i) assert(c[i] == pivotVal);
        for (int i = gt + 1; i < static_cast<int>(c.size()); ++i) assert(c[i] > pivotVal);
    }

    // 小區間閾值附近（<= 16 走插入排序分支）
    for (int n2 = 0; n2 < 40; ++n2) {
        vector<int> arr;
        for (int i = 0; i < n2; ++i) arr.push_back(randInt(-5, 5));
        vector<int> exp = arr;
        sort(exp.begin(), exp.end());
        vector<int> b = arr, c = arr;
        quickSort(b);
        quickSortLomuto(c);
        assert(b == exp);
        assert(c == exp);
    }

    // 隨機對拍：排序結果與 std::sort 比對，quickselect 與排序結果比對
    for (int t = 0; t < 500; ++t) {
        int n2 = randInt(0, 60);
        vector<int> arr;
        for (int i = 0; i < n2; ++i) arr.push_back(randInt(-9, 9));  // 值域小，製造大量重複
        vector<int> exp = arr;
        sort(exp.begin(), exp.end());
        vector<int> b = arr, c = arr;
        quickSort(b);
        quickSortLomuto(c);
        assert(b == exp);
        assert(c == exp);
        for (int k = 1; k <= n2; ++k) assert(kthSmallest(arr, k) == exp[k - 1]);
        assert(kthSmallest(arr, 0) == -1);
        assert(kthSmallest(arr, n2 + 1) == -1);
    }

    // 大數組：確認沒有爆棧
    {
        vector<int> big;
        for (int i = 0; i < 20000; ++i) big.push_back(randInt(-1000, 1000));
        vector<int> exp = big;
        sort(exp.begin(), exp.end());
        vector<int> b = big;
        quickSort(b);
        assert(b == exp);
        assert(kthSmallest(big, 12345) == exp[12344]);
    }

    cout << "all tests passed" << endl;
    return 0;
}
