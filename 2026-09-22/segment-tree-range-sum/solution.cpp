// 線段樹（區間和 + 懶標記區間加）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 題意：給定數組 a，支持（0-based，閉區間 [l, r]）：
//     query(l, r) 求 a[l..r] 的和
//     add(l, r, v) 把 a[l..r] 每個元素加 v
//   暴力單次 O(n)，m 次操作 O(nm)；線段樹把兩種操作都降到 O(log n)。
//
// 思路：線段樹是描述「區間」的二叉樹——葉子存單個元素，內部節點存左右兒子區間的合併和。
//   區間查詢：把目標區間拆成 O(log n) 個節點，完全覆蓋的節點直接返回。
//   區間修改：整段被覆蓋時只更新該節點 sum 並打懶標記 lazy，不再下遞歸；
//             下次必須深入時再 push_down 把標記分給兩個兒子。
//   懶標記的本質是「延遲執行」：沒人問細節就不必真的改到葉子。
//   堆式存儲：節點 p 的左兒子 2p、右兒子 2p+1，空間開 4n。
//
// 輸入：第一行 n q；第二行 n 個整數；接下來 q 行：1 l r（查詢）/ 2 l r v（區間加）
// 輸出：每個查詢操作輸出一行區間和
// 無 stdin 輸入時運行內置斷言測試。
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;
using ll = long long;

class SegmentTree {
public:
    explicit SegmentTree(const vector<ll>& arr) : n(static_cast<int>(arr.size())) {
        sum.assign(4 * n + 5, 0);   // sum[p]：節點 p 對應區間的元素和
        lazy.assign(4 * n + 5, 0);  // lazy[p]：待下推給子孫的「整體加多少」
        if (n > 0) build(1, 0, n - 1, arr);
    }

    int size() const { return n; }

    // 區間加：a[ql..qr] += v，O(log n)
    void add(int ql, int qr, ll v) {
        if (n == 0) return;
        add(1, 0, n - 1, ql, qr, v);
    }

    // 區間求和：返回 a[ql..qr] 的和，O(log n)
    ll query(int ql, int qr) {
        if (n == 0) return 0;
        return query(1, 0, n - 1, ql, qr);
    }

private:
    int n;
    vector<ll> sum, lazy;

    // 自底向上建樹，O(n)
    void build(int p, int l, int r, const vector<ll>& arr) {
        if (l == r) {
            sum[p] = arr[l];
            return;
        }
        int mid = (l + r) / 2;
        build(2 * p, l, mid, arr);
        build(2 * p + 1, mid + 1, r, arr);
        sum[p] = sum[2 * p] + sum[2 * p + 1];
    }

    // 給節點 p 整段加 v：更新 sum，並累積懶標記
    void apply(int p, int l, int r, ll v) {
        sum[p] += v * (r - l + 1);
        lazy[p] += v;
    }

    // 把 p 的懶標記下推給兩個兒子，只在需要深入時使用
    void pushDown(int p, int l, int r) {
        if (lazy[p] == 0 || l == r) return;
        int mid = (l + r) / 2;
        ll v = lazy[p];
        apply(2 * p, l, mid, v);
        apply(2 * p + 1, mid + 1, r, v);
        lazy[p] = 0;
    }

    void add(int p, int l, int r, int ql, int qr, ll v) {
        if (ql <= l && r <= qr) {  // 完全覆蓋，打標記後直接返回
            apply(p, l, r, v);
            return;
        }
        pushDown(p, l, r);
        int mid = (l + r) / 2;
        if (ql <= mid) add(2 * p, l, mid, ql, qr, v);
        if (qr > mid) add(2 * p + 1, mid + 1, r, ql, qr, v);
        sum[p] = sum[2 * p] + sum[2 * p + 1];
    }

    ll query(int p, int l, int r, int ql, int qr) {
        if (ql <= l && r <= qr) return sum[p];  // 完全覆蓋，直接返回整段和
        pushDown(p, l, r);
        int mid = (l + r) / 2;
        ll ans = 0;
        if (ql <= mid) ans += query(2 * p, l, mid, ql, qr);
        if (qr > mid) ans += query(2 * p + 1, mid + 1, r, ql, qr);
        return ans;
    }
};

int main() {
    int n, q;
    if (cin >> n >> q) {  // IO 模式
        vector<ll> arr(n);
        for (int i = 0; i < n; ++i) cin >> arr[i];
        SegmentTree st(arr);
        for (int i = 0; i < q; ++i) {
            int kind, l, r;
            cin >> kind >> l >> r;
            if (kind == 1) {
                cout << st.query(l, r) << "\n";
            } else {
                ll v;
                cin >> v;
                st.add(l, r, v);
            }
        }
        return 0;
    }

    {
        vector<ll> arr = {1, 3, 5, 7, 9, 11};
        SegmentTree st(arr);
        assert(st.query(0, 5) == 36);
        assert(st.query(1, 3) == 15);  // 3 + 5 + 7
        assert(st.query(2, 2) == 5);   // 單點查詢
        st.add(1, 3, 10);              // [1, 13, 15, 17, 9, 11]
        assert(st.query(1, 3) == 45);
        assert(st.query(0, 5) == 66);  // 36 + 3 * 10
        assert(st.query(0, 0) == 1);   // 區間外的點不受影響
        assert(st.query(4, 5) == 20);
    }

    // 懶標記疊加與部分覆蓋交叉驗證
    {
        SegmentTree st(vector<ll>(5, 0));
        // 結果數組爲 [2, 5, 5, 5, 2]
        st.add(0, 4, 2);
        st.add(1, 3, 3);
        assert(st.query(0, 4) == 19);
        assert(st.query(0, 0) == 2);
        assert(st.query(1, 3) == 15);  // 5 + 5 + 5
        assert(st.query(4, 4) == 2);
    }

    // 單元素與空數組邊界
    {
        SegmentTree st(vector<ll>{42});
        assert(st.query(0, 0) == 42);
        st.add(0, 0, -40);
        assert(st.query(0, 0) == 2);

        SegmentTree st0(vector<ll>{});
        assert(st0.query(0, 0) == 0);
    }

    cout << "all tests passed" << endl;
    return 0;
}
