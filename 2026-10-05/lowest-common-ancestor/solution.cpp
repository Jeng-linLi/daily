// LCA（最近公共祖先）：倍增 / 歐拉序 RMQ，並附帶樹上距離與 k 級祖先
//
// 題意：
//     給定 n 個節點的無向樹（邊表，0-indexed），以 0 號點為根，在線回答：
//       1. lca(u, v)；2. dist(u, v)（樹上距離 / 邊數）；
//       3. kth_ancestor(v, k)（超過根為 -1）；4. is_ancestor(u, v)（含 u == v）。
//     預處理 O(n log n)、查詢 O(log n)（倍增）；另有歐拉序 + Sparse Table 的 O(1) 查詢版
//     與 O(n) 爬 parent 的樸素版作交叉驗證。
//
// 思路：
//     ### 倍增（binary lifting）
//     up[k][v] = v 的第 2^k 級祖先（不存在為 -1），遞推 up[k][v] = up[k-1][ up[k-1][v] ]。
//     利用二進制拆分：
//       - 先把較深的點往上拉 diff = depth[u] - depth[v] 步（對 diff 的每個 1 位跳一次）；
//       - 同深度後兩點一起往上跳：從大到小枚舉 k，若 up[k][u] != up[k][v] 就同時跳
//         （相等代表跳上去會跳過頭）。最後 up[0][u] 即 LCA。
//     ### 歐拉序 + RMQ
//     DFS 進入節點記一次、每訪問完一個子節點返回父節點時再記一次，得到長度 2n-1 的序列。
//     u 與 v 的 LCA = 序列中 [first[u], first[v]] 區間內 depth 最小的節點。
//     區間最小值用 Sparse Table 預處理 → O(1) 查詢。這條路線不需要對齊深度，
//     與倍增結果必然一致，正好互相對拍。
//     ### 樹上距離
//     dist(u, v) = depth[u] + depth[v] - 2 * depth[lca(u, v)]。
//     ### 祖先判定
//     DFS 進出時間 tin / tout：u 是 v 的祖先 ⟺ tin[u] <= tin[v] 且 tout[v] <= tout[u]。
//
// 輸入格式（stdin，全部以空白分隔）：
//     n
//     n-1 行：u v        （無向邊；n == 1 時沒有邊）
//     q
//     q 行：lca u v  |  anc v k
// 輸出格式（stdout）：
//     每個查詢一行：
//       `lca u v` -> `<lca> <dist>`
//       `anc v k` -> `<ancestor>`（不存在為 -1）
//     節點越界或不在根 0 所在分量內時輸出 `-1 -1` / `-1`。
// 輸入被截斷時缺的部分按 0 個查詢處理。無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class LCA {
public:
    int n, root, LOG;
    vector<vector<int>> adj;
    vector<int> depth, tin, tout;
    vector<vector<int>> up;              // up[k][v]
    vector<int> euler, eulerDepth;
    vector<vector<int>> st;              // Sparse Table：存歐拉序下標
    vector<int> logTbl;
    int stK;

    LCA(int n_, const vector<pair<int, int>>& edges, int root_ = 0)
        : n(n_), root(root_), adj((size_t)max(n_, 0)), depth((size_t)max(n_, 0), -1),
          tin((size_t)max(n_, 0), -1), tout((size_t)max(n_, 0), -1) {
        for (const auto& e : edges) {
            int u = e.first, v = e.second;
            if (0 <= u && u < n && 0 <= v && v < n && u != v) {
                adj[(size_t)u].push_back(v);
                adj[(size_t)v].push_back(u);
            }
        }
        for (auto& lst : adj) sort(lst.begin(), lst.end());   // 排序保證 DFS 順序確定

        LOG = 1;
        while ((1 << LOG) <= max(n, 1)) LOG++;
        up.assign((size_t)LOG, vector<int>((size_t)max(n, 0), -1));

        buildTree();
        buildLifting();
        buildEuler();
        buildSparse();
    }

private:
    // ---------------------------------------------------------- 建樹（迭代 DFS）
    void buildTree() {
        if (n == 0 || !(0 <= root && root < n)) return;
        vector<pair<int, int>> stk;
        depth[(size_t)root] = 0;
        up[0][(size_t)root] = -1;
        stk.push_back({root, -1});
        while (!stk.empty()) {
            pair<int, int> cur = stk.back();
            stk.pop_back();
            int u = cur.first, p = cur.second;
            const vector<int>& nb = adj[(size_t)u];
            for (int i = (int)nb.size() - 1; i >= 0; --i) {   // 反序壓棧 → 彈出時為升序
                int v = nb[(size_t)i];
                if (v == p || depth[(size_t)v] >= 0) continue;
                depth[(size_t)v] = depth[(size_t)u] + 1;
                up[0][(size_t)v] = u;
                stk.push_back({v, u});
            }
        }
    }

    // ---------------------------------------------------------- 倍增表
    void buildLifting() {
        for (int k = 1; k < LOG; ++k) {
            const vector<int>& prev = up[(size_t)(k - 1)];
            vector<int>& cur = up[(size_t)k];
            for (int v = 0; v < n; ++v) {
                int p = prev[(size_t)v];
                cur[(size_t)v] = (p < 0) ? -1 : prev[(size_t)p];
            }
        }
    }

    // ---------------------------------------------------------- 歐拉序（迭代 DFS）
    void buildEuler() {
        if (n == 0 || depth[(size_t)root] < 0) return;
        // 注意：push_back 可能讓 vector 重新配置，這裡**不持有** stk.back() 的引用，
        // 一律用 stk.back()[...] 當場存取。
        vector<array<int, 3>> stk;                            // {節點, 父, 下一個孩子下標}
        stk.push_back({root, -1, 0});
        while (!stk.empty()) {
            int u = stk.back()[0], p = stk.back()[1], ci = stk.back()[2];
            if (ci == 0) {                                    // 第一次進入 u
                tin[(size_t)u] = (int)euler.size();
                euler.push_back(u);
                eulerDepth.push_back(depth[(size_t)u]);
            }
            if (ci < (int)adj[(size_t)u].size()) {
                stk.back()[2] = ci + 1;
                int v = adj[(size_t)u][(size_t)ci];
                if (v != p) stk.push_back({v, u, 0});
            } else {
                stk.pop_back();
                tout[(size_t)u] = (int)euler.size() - 1;
                if (p >= 0) {                                 // 從子節點返回，再記一次父節點
                    euler.push_back(p);
                    eulerDepth.push_back(depth[(size_t)p]);
                }
            }
        }
    }

    // ---------------------------------------------------------- Sparse Table
    void buildSparse() {
        int m = (int)eulerDepth.size();
        stK = 1;
        while ((1 << stK) <= max(m, 1)) stK++;
        logTbl.assign((size_t)(m + 1), 0);
        for (int i = 2; i <= m; ++i) logTbl[(size_t)i] = logTbl[(size_t)(i / 2)] + 1;
        st.clear();
        vector<int> base((size_t)m);
        for (int i = 0; i < m; ++i) base[(size_t)i] = i;
        st.push_back(base);
        for (int k = 1; k < stK; ++k) {
            const vector<int>& prev = st[(size_t)(k - 1)];
            int half = 1 << (k - 1);
            vector<int> nxt;
            for (int i = 0; i + (1 << k) <= m; ++i) {
                int a = prev[(size_t)i], b = prev[(size_t)(i + half)];
                // 深度小者勝；深度相同取歐拉序下標小者（保證確定）
                nxt.push_back(eulerDepth[(size_t)a] <= eulerDepth[(size_t)b] ? a : b);
            }
            st.push_back(nxt);
        }
    }

    int rmq(int l, int r) const {
        int k = logTbl[(size_t)(r - l + 1)];
        int a = st[(size_t)k][(size_t)l];
        int b = st[(size_t)k][(size_t)(r - (1 << k) + 1)];
        return eulerDepth[(size_t)a] <= eulerDepth[(size_t)b] ? a : b;
    }

public:
    bool valid(int v) const { return 0 <= v && v < n && depth[(size_t)v] >= 0; }

    // ---------------------------------------------------------- 倍增版 LCA，O(log n)
    int lca(int u, int v) const {
        if (!valid(u) || !valid(v)) return -1;
        if (depth[(size_t)u] < depth[(size_t)v]) swap(u, v);
        int diff = depth[(size_t)u] - depth[(size_t)v];
        int k = 0;
        while (diff) {
            if (diff & 1) {
                u = up[(size_t)k][(size_t)u];
                if (u < 0) return -1;
            }
            diff >>= 1;
            k++;
        }
        if (u == v) return u;
        for (int kk = LOG - 1; kk >= 0; --kk) {
            int pu = up[(size_t)kk][(size_t)u], pv = up[(size_t)kk][(size_t)v];
            if (pu != pv) {                                   // 相等代表跳上去會跳過頭
                if (pu >= 0) u = pu;
                if (pv >= 0) v = pv;
            }
        }
        return up[0][(size_t)u];
    }

    // ---------------------------------------------------------- 歐拉序 + RMQ 版，O(1)
    int lcaRmq(int u, int v) const {
        if (!valid(u) || !valid(v)) return -1;
        int l = tin[(size_t)u], r = tin[(size_t)v];
        if (l > r) swap(l, r);
        return euler[(size_t)rmq(l, r)];
    }

    // ---------------------------------------------------------- 樸素基準，O(n)
    int lcaNaive(int u, int v) const {
        if (!valid(u) || !valid(v)) return -1;
        int a = u, b = v;
        while (depth[(size_t)a] > depth[(size_t)b]) a = up[0][(size_t)a];
        while (depth[(size_t)b] > depth[(size_t)a]) b = up[0][(size_t)b];
        while (a != b) {
            a = up[0][(size_t)a];
            b = up[0][(size_t)b];
        }
        return a;
    }

    int dist(int u, int v) const {
        int w = lca(u, v);
        if (w < 0) return -1;
        return depth[(size_t)u] + depth[(size_t)v] - 2 * depth[(size_t)w];
    }

    int kthAncestor(int v, int k) const {
        if (!valid(v) || k < 0) return -1;
        if (k > depth[(size_t)v]) return -1;
        int i = 0;
        while (k) {
            if (k & 1) {
                v = up[(size_t)i][(size_t)v];
                if (v < 0) return -1;
            }
            k >>= 1;
            i++;
        }
        return v;
    }

    int kthAncestorNaive(int v, int k) const {
        if (!valid(v) || k < 0) return -1;
        if (k > depth[(size_t)v]) return -1;
        int cur = v;
        for (int t = 0; t < k; ++t) cur = up[0][(size_t)cur];
        return cur;
    }

    bool isAncestor(int u, int v) const {
        if (!valid(u) || !valid(v)) return false;
        return tin[(size_t)u] <= tin[(size_t)v] && tout[(size_t)v] <= tout[(size_t)u];
    }
};

// ---------------------------------------------------------------- 小工具
// 只接受 [+-]?digits，與 Python 版的 parse_int 完全一致；其它 token 一律當成「輸入結束」，
// 避免一個語言當場崩潰、另一個若無其事地繼續解析。
static bool tryLL(const string& s, long long& out) {
    if (s.empty()) return false;
    size_t i = 0;
    bool neg = false;
    if (s[0] == '-' || s[0] == '+') {
        neg = (s[0] == '-');
        i = 1;
    }
    if (i >= s.size()) return false;
    long long val = 0;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') return false;
        val = val * 10 + (s[i] - '0');
        if (val > 4000000000000000000LL) return false;      // 溢出保護
    }
    out = neg ? -val : val;
    return true;
}

static mt19937 rngEngine;
static int rndInt(int lo, int hi) { return lo + (int)(rngEngine() % (unsigned)(hi - lo + 1)); }

static vector<pair<int, int>> buildRandomTree(int n) {
    vector<pair<int, int>> edges;
    for (int v = 1; v < n; ++v) edges.push_back({rndInt(0, v - 1), v});
    return edges;
}

// ---------------------------------------------------------------- IO 模式
static void runIo(const vector<string>& toks) {
    size_t pos = 0;
    auto nxt = [&]() -> const string* {
        if (pos < toks.size()) return &toks[pos++];
        return nullptr;
    };
    // 取得下一個 token 並轉成整數；缺 token 或非數字都返回 nullopt
    auto nextInt = [&]() -> optional<long long> {
        const string* v = nxt();
        if (!v) return nullopt;
        long long out = 0;
        if (!tryLL(*v, out)) return nullopt;
        return out;
    };
    auto nextIntOr = [&](int defVal) -> int {
        optional<long long> v = nextInt();
        return v ? (int)*v : defVal;
    };

    int n = max(nextIntOr(0), 0);
    vector<pair<int, int>> edges;
    if (n > 1) {
        for (int i = 0; i < n - 1; ++i) {
            optional<long long> a = nextInt(), b = nextInt();
            if (!a || !b) break;                              // 輸入截斷或出現非數字 token
            edges.push_back({(int)*a, (int)*b});
        }
    }
    int q = nextIntOr(0);

    LCA solver(n, edges, 0);
    for (int i = 0; i < q; ++i) {
        const string* op = nxt();
        if (!op) break;
        if (*op == "lca") {
            optional<long long> u = nextInt(), v = nextInt();
            if (!u || !v) break;
            cout << solver.lca((int)*u, (int)*v) << ' ' << solver.dist((int)*u, (int)*v) << '\n';
        } else if (*op == "anc") {
            optional<long long> v = nextInt(), k = nextInt();
            if (!v || !k) break;
            cout << solver.kthAncestor((int)*v, (int)*k) << '\n';
        }
    }
}

// ---------------------------------------------------------------- 測試
static void runTests() {
    // 單點樹
    {
        LCA s(1, {}, 0);
        assert((s.depth == vector<int>{0}));
        assert(s.lca(0, 0) == 0);
        assert(s.dist(0, 0) == 0);
        assert(s.kthAncestor(0, 0) == 0);
        assert(s.kthAncestor(0, 1) == -1);
        assert(s.isAncestor(0, 0));
        assert(s.lcaRmq(0, 0) == 0);
        assert(s.lcaNaive(0, 0) == 0);
    }

    // README 示例樹：0-1, 0-2, 1-3, 1-4, 2-5, 2-6
    {
        vector<pair<int, int>> edges = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}};
        LCA s(7, edges, 0);
        assert((s.depth == vector<int>{0, 1, 1, 2, 2, 2, 2}));
        assert(s.lca(3, 4) == 1);
        assert(s.dist(3, 4) == 2);
        assert(s.lca(3, 5) == 0);
        assert(s.dist(3, 5) == 4);
        assert(s.lca(5, 6) == 2);
        assert(s.dist(5, 6) == 2);
        assert(s.lca(0, 3) == 0);
        assert(s.dist(0, 3) == 2);
        assert(s.lca(3, 3) == 3);
        assert(s.dist(3, 3) == 0);
        assert(s.kthAncestor(3, 0) == 3);
        assert(s.kthAncestor(3, 1) == 1);
        assert(s.kthAncestor(3, 2) == 0);
        assert(s.kthAncestor(3, 3) == -1);
        assert(s.isAncestor(0, 6));
        assert(!s.isAncestor(1, 6));
        assert(s.isAncestor(2, 2));
        for (int u = 0; u < 7; ++u)
            for (int v = 0; v < 7; ++v)
                assert(s.lca(u, v) == s.lcaRmq(u, v) && s.lca(u, v) == s.lcaNaive(u, v));
    }

    // 鏈狀樹
    {
        vector<pair<int, int>> chain;
        for (int i = 0; i < 9; ++i) chain.push_back({i, i + 1});
        LCA s(10, chain, 0);
        assert((s.depth == vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));
        assert(s.lca(0, 9) == 0);
        assert(s.dist(0, 9) == 9);
        assert(s.dist(3, 7) == 4);
        assert(s.kthAncestor(9, 9) == 0);
        assert(s.kthAncestor(9, 10) == -1);
        for (int u = 0; u < 10; ++u)
            for (int v = 0; v < 10; ++v) {
                assert(s.lca(u, v) == min(u, v));
                assert(s.dist(u, v) == abs(u - v));
                assert(s.lca(u, v) == s.lcaRmq(u, v) && s.lca(u, v) == s.lcaNaive(u, v));
            }
    }

    // 星狀樹
    {
        vector<pair<int, int>> star;
        for (int i = 1; i < 8; ++i) star.push_back({0, i});
        LCA s(8, star, 0);
        for (int u = 1; u < 8; ++u)
            for (int v = 1; v < 8; ++v) {
                assert(s.lca(u, v) == (u != v ? 0 : u));
                assert(s.dist(u, v) == (u != v ? 2 : 0));
            }
    }

    // 森林 / 不連通
    {
        vector<pair<int, int>> edges = {{0, 1}, {1, 2}, {4, 5}};
        LCA s(6, edges, 0);
        assert(s.depth[(size_t)3] == -1);
        assert(s.depth[(size_t)4] == -1);
        assert(s.lca(0, 4) == -1);
        assert(s.dist(0, 4) == -1);
        assert(s.lca(0, 2) == 0);
        assert(s.dist(0, 2) == 2);
        assert(s.kthAncestor(4, 0) == -1);
        assert(!s.isAncestor(0, 4));
    }

    // 自環與越界邊被忽略
    {
        vector<pair<int, int>> edges = {{0, 0}, {0, 1}, {1, 2}, {1, 99}};
        LCA s(3, edges, 0);
        assert((s.depth == vector<int>{0, 1, 2}));
        assert(s.lca(1, 2) == 1);
    }

    // 空樹
    {
        LCA s(0, {}, 0);
        assert(s.lca(0, 0) == -1);
        assert(s.dist(0, 0) == -1);
    }

    // 越界查詢
    {
        vector<pair<int, int>> edges = {{0, 1}, {1, 2}};
        LCA s(3, edges, 0);
        assert(s.lca(0, 7) == -1);
        assert(s.lca(-1, 1) == -1);
        assert(s.dist(0, 9) == -1);
        assert(s.kthAncestor(9, 1) == -1);
        assert(s.kthAncestor(1, -1) == -1);
        assert(!s.isAncestor(0, 9));
    }

    // 隨機對拍：三種路線一致 + k 級祖先 + 祖先判定
    rngEngine.seed(20261005);
    for (int t = 0; t < 120; ++t) {
        int n = rndInt(1, 14);
        vector<pair<int, int>> edges = buildRandomTree(n);
        LCA s(n, edges, 0);
        assert(s.depth[(size_t)0] == 0);
        for (int d : s.depth) assert(d >= 0);                 // 隨機 parent 建出來必連通
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                int w = s.lca(u, v);
                assert(w == s.lcaRmq(u, v) && w == s.lcaNaive(u, v));
                assert(s.dist(u, v) == s.depth[(size_t)u] + s.depth[(size_t)v] - 2 * s.depth[(size_t)w]);
            }
            for (int k = 0; k < n + 2; ++k) assert(s.kthAncestor(u, k) == s.kthAncestorNaive(u, k));
            for (int v = 0; v < n; ++v) {
                bool brute = false;
                int cur = v;
                while (cur >= 0) {
                    if (cur == u) { brute = true; break; }
                    cur = s.up[0][(size_t)cur];
                }
                assert(s.isAncestor(u, v) == brute);
            }
        }
    }

    cout << "all tests passed" << '\n';
}

int main() {
    vector<string> toks;
    string t;
    while (cin >> t) toks.push_back(t);
    if (toks.empty()) runTests();
    else runIo(toks);
    return 0;
}
