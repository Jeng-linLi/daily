// 樹形 DP 與換根 DP（最大獨立集 / 直徑 / 重心 / 各點最遠距離 / 各點距離和）
//
// 題意：
//     給定一棵樹（或森林，無向無環圖，n 個點 m 條邊），求：
//       1. **最大獨立集**：最大的點集，使集合中任意兩點之間沒有邊；
//       2. **樹的直徑**：最長路徑的長度（邊數）與具體路徑；
//       3. **樹的重心**：刪掉該點後，剩下的最大連通塊大小 ≤ n/2 的點（可能 1 或 2 個）；
//       4. **換根 DP（rerooting）之 maxdist**：對每個點 u，求 u 到「所在連通分量內最遠點」的距離；
//       5. **換根 DP 之 sumdist**：對每個點 u，求 u 到「所在連通分量內所有點」的距離之和。
//
// 思路：
//     ### 樹形 DP 的基本套路
//     任取一點為根，把「無根樹」轉成「有根樹」，於是可以自底向上（後序）做 DP：
//     每個點 u 的狀態只依賴它的子節點。`dp[u] = 合併(所有子節點 v 的 dp[v])`。
//     複雜度通常是 O(n)，因為每條邊只被處理常數次。
//
//     ### 最大獨立集
//     每個點有「選 / 不選」兩種狀態：
//         dp1[u] = 1 + Σ dp0[v]                （選 u ⟹ 所有子節點都不能選）
//         dp0[u] = Σ max(dp0[v], dp1[v])       （不選 u ⟹ 子節點隨意）
//     答案 max(dp0[root], dp1[root])。這是「樹上背包 / 狀態機 DP」最簡單的一例。
//     方案重建：父節點被選則子節點強制不選；否則取 dp1 ≥ dp0 者（平手優先選）。
//
//     ### 樹的直徑（兩次 BFS / DFS）
//     從任意點 s 出發走到最遠的點 a，再從 a 出發走到最遠的點 b，
//     則 a–b 必定是一條直徑。對森林則對每個連通分量各做一次並取最大。
//     （這個性質對**樹**成立；一般圖不成立，一般圖要用 Floyd 或 Johnson。）
//
//     ### 樹的重心
//     以任意點為根，記 sz[u] 為 u 的子樹大小。刪掉 u 後會分成若干塊：
//     每個子節點 v 對應一塊大小 sz[v]，父親方向對應一塊大小 compSize - sz[u]。
//     取這些塊的最大值，若 ≤ compSize/2 則 u 是重心。重心必存在且最多 2 個。
//
//     ### 換根 DP（本題重點）
//     上面那些量若「以每個點為根各算一次」是 O(n^2)。換根 DP 把它降到 O(n)：
//     先做一次自底向上（bottom-up）拿到「向下」的資訊，再做一次自頂向下（top-down）
//     把「向上」的資訊補給每個子節點，於是每個點都同時擁有「往下看」和「往上看」的視角。
//
//     - **maxdist**：對每個 u 維護向下最長的兩條路徑 `best1 / best2`（分別記是從哪個子節點來的）。
//       對子節點 v，其「向上的最遠距離」為 `1 + max(up[u], best_excluding_v[u])`，
//       其中 best_excluding_v 就是「不從 v 這個分支走」的最佳值（who1[u] == v 時取 best2，否則取 best1）。
//       答案 `maxdist[u] = max(best1[u], up[u])`。
//     - **sumdist**：記 `sz[u]` 與 `sub[u]`（u 到子樹內所有點的距離和）。
//       經典公式：把根從 u 換到子節點 v 時，v 這一側的 sz[v] 個點距離各減 1，
//       其餘 compSize - sz[v] 個點距離各加 1，所以
//           sumdist[v] = sumdist[u] - sz[v] + (compSize - sz[v])
//       O(1) 轉移，一趟 top-down 全部算完。
//
//     森林的情況：對每個連通分量獨立計算，compSize 用該分量的大小。
//
// 應用場景：
//     社交網路中「影響力最大」的節點（距離和最小 = 最接近所有人的點，即重心類指標）、
//     伺服器 / 倉庫選址（樹上最小化最大延遲 = 直徑中心）、
//     依賴樹上的任務調度（最大獨立集 = 互不衝突的最大任務集）、
//     文件目錄樹與組織架構的聚合統計。
//
// 複雜度：
//     最大獨立集 / 重心 / 換根 DP    O(n) 時間、O(n) 空間
//     直徑（兩次 BFS，對每個分量）   O(n) 時間（森林亦然）、O(n) 空間
//     暴力基準（從每個點 BFS）       O(n · (n + m))
//
// 輸入格式（stdin，全部以空白分隔）：
//     n m
//     m 行：u v        （0-indexed 無向邊；重邊與自環會被忽略；越界邊忽略）
//     n = 0 或輸入不足時視為空圖。遇到非整數 token 視為輸入結束。
// 輸出格式（stdout）：
//     第 1 行：最大獨立集大小
//     第 2 行：獨立集節點（升序，空格分隔）
//     第 3 行：直徑長度（邊數；空圖為 0）
//     第 4 行：直徑路徑（空格分隔；空圖輸出空行）
//     第 5 行：重心個數
//     第 6 行：重心（升序，空格分隔）
//     第 7 行：maxdist[0..n-1]（空格分隔）
//     第 8 行：sumdist[0..n-1]（空格分隔）
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <deque>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// 嚴格整數規則：只接受 [+-]?digits，與 Python 的 re 版完全一致
static bool tryLL(const string &tok, long long &out) {
    if (tok.empty()) return false;
    size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i >= tok.size()) return false;
    for (size_t k = i; k < tok.size(); ++k) {
        if (!isdigit((unsigned char)tok[k])) return false;
    }
    bool neg = (tok[0] == '-');
    string digits = tok.substr(i);
    while (digits.size() > 1 && digits[0] == '0') digits.erase(digits.begin());
    if (digits.size() > 18) digits = digits.substr(digits.size() - 18);   // 保底截斷，避免 UB
    long long v = 0;
    for (char c : digits) v = v * 10 + (c - '0');
    out = neg ? -v : v;
    return true;
}

static vector<vector<int>> build_adj(int n, const vector<pair<int, int>> &edges) {
    set<pair<int, int>> st;
    for (auto e : edges) {
        int u = e.first, v = e.second;
        if (u == v) continue;
        if (u < 0 || u >= n || v < 0 || v >= n) continue;
        if (u > v) swap(u, v);
        st.insert(make_pair(u, v));
    }
    vector<vector<int>> adj(n);
    for (auto pr : st) {
        adj[pr.first].push_back(pr.second);
        adj[pr.second].push_back(pr.first);
    }
    for (int i = 0; i < n; ++i) sort(adj[i].begin(), adj[i].end());
    return adj;
}

struct ForestData {
    vector<int> parent;   // parent[root] = -2
    vector<int> root;
    vector<int> order;    // 先序：父必在子之前
};

static ForestData forest_data(int n, const vector<vector<int>> &adj) {
    ForestData fd;
    fd.parent.assign(n, -1);
    fd.root.resize(n);
    for (int i = 0; i < n; ++i) fd.root[i] = i;
    vector<int> stack;
    for (int s = 0; s < n; ++s) {
        if (fd.parent[s] != -1) continue;
        fd.parent[s] = -2;
        stack.clear();
        stack.push_back(s);
        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            fd.order.push_back(u);
            for (int v : adj[u]) {
                if (fd.parent[v] == -1) {
                    fd.parent[v] = u;
                    fd.root[v] = s;
                    stack.push_back(v);
                }
            }
        }
    }
    return fd;
}

// ---------------------------------------------------------------- 最大獨立集

static pair<int, vector<int>> max_independent_set(int n, const vector<vector<int>> &adj,
                                                  const vector<int> &parent,
                                                  const vector<int> &order) {
    vector<int> dp0(n, 0), dp1(n, 0);
    for (int i = (int)order.size() - 1; i >= 0; --i) {
        int u = order[i];
        dp1[u] = 1;
        for (int v : adj[u]) {
            if (parent[v] == u) {
                dp1[u] += dp0[v];
                dp0[u] += (dp0[v] >= dp1[v] ? dp0[v] : dp1[v]);
            }
        }
    }
    vector<char> chosen(n, 0);
    for (int u : order) {
        int p = parent[u];
        bool blocked = (p >= 0 && chosen[p]);
        chosen[u] = (!blocked && dp1[u] >= dp0[u]) ? 1 : 0;
    }
    vector<int> nodes;
    for (int u = 0; u < n; ++u) if (chosen[u]) nodes.push_back(u);
    return make_pair((int)nodes.size(), nodes);
}

static int mis_bruteforce(int n, const vector<vector<int>> &adj) {
    int best = 0;
    for (int mask = 0; mask < (1 << n); ++mask) {
        bool ok = true;
        for (int u = 0; u < n && ok; ++u) {
            if (!(mask >> u & 1)) continue;
            for (int v : adj[u]) {
                if (mask >> v & 1) { ok = false; break; }
            }
        }
        if (ok) best = max(best, __builtin_popcount((unsigned)mask));
    }
    return best;
}

// ---------------------------------------------------------------- 直徑 / 重心

struct BfsResult {
    int far;
    vector<int> dist;
    vector<int> par;
};

static BfsResult bfs_far(int n, const vector<vector<int>> &adj, int s) {
    vector<int> dist(n, -1), par(n, -1);
    deque<int> q;
    dist[s] = 0;
    q.push_back(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        for (int v : adj[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                par[v] = u;
                q.push_back(v);
            }
        }
    }
    int best = s;
    for (int i = 0; i < n; ++i) if (dist[i] > dist[best]) best = i;
    BfsResult r;
    r.far = best;
    r.dist = dist;
    r.par = par;
    return r;
}

static pair<int, vector<int>> tree_diameter(int n, const vector<vector<int>> &adj,
                                            const vector<int> &root) {
    int best_len = -1;                     // 哨兵：確保第一個分量必定被採用
    vector<int> best_path;
    for (int s = 0; s < n; ++s) {
        if (root[s] != s) continue;
        int a = bfs_far(n, adj, s).far;
        BfsResult br = bfs_far(n, adj, a);
        int b = br.far;
        int length = br.dist[b];
        vector<int> path;
        int cur = b;
        while (cur != -1) {                // 沿 parent 從 b 一路走回 a
            path.push_back(cur);
            if (cur == a) break;
            cur = br.par[cur];
        }
        if (length > best_len) {
            best_len = length;
            best_path = path;
        }
    }
    return make_pair(best_len < 0 ? 0 : best_len, best_path);
}

static vector<int> centroids(int n, const vector<vector<int>> &adj, const vector<int> &parent,
                             const vector<int> &root, const vector<int> &sz,
                             const vector<int> &comp) {
    vector<int> res;
    for (int u = 0; u < n; ++u) {
        int mx = comp[root[u]] - sz[u];
        for (int v : adj[u]) {
            if (parent[v] == u && sz[v] > mx) mx = sz[v];
        }
        if (mx * 2 <= comp[root[u]]) res.push_back(u);
    }
    return res;
}

// 重心暴力版：真的把 u 刪掉，再數 u 所在分量內各連通塊的大小。
// 注意森林要**逐分量**判定：只統計與 u 同分量的塊，且與 comp[root[u]] 比較。
// 若把其它分量也算進來，塊大小會被無關的分量撐大，定義就錯了。
static vector<int> centroids_bruteforce(int n, const vector<vector<int>> &adj,
                                        const vector<int> &comp, const vector<int> &root) {
    vector<int> res;
    for (int u = 0; u < n; ++u) {
        int r = root[u];
        vector<char> seen(n, 0);
        seen[u] = 1;
        int mx = 0;
        for (int s = 0; s < n; ++s) {
            if (seen[s] || root[s] != r || s == u) continue;
            int cnt = 0;
            deque<int> q;
            seen[s] = 1;
            q.push_back(s);
            while (!q.empty()) {
                int x = q.front();
                q.pop_front();
                ++cnt;
                for (int y : adj[x]) {
                    if (!seen[y] && y != u && root[y] == r) {
                        seen[y] = 1;
                        q.push_back(y);
                    }
                }
            }
            if (cnt > mx) mx = cnt;
        }
        if (mx * 2 <= comp[r]) res.push_back(u);
    }
    return res;
}

// ---------------------------------------------------------------- 換根 DP

struct RerootResult {
    vector<int> maxdist;
    vector<long long> sumdist;
    vector<int> sz;
    vector<int> comp;
};

static RerootResult reroot(int n, const vector<vector<int>> &adj, const vector<int> &parent,
                           const vector<int> &root, const vector<int> &order) {
    RerootResult rr;
    rr.comp.assign(n, 0);
    for (int u = 0; u < n; ++u) rr.comp[root[u]] += 1;

    vector<int> best1(n, 0), best2(n, 0), who1(n, -1);
    rr.sz.assign(n, 1);
    vector<long long> sub(n, 0);

    for (int i = (int)order.size() - 1; i >= 0; --i) {
        int u = order[i];
        for (int v : adj[u]) {
            if (parent[v] == u) {
                rr.sz[u] += rr.sz[v];
                sub[u] += sub[v] + rr.sz[v];
                int val = best1[v] + 1;
                if (val > best1[u]) {
                    best2[u] = best1[u];
                    best1[u] = val;
                    who1[u] = v;
                } else if (val > best2[u]) {
                    best2[u] = val;
                }
            }
        }
    }

    vector<int> up(n, 0);
    rr.maxdist.assign(n, 0);
    rr.sumdist.assign(n, 0);
    for (int u : order) {
        rr.maxdist[u] = max(best1[u], up[u]);
        if (parent[u] == -2) rr.sumdist[u] = sub[u];
        for (int v : adj[u]) {
            if (parent[v] == u) {
                int excl = (who1[u] == v ? best2[u] : best1[u]);
                int via_up = max(up[u], excl);
                up[v] = 1 + via_up;
                rr.sumdist[v] = rr.sumdist[u] - rr.sz[v] + (long long)(rr.comp[root[u]] - rr.sz[v]);
            }
        }
    }
    return rr;
}

static pair<vector<int>, vector<long long>> brute_distances(int n, const vector<vector<int>> &adj) {
    vector<int> md(n, 0);
    vector<long long> sd(n, 0);
    for (int s = 0; s < n; ++s) {
        BfsResult br = bfs_far(n, adj, s);
        int mx = 0;
        long long total = 0;
        for (int i = 0; i < n; ++i) {
            if (br.dist[i] > 0) {
                if (br.dist[i] > mx) mx = br.dist[i];
                total += br.dist[i];
            }
        }
        md[s] = mx;
        sd[s] = total;
    }
    return make_pair(md, sd);
}

// ---------------------------------------------------------------- IO 與測試

static string join_ints(const vector<int> &v) {
    string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += " ";
        out += to_string(v[i]);
    }
    return out;
}

static string join_ll(const vector<long long> &v) {
    string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += " ";
        out += to_string(v[i]);
    }
    return out;
}

static void run_io(const string &raw) {
    vector<string> toks;
    {
        istringstream iss(raw);
        string t;
        while (iss >> t) toks.push_back(t);
    }
    size_t pos = 0;
    auto nxt = [&]() -> long long {
        long long v = 0;
        if (pos < toks.size()) {
            long long parsed = 0;
            if (tryLL(toks[pos], parsed)) v = parsed;
        }
        ++pos;
        return v;
    };

    long long nll = nxt(), mll = nxt();
    int n = (int)nll, m = (int)mll;
    if (n < 0) n = 0;
    if (m < 0) m = 0;
    vector<pair<int, int>> edges;
    for (int i = 0; i < m; ++i) {
        long long u = nxt(), v = nxt();
        edges.push_back(make_pair((int)u, (int)v));
    }

    vector<vector<int>> adj = build_adj(n, edges);
    ForestData fd = forest_data(n, adj);
    pair<int, vector<int>> mis = max_independent_set(n, adj, fd.parent, fd.order);
    pair<int, vector<int>> dia = tree_diameter(n, adj, fd.root);
    RerootResult rr = reroot(n, adj, fd.parent, fd.root, fd.order);
    vector<int> cens = centroids(n, adj, fd.parent, fd.root, rr.sz, rr.comp);

    vector<string> out;
    out.push_back(to_string(mis.first));
    out.push_back(join_ints(mis.second));
    out.push_back(to_string(dia.first));
    out.push_back(join_ints(dia.second));
    out.push_back(to_string((int)cens.size()));
    out.push_back(join_ints(cens));
    out.push_back(join_ints(rr.maxdist));
    out.push_back(join_ll(rr.sumdist));
    for (size_t i = 0; i < out.size(); ++i) cout << out[i] << "\n";
}

static void run_tests() {
    // ---- 固定用例：鏈 0-1-2-3 ----
    {
        vector<pair<int, int>> chain{{0, 1}, {1, 2}, {2, 3}};
        vector<vector<int>> adj = build_adj(4, chain);
        ForestData fd = forest_data(4, adj);
        pair<int, vector<int>> mis = max_independent_set(4, adj, fd.parent, fd.order);
        assert(mis.first == 2);
        assert(mis.second == vector<int>({0, 2}) || mis.second == vector<int>({1, 3}));
        assert(mis_bruteforce(4, adj) == 2);
        assert(tree_diameter(4, adj, fd.root) == make_pair(3, vector<int>({0, 1, 2, 3})));
        RerootResult rr = reroot(4, adj, fd.parent, fd.root, fd.order);
        assert(rr.maxdist == vector<int>({3, 2, 2, 3}));
        assert(rr.sumdist == vector<long long>({6, 4, 4, 6}));
        assert(brute_distances(4, adj) == make_pair(rr.maxdist, rr.sumdist));
        assert(centroids(4, adj, fd.parent, fd.root, rr.sz, rr.comp) == vector<int>({1, 2}));
        assert(centroids_bruteforce(4, adj, rr.comp, fd.root) == vector<int>({1, 2}));
    }

    // ---- 固定用例：星形（中心 0，葉 1..4）----
    {
        vector<pair<int, int>> star;
        for (int i = 1; i <= 4; ++i) star.push_back(make_pair(0, i));
        vector<vector<int>> adj = build_adj(5, star);
        ForestData fd = forest_data(5, adj);
        pair<int, vector<int>> mis = max_independent_set(5, adj, fd.parent, fd.order);
        assert(mis.first == 4 && mis.second == vector<int>({1, 2, 3, 4}));
        assert(tree_diameter(5, adj, fd.root) == make_pair(2, vector<int>({2, 0, 1})));
        RerootResult rr = reroot(5, adj, fd.parent, fd.root, fd.order);
        assert(rr.maxdist == vector<int>({1, 2, 2, 2, 2}));
        assert(rr.sumdist == vector<long long>({4, 7, 7, 7, 7}));
        assert(brute_distances(5, adj) == make_pair(rr.maxdist, rr.sumdist));
        assert(centroids(5, adj, fd.parent, fd.root, rr.sz, rr.comp) == vector<int>({0}));
    }

    // ---- 邊界：空圖 / 單點 / 森林 ----
    {
        vector<vector<int>> adj = build_adj(0, vector<pair<int, int>>());
        ForestData fd = forest_data(0, adj);
        assert(max_independent_set(0, adj, fd.parent, fd.order).first == 0);
        assert(tree_diameter(0, adj, fd.root) == make_pair(0, vector<int>()));
        assert(reroot(0, adj, fd.parent, fd.root, fd.order).maxdist.empty());
    }
    {
        vector<vector<int>> adj = build_adj(1, vector<pair<int, int>>());
        ForestData fd = forest_data(1, adj);
        assert(max_independent_set(1, adj, fd.parent, fd.order) == make_pair(1, vector<int>({0})));
        assert(tree_diameter(1, adj, fd.root) == make_pair(0, vector<int>({0})));
        RerootResult rr = reroot(1, adj, fd.parent, fd.root, fd.order);
        assert(rr.maxdist == vector<int>({0}) && rr.sumdist == vector<long long>({0}));
        assert(centroids(1, adj, fd.parent, fd.root, rr.sz, rr.comp) == vector<int>({0}));
    }
    {
        vector<pair<int, int>> forest{{0, 1}, {2, 3}};
        vector<vector<int>> adj = build_adj(4, forest);
        ForestData fd = forest_data(4, adj);
        assert(max_independent_set(4, adj, fd.parent, fd.order).first == 2);  // 兩條邊各取一端
        RerootResult rr = reroot(4, adj, fd.parent, fd.root, fd.order);
        assert(rr.maxdist == vector<int>({1, 1, 1, 1}));
        assert(rr.sumdist == vector<long long>({1, 1, 1, 1}));
        assert(brute_distances(4, adj) == make_pair(rr.maxdist, rr.sumdist));
    }

    mt19937 rng(20261007);

    // ---- 隨機對拍：隨機樹（n ≤ 12）----
    for (int t = 0; t < 400; ++t) {
        int n = 1 + (int)(rng() % 12);
        vector<int> perm(n);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), rng);
        vector<pair<int, int>> edges;
        for (int i = 1; i < n; ++i) {          // 隨機父節點 → 保證是一棵樹
            int j = (int)(rng() % i);
            edges.push_back(make_pair(perm[i], perm[j]));
        }
        vector<vector<int>> adj = build_adj(n, edges);
        ForestData fd = forest_data(n, adj);

        pair<int, vector<int>> mis = max_independent_set(n, adj, fd.parent, fd.order);
        assert(mis.first == mis_bruteforce(n, adj));
        {
            vector<char> in(n, 0);
            for (int x : mis.second) in[x] = 1;
            assert((int)mis.second.size() == mis.first);
            for (int u : mis.second)
                for (int v : adj[u]) assert(!in[v]);   // 驗證確實是獨立集
        }

        RerootResult rr = reroot(n, adj, fd.parent, fd.root, fd.order);
        pair<vector<int>, vector<long long>> bd = brute_distances(n, adj);
        assert(rr.maxdist == bd.first);
        assert(rr.sumdist == bd.second);

        assert(centroids(n, adj, fd.parent, fd.root, rr.sz, rr.comp) ==
               centroids_bruteforce(n, adj, rr.comp, fd.root));

        pair<int, vector<int>> dia = tree_diameter(n, adj, fd.root);
        assert(dia.first == (int)dia.second.size() - 1);
        {
            vector<int> uniq = dia.second;
            sort(uniq.begin(), uniq.end());
            uniq.erase(unique(uniq.begin(), uniq.end()), uniq.end());
            assert((int)uniq.size() == (int)dia.second.size());
            for (int i = 0; i + 1 < (int)dia.second.size(); ++i) {
                int a = dia.second[i], b = dia.second[i + 1];
                assert(find(adj[a].begin(), adj[a].end(), b) != adj[a].end());
            }
        }
        assert(dia.first == *max_element(rr.maxdist.begin(), rr.maxdist.end()));
    }

    // ---- 隨機對拍：隨機森林（含孤立點）----
    for (int t = 0; t < 200; ++t) {
        int n = 1 + (int)(rng() % 12);
        int m = (int)(rng() % (n + 1));
        vector<pair<int, int>> raw;
        for (int i = 0; i < m; ++i) {
            int u = (int)(rng() % n), v = (int)(rng() % n);
            if (u != v) raw.push_back(make_pair(u, v));
        }
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        function<int(int)> uf_find = [&](int x) { return p[x] == x ? x : p[x] = uf_find(p[x]); };
        vector<pair<int, int>> keep;
        for (auto e : raw) {
            int a = uf_find(e.first), b = uf_find(e.second);
            if (a != b) { p[a] = b; keep.push_back(e); }
        }
        vector<vector<int>> adj = build_adj(n, keep);
        ForestData fd = forest_data(n, adj);
        RerootResult rr = reroot(n, adj, fd.parent, fd.root, fd.order);
        pair<vector<int>, vector<long long>> bd = brute_distances(n, adj);
        assert(rr.maxdist == bd.first);
        assert(rr.sumdist == bd.second);
        assert(max_independent_set(n, adj, fd.parent, fd.order).first == mis_bruteforce(n, adj));
        assert(centroids(n, adj, fd.parent, fd.root, rr.sz, rr.comp) ==
               centroids_bruteforce(n, adj, rr.comp, fd.root));
    }

    cout << "all tests passed" << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ostringstream oss;
    oss << cin.rdbuf();
    string raw = oss.str();

    bool has_input = raw.find_first_not_of(" \t\r\n") != string::npos;
    if (has_input) run_io(raw);
    else run_tests();
    return 0;
}
