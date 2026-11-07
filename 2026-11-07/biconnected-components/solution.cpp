// 雙連通分量：邊雙連通分量（e-DCC）、點雙連通分量（v-DCC）、橋、割點
//
// 題意：給定一張簡單無向圖（n 點 m 邊，無自環、無重邊），求：
//         1. 橋（刪掉後圖不連通的邊）；
//         2. 割點（刪掉後圖不連通的點）；
//         3. 邊雙連通分量（= 刪去所有橋後的連通分量，孤立點各自成一個分量）；
//         4. 點雙連通分量（割點可屬於多個 v-DCC，橋兩端構成大小 2 的
//            v-DCC，孤立點自成一個 v-DCC）。
//
// 思路（Tarjan，迭代版避免深遞迴棧溢出）：
//   橋：一次 DFS 求 low/disc，樹邊 (p,u) 是橋 iff low[u] > disc[p]。
//   邊雙：刪橋後 DFS 劃分分量（含孤立點），分量按「最小點編號」重新編號
//         ——與 DFS 順序無關，保證跨語言輸出一致。
//   點雙：維護點棧，child c 完成時若 low[c] >= disc[u]，把棧中 c 以上的
//         點連同 u 彈出成一個分量；根若無子樹分量則自成單點分量。
//   割點：u 是割點 iff（u 是 DFS 根且子樹數 >= 2）或存在子樹 c 使
//         low[c] >= disc[u]（根單獨判斷）。
//
// IO 模式（stdin 有內容時）：
//   n m
//   m 行：u v
//   輸出（依序）：
//     行 1：邊雙分量數 e
//     行 2：n 個整數，每點的邊雙編號（按分量最小點升序編號）
//     行 3：點雙分量數 v
//     行 4..：v 行，每行「大小 + 升序點列」（分量按字典序升序）
//     行：橋數 b；接著 b 行「u v」（端點升序、橋之間字典序升序）
//     行：割點數 a；接著 1 行割點升序點列（a=0 時為空行）
//   無 stdin 時跑內建斷言並輸出 all tests passed。

// 註：不用 bits/stdc++.h——本機 MinGW 8.1 的 <filesystem> 頭有已知 bug，
// 顯式列出所需頭文件更穩。
#include <cstdio>
#include <cstring>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <utility>
using namespace std;

struct Tarjan {
    int n;
    vector<vector<pair<int, int>>> adj;  // (鄰點, 邊編號)
    vector<int> disc, low, ptr;
    int timer;

    explicit Tarjan(int n_) : n(n_), adj(n_), disc(n_, -1), low(n_, 0),
                              ptr(n_, 0), timer(0) {}

    void addEdge(int u, int v, int id) {
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
    }

    // 一次迭代 DFS：同時求橋（low[u] > disc[p]）與割點標記。
    // isBridgeEdge / isCut 由呼叫端按邊數 / 點數預先 resize。
    void bridgesAndCuts(vector<int>& isBridgeEdge, vector<char>& isCut) {
        for (int s = 0; s < n; ++s) {
            if (disc[s] != -1) continue;
            disc[s] = low[s] = timer++;
            int rootChildren = 0;
            vector<pair<int, int>> stk;  // (點, 進入邊編號)
            stk.push_back({s, -1});
            while (!stk.empty()) {
                auto [u, pe] = stk.back();
                bool advanced = false;
                while (ptr[u] < (int)adj[u].size()) {
                    auto [v, eid] = adj[u][ptr[u]];
                    ++ptr[u];
                    if (eid == pe) continue;
                    if (disc[v] == -1) {
                        disc[v] = low[v] = timer++;
                        stk.push_back({v, eid});
                        advanced = true;
                        if (u == s) ++rootChildren;
                        break;
                    }
                    if (disc[v] < low[u]) low[u] = disc[v];
                }
                if (!advanced) {
                    stk.pop_back();
                    if (!stk.empty()) {
                        int p = stk.back().first;
                        if (low[u] < low[p]) low[p] = low[u];
                        if (low[u] > disc[p]) isBridgeEdge[pe] = 1;
                        if (low[u] >= disc[p] && p != s) isCut[p] = 1;
                    }
                }
            }
            if (rootChildren >= 2) isCut[s] = 1;
        }
    }
};

// 邊雙：刪橋後 DFS 劃分分量 + 按最小點編號重編號
vector<int> edgeDcc(int n, const vector<pair<int, int>>& edges,
                    const vector<int>& isBridgeEdge) {
    vector<vector<pair<int, int>>> adj(n);
    for (size_t i = 0; i < edges.size(); ++i) {
        if (isBridgeEdge[i]) continue;
        adj[edges[i].first].push_back({edges[i].second, (int)i});
        adj[edges[i].second].push_back({edges[i].first, (int)i});
    }
    vector<int> comp(n, -1);
    int nc = 0;
    for (int s = 0; s < n; ++s) {
        if (comp[s] != -1) continue;
        comp[s] = nc;
        vector<int> st{s};
        while (!st.empty()) {
            int u = st.back();
            st.pop_back();
            for (auto [v, eid] : adj[u]) {
                if (comp[v] == -1) {
                    comp[v] = nc;
                    st.push_back(v);
                }
            }
        }
        ++nc;
    }
    vector<int> minv(nc, INT_MAX);
    for (int v = 0; v < n; ++v) minv[comp[v]] = min(minv[comp[v]], v);
    vector<int> order(nc);
    for (int c = 0; c < nc; ++c) order[c] = c;
    sort(order.begin(), order.end(),
         [&](int a, int b) { return minv[a] < minv[b]; });
    vector<int> remap(nc);
    for (int i = 0; i < nc; ++i) remap[order[i]] = i;
    vector<int> res(n);
    for (int v = 0; v < n; ++v) res[v] = remap[comp[v]];
    return res;
}

// 點雙（點棧 Tarjan）：回傳分量列表（分量內升序、分量間字典序升序）
vector<vector<int>> vertexDcc(int n, const vector<pair<int, int>>& edges) {
    vector<vector<pair<int, int>>> adj(n);
    for (size_t i = 0; i < edges.size(); ++i) {
        adj[edges[i].first].push_back({edges[i].second, (int)i});
        adj[edges[i].second].push_back({edges[i].first, (int)i});
    }
    vector<int> disc(n, -1), low(n, 0), ptr(n, 0);
    int timer = 0;
    vector<vector<int>> comps;
    vector<int> vstk;
    for (int s = 0; s < n; ++s) {
        if (disc[s] != -1) continue;
        size_t created = comps.size();
        disc[s] = low[s] = timer++;
        vstk.push_back(s);
        vector<pair<int, int>> stk{{s, -1}};
        while (!stk.empty()) {
            auto [u, pe] = stk.back();
            bool advanced = false;
            while (ptr[u] < (int)adj[u].size()) {
                auto [v, eid] = adj[u][ptr[u]];
                ++ptr[u];
                if (eid == pe) continue;
                if (disc[v] == -1) {
                    disc[v] = low[v] = timer++;
                    vstk.push_back(v);
                    stk.push_back({v, eid});
                    advanced = true;
                    break;
                }
                if (disc[v] < low[u]) low[u] = disc[v];
            }
            if (!advanced) {
                stk.pop_back();
                if (!stk.empty()) {
                    int p = stk.back().first;
                    if (low[u] < low[p]) low[p] = low[u];
                    if (low[u] >= disc[p]) {
                        // p 是割點（或根）：彈出 u 及其以上 + p 成一個分量
                        vector<int> comp;
                        while (true) {
                            int w = vstk.back();
                            vstk.pop_back();
                            comp.push_back(w);
                            if (w == u) break;
                        }
                        comp.push_back(p);
                        sort(comp.begin(), comp.end());
                        comps.push_back(comp);
                    }
                }
            }
        }
        if (comps.size() == created) {
            // 根孤立（整棵 DFS 沒產生任何分量）：自成單點分量
            comps.push_back({vstk.back()});
            vstk.pop_back();
        } else {
            // 根已含在子分量中，棧底剩餘的根是重複，丟棄
            vstk.pop_back();
        }
    }
    sort(comps.begin(), comps.end());
    return comps;
}

// ---------- 內建測試（與指數暴力對拍） ----------

struct Brute {
    int n;
    vector<pair<int, int>> edges;
    vector<vector<int>> adj;

    Brute(int n_, const vector<pair<int, int>>& es) : n(n_), edges(es), adj(n_) {
        for (auto& e : es) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
    }

    // vs 內部連通分量數（只走 vs 內部邊，不經過 vs 以外的點）
    int compCount(const vector<int>& vs) const {
        vector<char> in(n, 0);
        for (int v : vs) in[v] = 1;
        vector<char> seen(n, 0);
        int cnt = 0;
        for (int s : vs) {
            if (seen[s]) continue;
            ++cnt;
            vector<int> st{s};
            seen[s] = 1;
            while (!st.empty()) {
                int u = st.back();
                st.pop_back();
                for (int v : adj[u]) {
                    if (in[v] && !seen[v]) {
                        seen[v] = 1;
                        st.push_back(v);
                    }
                }
            }
        }
        return cnt;
    }

    vector<pair<int, int>> bridges() const {
        vector<pair<int, int>> res;
        for (size_t i = 0; i < edges.size(); ++i) {
            vector<vector<int>> ra(n);
            for (size_t j = 0; j < edges.size(); ++j) {
                if (j == i) continue;
                ra[edges[j].first].push_back(edges[j].second);
                ra[edges[j].second].push_back(edges[j].first);
            }
            int u = edges[i].first, v = edges[i].second;
            vector<char> seen(n, 0);
            vector<int> st{u};
            seen[u] = 1;
            while (!st.empty()) {
                int x = st.back();
                st.pop_back();
                for (int y : ra[x])
                    if (!seen[y]) { seen[y] = 1; st.push_back(y); }
            }
            if (!seen[v]) res.push_back({min(u, v), max(u, v)});
        }
        sort(res.begin(), res.end());
        return res;
    }

    vector<int> articulation() const {
        vector<int> all(n);
        for (int i = 0; i < n; ++i) all[i] = i;
        int full = compCount(all);
        vector<int> res;
        for (int x = 0; x < n; ++x) {
            vector<int> rest;
            for (int v = 0; v < n; ++v)
                if (v != x) rest.push_back(v);
            // x 孤立時刪去後分量數減一，閾值也要減一
            if (compCount(rest) > full - (adj[x].empty() ? 1 : 0))
                res.push_back(x);
        }
        return res;
    }

    // 邊雙：刪所有橋後的連通分量 + 孤立點，按最小點編號重編號
    vector<int> edgeDccBrute() const {
        auto br = bridges();
        set<pair<int, int>> bset(br.begin(), br.end());
        vector<int> comp(n, -1);
        int nc = 0;
        for (int s = 0; s < n; ++s) {
            if (comp[s] != -1) continue;
            comp[s] = nc;
            vector<int> st{s};
            while (!st.empty()) {
                int u = st.back();
                st.pop_back();
                for (int v : adj[u]) {
                    pair<int, int> key{min(u, v), max(u, v)};
                    if (!bset.count(key) && comp[v] == -1) {
                        comp[v] = nc;
                        st.push_back(v);
                    }
                }
            }
            ++nc;
        }
        vector<int> minv(nc, INT_MAX);
        for (int v = 0; v < n; ++v) minv[comp[v]] = min(minv[comp[v]], v);
        vector<int> order(nc);
        for (int c = 0; c < nc; ++c) order[c] = c;
        sort(order.begin(), order.end(),
             [&](int a, int b) { return minv[a] < minv[b]; });
        vector<int> remap(nc);
        for (int i = 0; i < nc; ++i) remap[order[i]] = i;
        vector<int> res(n);
        for (int v = 0; v < n; ++v) res[v] = remap[comp[v]];
        return res;
    }

    // 點雙定義暴力：極大的「2-點連通 / 橋對 / 孤立點」子集（n <= 8 枚舉子集）
    vector<vector<int>> vertexDccBrute() const {
        set<pair<int, int>> eset;
        vector<int> deg(n, 0);
        for (auto& e : edges) {
            eset.insert({min(e.first, e.second), max(e.first, e.second)});
            ++deg[e.first];
            ++deg[e.second];
        }
        auto valid = [&](const vector<int>& sub) {
            int k = (int)sub.size();
            if (k == 1) return deg[sub[0]] == 0;
            if (k == 2)
                return eset.count({min(sub[0], sub[1]), max(sub[0], sub[1])}) > 0;
            for (int skip : sub) {
                // 檢查 sub 去掉 skip 後內部連通（只用 sub 內部邊）
                vector<char> in(n, 0);
                for (int v : sub) in[v] = 1;
                vector<vector<int>> sa(n);
                for (auto& e : edges) {
                    if (in[e.first] && in[e.second] && e.first != skip &&
                        e.second != skip) {
                        sa[e.first].push_back(e.second);
                        sa[e.second].push_back(e.first);
                    }
                }
                int start = -1;
                for (int v : sub)
                    if (v != skip) { start = v; break; }
                vector<char> seen(n, 0);
                vector<int> st{start};
                seen[start] = 1;
                int cnt = 1;
                while (!st.empty()) {
                    int u = st.back();
                    st.pop_back();
                    for (int v : sa[u])
                        if (!seen[v]) { seen[v] = 1; ++cnt; st.push_back(v); }
                }
                if (cnt != k - 1) return false;
            }
            return true;
        };
        vector<vector<int>> vsets;
        for (int mask = 1; mask < (1 << n); ++mask) {
            vector<int> sub;
            for (int i = 0; i < n; ++i)
                if (mask >> i & 1) sub.push_back(i);
            if (valid(sub)) vsets.push_back(sub);
        }
        vector<vector<int>> maximal;
        for (auto& a : vsets) {
            set<int> sa(a.begin(), a.end());
            bool dominated = false;
            for (auto& b : vsets) {
                set<int> sb(b.begin(), b.end());
                if (sb.size() > sa.size() &&
                    includes(sb.begin(), sb.end(), sa.begin(), sa.end())) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) maximal.push_back(a);
        }
        sort(maximal.begin(), maximal.end());
        return maximal;
    }
};

static unsigned long long g_seed = 20261107ULL;
static int rnd(int mv) {
    g_seed = g_seed * 1103515245ULL + 12345ULL;
    return (int)((g_seed >> 1) % (unsigned long long)mv);
}

void runAllTests() {
    struct Case { int n; vector<pair<int, int>> es; };
    vector<Case> fixed = {
        {6, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {2, 3}}},  // 兩環一橋
        {4, {{0, 1}, {1, 2}, {2, 3}}},            // 鏈：全橋
        {4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}},    // 四環
        {5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}}},    // 星：中心割點
        {5, {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 2}}},  // 8 字
        {3, {}},                                   // 全孤立
        {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 2}, {0, 3}}},
    };

    auto check = [&](int n, const vector<pair<int, int>>& es) {
        // 橋
        Tarjan tj(n);
        for (size_t i = 0; i < es.size(); ++i) tj.addEdge(es[i].first, es[i].second, (int)i);
        vector<int> isBridge(es.size(), 0);
        vector<char> isCut(n, 0);
        tj.bridgesAndCuts(isBridge, isCut);
        Brute br(n, es);
        auto wantB = br.bridges();
        vector<pair<int, int>> gotB;
        for (size_t i = 0; i < es.size(); ++i)
            if (isBridge[i])
                gotB.push_back({min(es[i].first, es[i].second), max(es[i].first, es[i].second)});
        sort(gotB.begin(), gotB.end());
        assert(gotB == wantB);
        // 割點
        vector<int> cuts;
        for (int v = 0; v < n; ++v)
            if (isCut[v]) cuts.push_back(v);
        assert(cuts == br.articulation());
        // 邊雙
        assert(edgeDcc(n, es, isBridge) == br.edgeDccBrute());
        // 點雙
        assert(vertexDcc(n, es) == br.vertexDccBrute());
    };

    for (auto& c : fixed) check(c.n, c.es);

    // 隨機小圖（n <= 7）與全暴力對拍
    for (int t = 0; t < 250; ++t) {
        int n = rnd(7) + 1;
        vector<pair<int, int>> pairs;
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b) pairs.push_back({a, b});
        int m = rnd((int)pairs.size() + 1);
        vector<int> idx(pairs.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = (int)i;
        for (int i = (int)idx.size() - 1; i > 0; --i) {
            int j = rnd(i + 1);
            swap(idx[i], idx[j]);
        }
        vector<pair<int, int>> es;
        for (int i = 0; i < m; ++i) es.push_back(pairs[idx[i]]);
        check(n, es);
    }

    // 中等隨機圖冒煙：邊雙分量內部連通、點雙覆蓋所有點
    for (int t = 0; t < 30; ++t) {
        int n = rnd(60) + 10;
        int m = rnd(2 * n) + n / 2;
        vector<pair<int, int>> es;
        set<pair<int, int>> seen;
        int guard = 0;
        while ((int)es.size() < m && guard < m * 20) {
            ++guard;
            int a = rnd(n), b = rnd(n);
            if (a == b) continue;
            pair<int, int> key{min(a, b), max(a, b)};
            if (seen.count(key)) continue;
            seen.insert(key);
            es.push_back(key);
        }
        Tarjan tj(n);
        for (size_t i = 0; i < es.size(); ++i) tj.addEdge(es[i].first, es[i].second, (int)i);
        vector<int> isBridge(es.size(), 0);
        vector<char> isCut(n, 0);
        tj.bridgesAndCuts(isBridge, isCut);
        auto ecomp = edgeDcc(n, es, isBridge);
        map<int, vector<int>> groups;
        for (int v = 0; v < n; ++v) groups[ecomp[v]].push_back(v);
        for (auto& kv : groups) {
            if (kv.second.size() < 2) continue;
            // 檢查分量內部連通（只用非橋邊）
            vector<char> in(n, 0);
            for (int v : kv.second) in[v] = 1;
            vector<vector<int>> sa(n);
            for (size_t i = 0; i < es.size(); ++i)
                if (!isBridge[i] && in[es[i].first] && in[es[i].second]) {
                    sa[es[i].first].push_back(es[i].second);
                    sa[es[i].second].push_back(es[i].first);
                }
            vector<char> seen(n, 0);
            vector<int> st{kv.second[0]};
            seen[kv.second[0]] = 1;
            int cnt = 1;
            while (!st.empty()) {
                int u = st.back();
                st.pop_back();
                for (int v : sa[u])
                    if (!seen[v]) { seen[v] = 1; ++cnt; st.push_back(v); }
            }
            assert(cnt == (int)kv.second.size());
        }
        auto vcomps = vertexDcc(n, es);
        vector<char> cov(n, 0);
        for (auto& c : vcomps)
            for (int v : c) cov[v] = 1;
        for (int v = 0; v < n; ++v) assert(cov[v]);
    }

    printf("all tests passed\n");
}

int main() {
    // 一次性讀入全部 stdin；無內容時跑內建測試
    std::stringstream buf;
    buf << cin.rdbuf();
    string data = buf.str();
    if (data.find_first_not_of(" \t\r\n") == string::npos) {
        runAllTests();
        return 0;
    }
    std::istringstream in(data);
    int n, m;
    in >> n >> m;
    vector<pair<int, int>> edges(m);
    for (int i = 0; i < m; ++i) in >> edges[i].first >> edges[i].second;

    Tarjan tj(n);
    for (int i = 0; i < m; ++i) tj.addEdge(edges[i].first, edges[i].second, i);
    vector<int> isBridge(m, 0);
    vector<char> isCut(n, 0);
    tj.bridgesAndCuts(isBridge, isCut);

    auto ecomp = edgeDcc(n, edges, isBridge);
    auto vcomps = vertexDcc(n, edges);
    vector<pair<int, int>> bl;
    for (int i = 0; i < m; ++i)
        if (isBridge[i])
            bl.push_back({min(edges[i].first, edges[i].second), max(edges[i].first, edges[i].second)});
    sort(bl.begin(), bl.end());
    vector<int> cuts;
    for (int v = 0; v < n; ++v)
        if (isCut[v]) cuts.push_back(v);

    ostringstream out;
    set<int> distinct(ecomp.begin(), ecomp.end());
    out << distinct.size() << "\n";
    for (int v = 0; v < n; ++v) out << ecomp[v] << " \n"[v == n - 1];
    out << vcomps.size() << "\n";
    for (auto& c : vcomps) {
        out << c.size();
        for (int v : c) out << " " << v;
        out << "\n";
    }
    out << bl.size() << "\n";
    for (auto& p : bl) out << p.first << " " << p.second << "\n";
    out << cuts.size() << "\n";
    for (size_t i = 0; i < cuts.size(); ++i)
        out << cuts[i] << (i + 1 == cuts.size() ? "\n" : " ");
    if (cuts.empty()) out << "\n";
    fputs(out.str().c_str(), stdout);
    return 0;
}
