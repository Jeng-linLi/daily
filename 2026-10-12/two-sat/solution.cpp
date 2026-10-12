// 2-SAT（2-可滿足性）：蘊含圖 + 強連通分量 + 構造解
//
// 與 solution.py 完全同構：同樣的演算法、同樣的確定性規則、同樣的輸入輸出格式。
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef long long ll;
typedef pair<int, int> Clause;

// ---------------------------------------------------------------- SCC 兩種寫法

// Kosaraju 迭代版：兩次 DFS，第二次在反圖上按完成時間由晚到早走。
static vector<int> kosarajuSCC(int n, const vector<vector<int>> &adj) {
    vector<vector<int>> radj(n);
    for (int u = 0; u < n; u++)
        for (int v : adj[u]) radj[v].push_back(u);

    vector<char> visited(n, 0);
    vector<int> order;
    for (int s = 0; s < n; s++) {
        if (visited[s]) continue;
        visited[s] = 1;
        vector<pair<int, int>> stk;
        stk.push_back({s, 0});
        while (!stk.empty()) {
            int u = stk.back().first;
            int i = stk.back().second;
            if (i < (int)adj[u].size()) {
                stk.back().second = i + 1;   // 顯式棧不持有引用，避免 vector 重配置失效
                int v = adj[u][i];
                if (!visited[v]) {
                    visited[v] = 1;
                    stk.push_back({v, 0});
                }
            } else {
                order.push_back(u);
                stk.pop_back();
            }
        }
    }

    vector<int> comp(n, -1);
    int c = 0;
    for (int idx = n - 1; idx >= 0; idx--) {   // order.size() == n
        int s = order[idx];
        if (comp[s] != -1) continue;
        comp[s] = c;
        vector<int> stk;
        stk.push_back(s);
        while (!stk.empty()) {
            int u = stk.back();
            stk.pop_back();
            for (int v : radj[u])
                if (comp[v] == -1) {
                    comp[v] = c;
                    stk.push_back(v);
                }
        }
        c++;
    }
    return comp;
}

// Tarjan 迭代版：dfn / low + 顯式棧，返回 comp[v]。
static vector<int> tarjanSCC(int n, const vector<vector<int>> &adj) {
    vector<int> dfn(n, 0), low(n, 0), comp(n, -1);
    vector<char> on(n, 0);
    vector<int> stk;
    int timer = 1, nc = 0;

    for (int s = 0; s < n; s++) {
        if (dfn[s]) continue;
        dfn[s] = low[s] = timer++;
        stk.push_back(s);
        on[s] = 1;
        vector<pair<int, int>> work;
        work.push_back({s, 0});
        while (!work.empty()) {
            int u = work.back().first;
            int i = work.back().second;
            if (i < (int)adj[u].size()) {
                work.back().second = i + 1;
                int v = adj[u][i];
                if (!dfn[v]) {
                    dfn[v] = low[v] = timer++;
                    stk.push_back(v);
                    on[v] = 1;
                    work.push_back({v, 0});
                } else if (on[v]) {
                    low[u] = min(low[u], dfn[v]);
                }
            } else {
                work.pop_back();
                if (low[u] == dfn[u]) {      // u 是分量的根，彈到 u 為止
                    while (true) {
                        int w = stk.back();
                        stk.pop_back();
                        on[w] = 0;
                        comp[w] = nc;
                        if (w == u) break;
                    }
                    nc++;
                }
                if (!work.empty()) {         // 回到父節點，傳遞 low
                    int p = work.back().first;
                    low[p] = min(low[p], low[u]);
                }
            }
        }
    }
    return comp;
}

// 把分量編號換成「該分量中最小點號」——與編號方式無關的正規形式。
static vector<int> compCanonical(const vector<int> &comp) {
    int mx = -1;
    for (int c : comp) mx = max(mx, c);
    vector<int> smallest(mx + 1, -1);
    for (int v = 0; v < (int)comp.size(); v++) {
        int c = comp[v];
        if (smallest[c] == -1 || v < smallest[c]) smallest[c] = v;
    }
    vector<int> res(comp.size());
    for (int v = 0; v < (int)comp.size(); v++) res[v] = smallest[comp[v]];
    return res;
}

// ---------------------------------------------------------------- 2-SAT

// 文字 → 點號：+(i+1) → 2i（真），-(i+1) → 2i+1（假）。
static int litNode(int lit) { return 2 * (abs(lit) - 1) + (lit > 0 ? 0 : 1); }

// 建蘊含圖：(a ∨ b) ⟹ 邊 ¬a→b 與 ¬b→a。
static vector<vector<int>> buildImplication(int n, const vector<Clause> &clauses) {
    vector<vector<int>> adj(2 * n);
    for (const Clause &cl : clauses) {
        int na = litNode(cl.first), nb = litNode(cl.second);
        adj[na ^ 1].push_back(nb);
        adj[nb ^ 1].push_back(na);
    }
    return adj;
}

// 縮點 DAG 的拓撲序排名：以「分量的最小點號」為鍵，鍵唯一 → 與 SCC 編號無關。
static vector<int> condensationRank(int n, const vector<vector<int>> &adj,
                                   const vector<int> &comp, int ncomp) {
    vector<int> cmin(ncomp, n);
    for (int v = 0; v < (int)comp.size(); v++) cmin[comp[v]] = min(cmin[comp[v]], v);

    vector<set<int>> cedge(ncomp);
    vector<int> indeg(ncomp, 0);
    for (int u = 0; u < (int)comp.size(); u++)
        for (int v : adj[u]) {
            int cu = comp[u], cv = comp[v];
            if (cu != cv && !cedge[cu].count(cv)) {
                cedge[cu].insert(cv);
                indeg[cv]++;
            }
        }

    // min-heap，鍵為 (分量最小點號, 分量編號)，與 Python 的 heapq 完全一致
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    for (int c = 0; c < ncomp; c++)
        if (indeg[c] == 0) pq.push({cmin[c], c});

    vector<int> rank(ncomp, 0);
    int t = 0;
    while (!pq.empty()) {
        int c = pq.top().second;
        pq.pop();
        rank[c] = t++;
        for (int d : cedge[c])
            if (--indeg[d] == 0) pq.push({cmin[d], d});
    }
    assert(t == ncomp);                        // 縮點後必為 DAG
    return rank;
}

// 解 2-SAT，返回 (是否有解, 一組解)。無解時解為空。
static pair<bool, vector<int>> solve2SAT(int n, const vector<Clause> &clauses) {
    vector<vector<int>> adj = buildImplication(n, clauses);
    vector<int> comp = kosarajuSCC(2 * n, adj);
    for (int v = 0; v < n; v++)
        if (comp[2 * v] == comp[2 * v + 1]) return {false, {}};   // x 與 ¬x 互通 → 矛盾
    int ncomp = comp.empty() ? 0 : *max_element(comp.begin(), comp.end()) + 1;
    vector<int> rank = condensationRank(2 * n, adj, comp, ncomp);
    // x 為真 ⟺ 真節點所在分量的拓撲序比較晚
    vector<int> assign(n);
    for (int v = 0; v < n; v++) assign[v] = rank[comp[2 * v]] > rank[comp[2 * v + 1]] ? 1 : 0;
    return {true, assign};
}

// ---------------------------------------------------------------- 對拍用暴力

// 文字在給定取值下是否為真。
static bool litValue(int lit, const vector<int> &assign) {
    int v = abs(lit) - 1;
    return assign[v] == (lit > 0 ? 1 : 0);
}

// 驗證一組取值是否滿足所有子句。
static bool checkAssign(int n, const vector<Clause> &clauses, const vector<int> &assign) {
    for (const Clause &cl : clauses)
        if (!(litValue(cl.first, assign) || litValue(cl.second, assign))) return false;
    return true;
}

// 暴力枚舉 2ⁿ 種取值，僅供對拍（n ≤ 12）。
static bool bf2SAT(int n, const vector<Clause> &clauses) {
    for (int mask = 0; mask < (1 << n); mask++) {
        vector<int> assign(n);
        for (int i = 0; i < n; i++) assign[i] = (mask >> i) & 1;
        if (checkAssign(n, clauses, assign)) return true;
    }
    return false;
}

// ---------------------------------------------------------------- 應用：兩時段排程

// 兩個半開區間 [s,e) 是否重疊（長度 0 的區間不與任何區間重疊）。
static bool overlap(pair<int, int> p, pair<int, int> q) {
    return max(p.first, q.first) < min(p.second, q.second);
}

// slots[i] = (a1, b1, a2, b2)：任務 i 可選 [a1,b1) 或 [a2,b2)，任務之間不得重疊。
// 返回 (是否可排, choice)，choice[i] = 0 選第一個、1 選第二個。
static pair<bool, vector<int>> scheduleTwoSlots(const vector<array<int, 4>> &slots) {
    int n = (int)slots.size();
    vector<Clause> clauses;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            pair<int, int> ai1 = {slots[i][0], slots[i][1]};
            pair<int, int> ai2 = {slots[i][2], slots[i][3]};
            pair<int, int> aj1 = {slots[j][0], slots[j][1]};
            pair<int, int> aj2 = {slots[j][2], slots[j][3]};
            if (overlap(ai1, aj1)) clauses.push_back({-(i + 1), -(j + 1)});
            if (overlap(ai1, aj2)) clauses.push_back({-(i + 1), (j + 1)});
            if (overlap(ai2, aj1)) clauses.push_back({(i + 1), -(j + 1)});
            if (overlap(ai2, aj2)) clauses.push_back({(i + 1), (j + 1)});
        }
    auto res = solve2SAT(n, clauses);
    if (!res.first) return {false, {}};
    vector<int> choice(n);
    for (int i = 0; i < n; i++) choice[i] = res.second[i] == 1 ? 0 : 1;
    return {true, choice};
}

// 暴力枚舉所有選擇，僅供對拍。
static pair<bool, vector<int>> bfSchedule(const vector<array<int, 4>> &slots) {
    int n = (int)slots.size();
    for (int mask = 0; mask < (1 << n); mask++) {
        vector<pair<int, int>> chosen(n);
        for (int i = 0; i < n; i++)
            chosen[i] = ((mask >> i) & 1) == 0 ? pair<int, int>{slots[i][0], slots[i][1]}
                                               : pair<int, int>{slots[i][2], slots[i][3]};
        bool ok = true;
        for (int i = 0; i < n && ok; i++)
            for (int j = i + 1; j < n; j++)
                if (overlap(chosen[i], chosen[j])) { ok = false; break; }
        if (ok) {
            vector<int> c(n);
            for (int i = 0; i < n; i++) c[i] = (mask >> i) & 1;
            return {true, c};
        }
    }
    return {false, {}};
}

// ---------------------------------------------------------------- 輸出

static string run(int n, const vector<Clause> &clauses) {
    auto res = solve2SAT(n, clauses);
    vector<vector<int>> adj = buildImplication(n, clauses);
    vector<int> comp = kosarajuSCC(2 * n, adj);
    int ncomp = comp.empty() ? 0 : *max_element(comp.begin(), comp.end()) + 1;
    vector<int> canon = compCanonical(comp);

    ostringstream oss;
    oss << (res.first ? "SATISFIABLE" : "UNSATISFIABLE") << "\n";
    if (res.first)
        for (int i = 0; i < n; i++) {
            if (i) oss << " ";
            oss << res.second[i];
        }
    oss << "\n" << ncomp << "\n";
    for (int i = 0; i < (int)canon.size(); i++) {
        if (i) oss << " ";
        oss << canon[i];
    }
    return oss.str();
}

// ---------------------------------------------------------------- 測試

static void check2SAT(int n, const vector<Clause> &clauses) {
    vector<vector<int>> adj = buildImplication(n, clauses);
    vector<int> c1 = kosarajuSCC(2 * n, adj);
    vector<int> c2 = tarjanSCC(2 * n, adj);
    assert(compCanonical(c1) == compCanonical(c2));
    int n1 = c1.empty() ? 0 : *max_element(c1.begin(), c1.end()) + 1;
    int n2 = c2.empty() ? 0 : *max_element(c2.begin(), c2.end()) + 1;
    assert(n1 == n2);

    auto res = solve2SAT(n, clauses);
    assert(res.first == bf2SAT(n, clauses));
    if (res.first) {
        assert((int)res.second.size() == n);
        assert(checkAssign(n, clauses, res.second));
    }
}

static void runTests() {
    mt19937 rng(20261012);

    // ---- 固定用例 ----
    check2SAT(2, {{1, 2}, {-1, 2}, {1, -2}});
    assert((solve2SAT(2, {{1, 2}, {-1, 2}, {1, -2}}) == pair<bool, vector<int>>{true, {1, 1}}));

    check2SAT(2, {{1, 2}, {-1, 2}, {1, -2}, {-1, -2}});
    assert(!solve2SAT(2, {{1, 2}, {-1, 2}, {1, -2}, {-1, -2}}).first);

    check2SAT(0, {});
    check2SAT(3, {});
    assert(run(0, {}) == "SATISFIABLE\n\n0\n");

    check2SAT(1, {{1, 1}});
    assert((solve2SAT(1, {{1, 1}}) == pair<bool, vector<int>>{true, {1}}));
    check2SAT(1, {{1, 1}, {-1, -1}});
    assert(!solve2SAT(1, {{1, 1}, {-1, -1}}).first);

    assert(run(2, {{1, 2}, {-1, 2}, {1, -2}}) == "SATISFIABLE\n1 1\n2\n0 1 0 1");
    assert(run(1, {{1, 1}}) == "SATISFIABLE\n1\n2\n0 1");

    // ---- 隨機對拍：小規模全枚舉 ----
    for (int iter = 0; iter < 600; iter++) {
        int n = 1 + (int)(rng() % 7);
        int m = (int)(rng() % 11);
        vector<Clause> clauses;
        for (int i = 0; i < m; i++) {
            int a = (rng() % 2 ? 1 : -1) * (1 + (int)(rng() % n));
            int b = (rng() % 2 ? 1 : -1) * (1 + (int)(rng() % n));
            clauses.push_back({a, b});
        }
        check2SAT(n, clauses);
    }

    // ---- 隨機對拍：子句數少，讓大 n 的有解/無解都出現 ----
    for (int iter = 0; iter < 300; iter++) {
        int n = 8 + (int)(rng() % 5);
        int m = (int)(rng() % 7);
        vector<Clause> clauses;
        for (int i = 0; i < m; i++) {
            int a = (rng() % 2 ? 1 : -1) * (1 + (int)(rng() % n));
            int b = (rng() % 2 ? 1 : -1) * (1 + (int)(rng() % n));
            clauses.push_back({a, b});
        }
        auto res = solve2SAT(n, clauses);
        assert(res.first == bf2SAT(n, clauses));
        if (res.first) assert(checkAssign(n, clauses, res.second));
    }

    // ---- 隨機圖：驗證兩種 SCC 寫法一致 ----
    for (int iter = 0; iter < 300; iter++) {
        int v = 1 + (int)(rng() % 20);
        int e = (int)(rng() % 41);
        vector<vector<int>> adj(v);
        for (int i = 0; i < e; i++) {
            int a = (int)(rng() % v), b = (int)(rng() % v);
            adj[a].push_back(b);
        }
        assert(compCanonical(kosarajuSCC(v, adj)) == compCanonical(tarjanSCC(v, adj)));
    }

    // ---- 應用：兩時段排程對拍 ----
    for (int iter = 0; iter < 400; iter++) {
        int k = (int)(rng() % 9);
        vector<array<int, 4>> slots;
        for (int i = 0; i < k; i++) {
            int a1 = (int)(rng() % 9), b1 = a1 + (int)(rng() % 4);
            int a2 = (int)(rng() % 9), b2 = a2 + (int)(rng() % 4);
            slots.push_back({{a1, b1, a2, b2}});
        }
        auto r1 = scheduleTwoSlots(slots);
        auto r2 = bfSchedule(slots);
        assert(r1.first == r2.first);
        if (r1.first) {
            vector<pair<int, int>> chosen(k);
            for (int i = 0; i < k; i++)
                chosen[i] = r1.second[i] == 0 ? pair<int, int>{slots[i][0], slots[i][1]}
                                              : pair<int, int>{slots[i][2], slots[i][3]};
            for (int i = 0; i < k; i++)
                for (int j = i + 1; j < k; j++) assert(!overlap(chosen[i], chosen[j]));
        }
    }

    // 經典不可排範例：三個任務，兩兩時段組合都撞
    vector<array<int, 4>> s3;
    s3.push_back({{0, 2, 0, 2}});
    s3.push_back({{0, 2, 0, 2}});
    s3.push_back({{0, 2, 0, 2}});
    assert(!scheduleTwoSlots(s3).first);
    // 可排：錯開即可
    vector<array<int, 4>> s2;
    s2.push_back({{0, 1, 5, 6}});
    s2.push_back({{1, 2, 6, 7}});
    assert(scheduleTwoSlots(s2).first);
}

// stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。
static bool isInt(const string &s) {
    if (s.empty()) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= (int)s.size()) return false;
    for (int j = i; j < (int)s.size(); j++)
        if (!isdigit((unsigned char)s[j])) return false;
    return true;
}

static void runIO(const vector<string> &tokens) {
    auto readInt = [&](size_t i) -> pair<bool, int> {
        if (i >= tokens.size() || !isInt(tokens[i])) return {false, 0};
        return {true, stoi(tokens[i])};
    };
    auto p0 = readInt(0);
    auto p1 = readInt(1);
    int n = p0.first ? max(0, p0.second) : 0;
    int m = p1.first ? max(0, p1.second) : 0;

    vector<Clause> clauses;
    for (int t = 0; t < m; t++) {
        auto a = readInt(2 + 2 * (size_t)t);
        auto b = readInt(3 + 2 * (size_t)t);
        if (!a.first || !b.first) break;
        if (a.second == 0 || b.second == 0 || abs(a.second) > n || abs(b.second) > n)
            continue;                            // 非法文字：整條子句忽略
        clauses.push_back({a.second, b.second});
    }
    cout << run(n, clauses) << "\n";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<string> tokens;
    bool any = false;
    string cur;
    int c;
    while ((c = cin.get()) != EOF) {
        if (isspace(c)) {
            if (!cur.empty()) { tokens.push_back(cur); cur.clear(); }
        } else {
            cur.push_back((char)c);
            any = true;
        }
    }
    if (!cur.empty()) tokens.push_back(cur);

    if (any) {
        runIO(tokens);
    } else {
        runTests();
        cout << "all tests passed" << endl;
    }
    return 0;
}
