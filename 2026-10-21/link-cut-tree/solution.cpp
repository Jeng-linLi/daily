// Link-Cut Tree（LCT）：用 splay 維護「動態森林」，支援加邊 / 刪邊 / 路徑查詢 / 路徑修改
//
// 題意：
//     有 n 個節點（編號 1..n），每個節點有一個整數權值。一開始整張圖**沒有邊**。
//     依序處理 q 個操作，操作會動態地連邊、斷邊（整個過程圖始終保持是森林），
//     並在任意時刻查詢／修改「u 到 v 這條路徑」上的資訊。
//
//     操作碼（每行一個操作）：
//         1  u v        連邊（若 u、v 已連通則忽略，避免產生環）
//         2  u v        斷邊（若邊 (u,v) 不存在則忽略）
//         3  u v        路徑查詢，輸出「權值和 最大值 最小值 邊數」；不連通輸出 -1
//         4  u v w      路徑修改：u 到 v 路徑上每個點權都加 w（不連通則忽略）
//         5  u v        連通性判定，輸出 1 / 0
//         6  u x        單點改值：把 u 的權值設成 x
//         7  u v        輸出路徑上的點權（由 u 走到 v 的順序）；不連通輸出 -1
//
// 思路：
//     ### 1. 為什麼需要 LCT
//     靜態樹上做路徑查詢，樹鏈剖分（heavy-light）就夠了；但它**不能加邊刪邊**。
//     LCT 用「**實鏈剖分**」解決這件事：把每個點的兒子中**至多一個**標成「實邊」，
//     其餘是「虛邊」。實邊串成的鏈用一棵 **splay** 維護（鍵是中序對應的樹上深度），
//     虛邊只記錄「這棵 splay 的根的父親是誰」（path-parent，單向）。
//
//     於是「森林」被拆成一棵棵 splay。兩個關鍵操作：
//
//     - **`access(x)`**：把根到 x 的路徑整條變成**實鏈**。
//       反覆「splay(x) → 把 x 的右子樹換成上一輪處理過的鏈 → 往上跳 path-parent」。
//       做完 x 所在的那棵 splay 就恰好含根到 x 的全部節點，且 x 是 splay 根。
//     - **`makeroot(x)`**：`access(x)` 之後翻轉整棵 splay（懶標記 rev），
//       把 x 變成整棵樹的**代表根**。這正是 LCT 能做到「無根樹任意路徑」的原因。
//
//     有了 `split(u,v) = makeroot(u); access(v)`，v 的那棵 splay 就**恰好**是 u..v 的路徑，
//     路徑查詢／修改全部退化成「對一個 splay 根做一次 O(1) 讀取 / O(log n) 懶標記」。
//
//     ### 2. 各操作的複雜度
//     全部都是 **均攤 O(log n)**（splay 的攤還分析）。這是目前已知動態樹的主流最優解。
//
//     ### 3. 兩語言必須一致的五個細節
//     - `access` 的寫法統一用 `for y in ...` 迴圈，且**不要**額外維護 `fa[last]`：
//       `x = fa[y = x]` 這個賦值順序本身就保證了 `fa[y] == x`，多寫一次反而可能寫錯。
//     - `isRoot(x)` 的定義：`fa[x] == 0` 或 x 不是 `fa[x]` 的**任何一個**兒子。
//       （splay 的父親 vs 虛邊的 path-parent 靠這個區分。）
//     - `cut(u,v)` 的成立條件是 `ch[v][0] == u && ch[u][1] == 0`：
//       前者表示 u 是 v 的前驅，後者表示 u..v 之間沒有第三個點 —— 兩者合起來才是「直接相鄰」。
//     - `rev` 與 `add` 兩個懶標記的**下推順序**固定為「先 rev 再 add」。
//     - 路徑輸出用 splay 的**中序遍歷**：`makeroot(u); access(v)` 之後中序就是 u → v。
//       遍歷時每一步都要先 `pushdown` 再往下走。
//
//     ### 4. 容易踩的坑
//     - 哨兵節點 0 的 `mx[0]` 必須是 -INF、`mn[0]` 必須是 +INF，否則取 max/min 會被 0 污染。
//     - `pushadd` 要乘的是 **splay 的子樹大小**（不是鏈長），所以 `sz` 必須維護。
//     - `findroot` 要先 `access(x)`，往左下角走到底，中途每一層都要 `pushdown`，
//       否則 rev 標記沒解開會走錯方向。
//
// 應用場景：
//     - 動態連通性：網路拓撲變動（連線／斷線）下的即時連通查詢
//     - 線上維護生成樹：邊權會改、邊會加刪的最小生成樹（動態 MST）
//     - 交通路網：施工封路（cut）與新路開通（link）後的兩地最長／最短瓶頸路段
//     - 版本管理：分支合併圖上查詢兩版本間的修改總量
//     - 遊戲伺服器：動態變化的地圖傳送網路上的路徑統計
//
// 複雜度（n 個點，q 次操作）：
//     link / cut / 路徑查詢 / 路徑修改     均攤 O(log n)
//     單點改值                             均攤 O(log n)
//     空間                                 O(n)

#include <algorithm>
#include <array>
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;

static const long long INFLL = (1LL << 60);
static const int NMAX = 300;
static const int QMAX = 400;
static const long long VMAX = 1000000LL;

// ------------------------------------------------------------------ LCT

class LCT {
public:
    int n;
    vector<array<int, 2>> ch;
    vector<int> fa, rev, sz;
    vector<long long> lazy, val, sm, mx, mn;

    explicit LCT(const vector<long long>& vals) {
        n = (int)vals.size();
        ch.assign(n + 1, {0, 0});
        fa.assign(n + 1, 0);
        rev.assign(n + 1, 0);
        sz.assign(n + 1, 0);
        lazy.assign(n + 1, 0);
        val.assign(n + 1, 0);
        sm.assign(n + 1, 0);
        mx.assign(n + 1, 0);
        mn.assign(n + 1, 0);
        mx[0] = -INFLL;
        mn[0] = INFLL;
        for (int i = 1; i <= n; ++i) {
            long long v = vals[i - 1];
            val[i] = v;
            sm[i] = v;
            mx[i] = v;
            mn[i] = v;
            sz[i] = 1;
        }
    }

    void pushup(int x) {
        int l = ch[x][0], r = ch[x][1];
        sz[x] = 1 + sz[l] + sz[r];
        sm[x] = val[x] + sm[l] + sm[r];
        long long m = val[x];
        if (mx[l] > m) m = mx[l];
        if (mx[r] > m) m = mx[r];
        mx[x] = m;
        m = val[x];
        if (mn[l] < m) m = mn[l];
        if (mn[r] < m) m = mn[r];
        mn[x] = m;
    }

    void pushrev(int x) {
        if (x == 0) return;
        swap(ch[x][0], ch[x][1]);
        rev[x] ^= 1;
    }

    void pushadd(int x, long long w) {
        if (x == 0) return;
        val[x] += w;
        sm[x] += w * (long long)sz[x];
        mx[x] += w;
        mn[x] += w;
        lazy[x] += w;
    }

    void pushdown(int x) {
        if (rev[x]) {
            pushrev(ch[x][0]);
            pushrev(ch[x][1]);
            rev[x] = 0;
        }
        if (lazy[x] != 0) {
            pushadd(ch[x][0], lazy[x]);
            pushadd(ch[x][1], lazy[x]);
            lazy[x] = 0;
        }
    }

    bool isRoot(int x) {
        int f = fa[x];
        if (f == 0) return true;
        return ch[f][0] != x && ch[f][1] != x;
    }

    void rotate(int x) {
        int y = fa[x], z = fa[y];
        int k = (ch[y][1] == x) ? 1 : 0;
        if (!isRoot(y)) {
            if (ch[z][1] == y) ch[z][1] = x;
            else ch[z][0] = x;
        }
        fa[x] = z;
        ch[y][k] = ch[x][k ^ 1];
        if (ch[x][k ^ 1]) fa[ch[x][k ^ 1]] = y;
        ch[x][k ^ 1] = y;
        fa[y] = x;
        pushup(y);
        pushup(x);
    }

    void splay(int x) {
        vector<int> st;
        int y = x;
        st.push_back(y);
        while (!isRoot(y)) {
            y = fa[y];
            st.push_back(y);
        }
        for (int i = (int)st.size() - 1; i >= 0; --i) pushdown(st[i]);
        while (!isRoot(x)) {
            int y = fa[x], z = fa[y];
            if (!isRoot(y)) {
                bool a = (ch[y][0] == x);
                bool b = (ch[z][0] == y);
                if (a != b) rotate(x);
                else rotate(y);
            }
            rotate(x);
        }
        pushup(x);
    }

    void access(int x) {
        int last = 0;
        int y = x;
        while (y != 0) {
            splay(y);
            ch[y][1] = last;
            pushup(y);
            last = y;
            y = fa[y];
        }
        splay(x);
    }

    void makeroot(int x) {
        access(x);
        pushrev(x);
    }

    int findroot(int x) {
        access(x);
        int y = x;
        while (ch[y][0] != 0) {
            pushdown(y);
            y = ch[y][0];
        }
        splay(y);
        return y;
    }

    void split(int u, int v) {
        makeroot(u);
        access(v);
    }

    bool connected(int u, int v) {
        if (u == v) return true;
        makeroot(u);
        return findroot(v) == u;
    }

    void link(int u, int v) {
        if (u == v) return;
        makeroot(u);
        if (findroot(v) == u) return;   // 已連通，連了會成環
        fa[u] = v;
    }

    void cut(int u, int v) {
        if (u == v) return;
        makeroot(u);
        access(v);
        // v 已是 splay 根；u 必須緊貼在 v 左邊、且 u 右邊沒有別的點
        if (ch[v][0] == u && ch[u][1] == 0) {
            ch[v][0] = 0;
            fa[u] = 0;
            pushup(v);
        }
    }

    bool pathQuery(int u, int v, vector<long long>& res) {
        if (!connected(u, v)) return false;
        split(u, v);
        res.clear();
        res.push_back(sm[v]);
        res.push_back(mx[v]);
        res.push_back(mn[v]);
        res.push_back((long long)sz[v] - 1);
        return true;
    }

    void pathAdd(int u, int v, long long w) {
        if (!connected(u, v)) return;
        split(u, v);
        pushadd(v, w);
    }

    void setValue(int u, long long x) {
        access(u);          // u 變成 splay 根，直接改權值即可
        val[u] = x;
        pushup(u);
    }

    bool pathValues(int u, int v, vector<long long>& res) {
        if (!connected(u, v)) return false;
        split(u, v);
        res.clear();
        vector<int> st;
        int cur = v;
        while (!st.empty() || cur != 0) {
            while (cur != 0) {
                pushdown(cur);
                st.push_back(cur);
                cur = ch[cur][0];
            }
            cur = st.back();
            st.pop_back();
            res.push_back(val[cur]);
            cur = ch[cur][1];
        }
        return true;
    }
};

// ------------------------------------------------------- 暴力對照（一般圖森林）

class Brute {
public:
    vector<long long> val;
    vector<set<int>> adj;

    explicit Brute(const vector<long long>& vals) {
        val.assign(1, 0);
        for (long long v : vals) val.push_back(v);
        adj.assign((int)vals.size() + 1, set<int>());
    }

    vector<int> path(int u, int v) {
        vector<int> empty;
        if (u == v) {
            empty.push_back(u);
            return empty;
        }
        int n = (int)adj.size() - 1;
        vector<int> par(n + 1, 0);
        par[u] = -1;
        deque<int> q;
        q.push_back(u);
        while (!q.empty()) {
            int x = q.front();
            q.pop_front();
            for (int y : adj[x]) {
                if (par[y] == 0) {
                    par[y] = x;
                    if (y == v) {
                        vector<int> p;
                        p.push_back(v);
                        while (p.back() != u) p.push_back(par[p.back()]);
                        reverse(p.begin(), p.end());
                        return p;
                    }
                    q.push_back(y);
                }
            }
        }
        return empty;    // 空集合代表不連通
    }

    bool connected(int u, int v) {
        return !path(u, v).empty();
    }

    void link(int u, int v) {
        if (connected(u, v)) return;
        adj[u].insert(v);
        adj[v].insert(u);
    }

    void cutEdge(int u, int v) {
        adj[u].erase(v);
        adj[v].erase(u);
    }

    bool pathQuery(int u, int v, vector<long long>& res) {
        vector<int> p = path(u, v);
        if (p.empty()) return false;
        long long s = 0, mxv = val[p[0]], mnv = val[p[0]];
        for (int i : p) {
            s += val[i];
            if (val[i] > mxv) mxv = val[i];
            if (val[i] < mnv) mnv = val[i];
        }
        res.clear();
        res.push_back(s);
        res.push_back(mxv);
        res.push_back(mnv);
        res.push_back((long long)p.size() - 1);
        return true;
    }

    void pathAdd(int u, int v, long long w) {
        vector<int> p = path(u, v);
        if (p.empty()) return;
        for (int i : p) val[i] += w;
    }

    void setValue(int u, long long x) { val[u] = x; }

    bool pathValues(int u, int v, vector<long long>& res) {
        vector<int> p = path(u, v);
        if (p.empty()) return false;
        res.clear();
        for (int i : p) res.push_back(val[i]);
        return true;
    }
};

// ------------------------------------------------------------------ 驅動

static string run(const vector<long long>& vals, const vector<vector<long long>>& ops) {
    int n = (int)vals.size();
    LCT lct(vals);
    Brute bf(vals);
    ostringstream oss;

    auto ok = [&](long long u) -> bool { return u >= 1 && u <= n; };

    for (const auto& op : ops) {
        if (op.empty()) continue;
        long long code = op[0];
        if (code == 1) {
            int u = (int)op[1], v = (int)op[2];
            if (ok(u) && ok(v)) {
                lct.link(u, v);
                bf.link(u, v);
            }
        } else if (code == 2) {
            int u = (int)op[1], v = (int)op[2];
            if (ok(u) && ok(v)) {
                lct.cut(u, v);
                bf.cutEdge(u, v);
            }
        } else if (code == 3) {
            int u = (int)op[1], v = (int)op[2];
            if (ok(u) && ok(v)) {
                vector<long long> r, e;
                bool br = lct.pathQuery(u, v, r);
                bool be = bf.pathQuery(u, v, e);
                assert(br == be);
                if (br) assert(r == e);
                if (!br) oss << "-1\n";
                else oss << r[0] << " " << r[1] << " " << r[2] << " " << r[3] << "\n";
            }
        } else if (code == 4) {
            int u = (int)op[1], v = (int)op[2];
            if (ok(u) && ok(v)) {
                lct.pathAdd(u, v, op[3]);
                bf.pathAdd(u, v, op[3]);
            }
        } else if (code == 5) {
            int u = (int)op[1], v = (int)op[2];
            if (ok(u) && ok(v)) {
                bool r = lct.connected(u, v);
                bool e = bf.connected(u, v);
                assert(r == e);
                oss << (r ? "1\n" : "0\n");
            }
        } else if (code == 6) {
            int u = (int)op[1];
            if (ok(u)) {
                lct.setValue(u, op[2]);
                bf.setValue(u, op[2]);
            }
        } else if (code == 7) {
            int u = (int)op[1], v = (int)op[2];
            if (ok(u) && ok(v)) {
                vector<long long> r, e;
                bool br = lct.pathValues(u, v, r);
                bool be = bf.pathValues(u, v, e);
                assert(br == be);
                if (br) assert(r == e);
                if (!br) {
                    oss << "-1\n";
                } else {
                    oss << r.size();
                    for (long long x : r) oss << " " << x;
                    oss << "\n";
                }
            }
        }
    }
    return oss.str();
}

static void checkCase(const vector<long long>& vals, const vector<vector<long long>>& ops) {
    run(vals, ops);
}

static void runTests() {
    checkCase({1, 2, 3, 4}, {
        {1, 1, 2}, {1, 2, 3}, {1, 3, 4},
        {3, 1, 4}, {7, 1, 4}, {5, 1, 4}, {5, 1, 3},
        {4, 2, 3, 10}, {3, 1, 4},
        {2, 2, 3}, {5, 1, 4}, {3, 1, 4},
    });

    checkCase({5, -3, 7, 0, 2}, {
        {1, 1, 2}, {1, 3, 4},
        {5, 1, 3}, {3, 3, 4},
        {1, 2, 3}, {5, 1, 5}, {3, 1, 5},
        {1, 4, 5}, {3, 1, 5}, {7, 5, 1},
        {6, 5, 100}, {3, 1, 5},
        {2, 3, 4}, {5, 3, 4},
    });

    checkCase({42}, {
        {3, 1, 1}, {7, 1, 1}, {5, 1, 1},
        {1, 1, 1}, {5, 1, 1},
        {2, 1, 1}, {3, 1, 2}, {5, 0, 1}, {3, 1, 0},
    });

    checkCase({}, {{3, 1, 1}, {5, 1, 2}});

    mt19937 rng(20261021);
    auto ri = [&](int a, int b) -> int {
        return (int)(rng() % (unsigned)(b - a + 1)) + a;
    };
    for (int t = 0; t < 220; ++t) {
        int n = ri(1, 16);
        vector<long long> vals(n);
        for (int i = 0; i < n; ++i) vals[i] = ri(-40, 40);
        vector<vector<long long>> ops;
        int m = ri(1, 40);
        for (int i = 0; i < m; ++i) {
            int code = ri(1, 9);
            if (code == 9) code = 6;
            else if (code == 8) code = 7;
            int u = ri(1, n), v = ri(1, n);
            if (code == 4) ops.push_back({4, u, v, ri(-20, 20)});
            else if (code == 6) ops.push_back({6, u, ri(-50, 50)});
            else ops.push_back({code, u, v});
        }
        checkCase(vals, ops);
    }

    for (int t = 0; t < 40; ++t) {
        int n = ri(2, 30);
        vector<long long> vals(n);
        for (int i = 0; i < n; ++i) vals[i] = ri(-100, 100);
        vector<vector<long long>> ops;
        int m = ri(1, 60);
        for (int i = 0; i < m; ++i) {
            int code = ri(1, 7);
            int u = ri(1, n), v = ri(1, n);
            if (code == 4) ops.push_back({4, u, v, ri(-30, 30)});
            else ops.push_back({code, u, v});
        }
        checkCase(vals, ops);
    }
}

// -- IO 解析：嚴格規則，非數字 token 一律視為「輸入到此為止」 --

static bool isIntTok(const string& s) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= s.size()) return false;
    for (size_t k = i; k < s.size(); ++k) {
        if (s[k] < '0' || s[k] > '9') return false;
    }
    return true;
}

static long long clampBig(long long v) {
    const long long LIM = 1000000000000000000LL;
    if (v > LIM) v = LIM;
    if (v < -LIM) v = -LIM;
    return v;
}

static long long clampV(long long v) {
    if (v > VMAX) v = VMAX;
    if (v < -VMAX) v = -VMAX;
    return v;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    string raw, line;
    ostringstream buf;
    while (getline(cin, line)) {
        buf << line << "\n";
    }
    raw = buf.str();

    bool hasInput = false;
    for (char c : raw) {
        if (!isspace((unsigned char)c)) { hasInput = true; break; }
    }

    if (!hasInput) {
        runTests();
        fputs("all tests passed\n", stdout);
        return 0;
    }

    vector<string> toks;
    {
        string cur;
        for (char c : raw) {
            if (isspace((unsigned char)c)) {
                if (!cur.empty()) { toks.push_back(cur); cur.clear(); }
            } else {
                cur.push_back(c);
            }
        }
        if (!cur.empty()) toks.push_back(cur);
    }

    size_t pos = 0;
    auto nxt = [&](long long defVal) -> long long {
        if (pos >= toks.size()) return defVal;
        string t = toks[pos++];
        if (!isIntTok(t)) return defVal;
        errno = 0;
        char* endp = nullptr;
        long long v = strtoll(t.c_str(), &endp, 10);
        if (endp == t.c_str()) return defVal;
        return clampBig(v);
    };

    long long nl = nxt(0);
    if (nl < 0) nl = 0;
    if (nl > NMAX) nl = NMAX;
    int n = (int)nl;
    vector<long long> vals(n);
    for (int i = 0; i < n; ++i) vals[i] = clampV(nxt(0));

    long long ql = nxt(0);
    if (ql < 0) ql = 0;
    if (ql > QMAX) ql = QMAX;

    vector<vector<long long>> ops;
    for (long long i = 0; i < ql; ++i) {
        long long code = nxt(0);
        if (code == 1 || code == 2 || code == 3 || code == 5 || code == 7) {
            ops.push_back({code, nxt(1), nxt(1)});
        } else if (code == 4) {
            long long u = nxt(1), v = nxt(1);
            ops.push_back({4, u, v, clampV(nxt(0))});
        } else if (code == 6) {
            ops.push_back({6, nxt(1), clampV(nxt(0))});
        }
    }

    string out = run(vals, ops);
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
