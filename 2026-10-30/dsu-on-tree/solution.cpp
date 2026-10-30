// 樹上啟發式合併（DSU on Tree / Sack）：O(n log n) 做子樹顏色統計
//
// 題意：
//     給定一棵有 n 個節點的根樹（根為 1），節點 u 有顏色 c[u]（1 <= c[u] <= n）。
//     有 q 個詢問，每個詢問給出一對 (u, f)，要求回答 u 的**子樹**中：
//         distinct —— 出現過的顏色數
//         eq       —— 出現次數**恰好等於 f** 的顏色數
//         shared   —— 出現過、而且在 u 的子樹**之外**也出現過的顏色數
//     限制：1 <= n, q <= 100000，0 <= f <= n。
//
// 思路：
//     ### 1. 為什麼需要這招
//     樸素做法是對每個節點 u 單獨遍歷一次子樹，O(n^2)。
//     「樹上啟發式合併」（又稱 DSU on Tree / Sack / 同步 bags）
//     藉由**重輕劃分**把總複雜度降到 O(n log n)。
//
//     ### 2. 重輕劃分（heavy-light by subtree size）
//     對每個節點 u，取子樹最大的那個兒子當「重兒子」（heavy child），
//     其餘都是「輕兒子」（light child）。
//     關鍵性質：從根往下走，**輕邊最多只有 O(log n) 條** ——
//     因為走過一條輕邊 (u,v) 時 size[v] <= size[u]/2，
//     子樹大小每走一條輕邊就至少減半，所以輕邊數是 O(log n)。
//
//     ### 3. Sack 的流程
//     維護一個「全域狀態」（bag），記錄目前被加入的節點的顏色統計。
//     對節點 u 的處理分成三段：
//         (a) 先處理所有**輕兒子**，處理完立刻清空它們的貢獻（keep = false）
//         (b) 再處理**重兒子**，但**保留**它的貢獻（keep = true）
//         (c) 把所有輕兒子的子樹重新加回來，再加上 u 自己，
//             此刻 bag 恰好就是 u 的子樹 —— 據此回答 u 的詢問
//         (d) 若 keep = false，把 u 的整個子樹移出 bag
//     輕兒子的貢獻被加入又移除了兩次、重兒子只加入一次，
//     而每個節點只會是它祖先的 O(log n) 個輕兒子之一，
//     所以總加入次數是 O(n log n)。
//
//     不變式：**每次進入 dfs(u, ·) 時 bag 都是空的**，
//     而 dfs(u, true) 返回時 bag 恰好是 u 的子樹。
//     所以 (c) 段之後 bag = u 的子樹，可以直接回答 u 的詢問。
//
//     ### 4. 維護哪些統計量
//     cnt[c]   = 顏色 c 在 bag 中的出現次數
//     freq[f]  = 出現次數恰好為 f 的顏色數（f >= 0）
//     distinct = 出現過的顏色數
//     shared   = 出現過、且在子樹之外也出現過的顏色數
//     total[c] = 顏色 c 在整棵樹中的總出現次數（一次預處理）
//     palette  = 整棵樹用到的不同顏色數
//
//     加入一個顏色 c（cnt 從 a 變成 a+1）：
//         - freq[a] -= 1；freq[a+1] += 1
//         - distinct：a == 0 時 +1
//         - shared：判斷 (0 < cnt[c] < total[c]) 在變化前後是否成立，
//           成立狀態改變就相應 +1 / -1
//     移出一個顏色 c 完全對稱。於是每次加入 / 移出都是 O(1)。
//
//     顏色的「宇宙」取整棵樹用到的 palette 種顏色，
//     因此 freq[0] 初始值 = palette，eq(0) = freq[0] = palette - distinct
//     （也就是「在整棵樹出現過、但沒出現在這個子樹裡」的顏色數），
//     與 shared 的語意保持一致。
//
// 複雜度：
//     預處理        O(n)
//     單次詢問      O(1)（在 DFS 到 u 的當下直接算）
//     總計          O(n log n + q)，空間 O(n)
//
// 應用場景：
//     - 專案目錄的檔案副檔名統計（每個子目錄有多少種副檔名、
//       有多少副檔名只出現一次、哪些副檔名與其他目錄共用）
//     - 基因組 / 蛋白質序列的子樹特徵統計
//     - 程式碼倉庫的語言分布分析（每個資料夾用了幾種語言）
//     - 線上統計「某個資料夾內有多少種唯一的檔案類型」

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr int QMAX = 100000;
constexpr long long PARSE_LIMIT = 1000000000000000000LL;

struct Answer {
    int distinct = 0;
    int eq = 0;
    int shared = 0;

    bool operator==(const Answer& o) const {
        return distinct == o.distinct && eq == o.eq && shared == o.shared;
    }
};

class Solver {
public:
    int n = 0;
    std::vector<int> color;
    std::vector<int> parent;
    // u -> [(f, query_index), ...]
    std::map<int, std::vector<std::pair<int, int>>> queries;

    std::vector<std::vector<int>> children;
    std::vector<std::vector<int>> light;
    std::vector<int> order;   // DFS 先序
    std::vector<int> tin;     // 先序中的位置
    std::vector<int> sz;      // 子樹大小
    std::vector<int> heavy;   // 重兒子（0 表示沒有）

    std::vector<int> total;   // 顏色在整棵樹的總出現次數
    int palette = 0;

    // 全域狀態
    std::vector<int> cnt;
    std::vector<int> freq;
    int distinct = 0;
    int shared = 0;

    std::vector<Answer> answers;

    void build(int n_, const std::vector<int>& color_,
               const std::vector<int>& parent_,
               const std::map<int, std::vector<std::pair<int, int>>>& queries_) {
        n = n_;
        color = color_;
        parent = parent_;
        queries = queries_;

        children.assign(n + 1, std::vector<int>());
        for (int u = 2; u <= n; ++u) {
            children[parent[u]].push_back(u);
        }

        // ---- 迭代 DFS：取得先序序列與 tin ----
        // 先序序列中，u 的子樹恰好佔據連續區間 [tin[u], tin[u]+sz[u])
        order.clear();
        order.reserve(n);
        tin.assign(n + 1, 0);
        std::vector<int> st;
        st.push_back(1);
        while (!st.empty()) {
            int u = st.back();
            st.pop_back();
            tin[u] = static_cast<int>(order.size());
            order.push_back(u);
            for (int i = static_cast<int>(children[u].size()) - 1; i >= 0; --i) {
                st.push_back(children[u][i]);
            }
        }

        // ---- 子樹大小、重兒子、輕兒子 ----
        sz.assign(n + 1, 1);
        for (int i = n - 1; i >= 0; --i) {
            int u = order[i];
            for (int c : children[u]) {
                sz[u] += sz[c];
            }
        }
        heavy.assign(n + 1, 0);
        for (int u = 2; u <= n; ++u) {
            int p = parent[u];
            if (heavy[p] == 0 || sz[u] > sz[heavy[p]]) {
                heavy[p] = u;
            }
        }
        light.assign(n + 1, std::vector<int>());
        for (int u = 1; u <= n; ++u) {
            int h = heavy[u];
            for (int c : children[u]) {
                if (c != h) light[u].push_back(c);
            }
        }

        // ---- 整棵樹的顏色總數 ----
        total.assign(n + 1, 0);
        for (int u = 1; u <= n; ++u) {
            total[color[u]] += 1;
        }
        palette = 0;
        for (int c = 1; c <= n; ++c) {
            if (total[c] > 0) ++palette;
        }

        // ---- 全域狀態 ----
        cnt.assign(n + 1, 0);
        freq.assign(n + 2, 0);
        freq[0] = palette;   // 顏色的宇宙 = palette 種
        distinct = 0;
        shared = 0;

        int nq = 0;
        for (const auto& kv : queries) {
            nq += static_cast<int>(kv.second.size());
        }
        answers.assign(nq, Answer());
    }

    // ------------------------------------------------------------ 狀態操作

    void addColor(int c) {
        int a = cnt[c];
        int b = a + 1;
        bool before = (0 < a && a < total[c]);
        bool after = (0 < b && b < total[c]);
        if (after && !before) {
            ++shared;
        } else if (before && !after) {
            --shared;
        }
        cnt[c] = b;
        --freq[a];
        ++freq[b];
        if (a == 0) ++distinct;
    }

    void removeColor(int c) {
        int a = cnt[c];
        int b = a - 1;
        bool before = (0 < a && a < total[c]);
        bool after = (0 < b && b < total[c]);
        if (after && !before) {
            ++shared;
        } else if (before && !after) {
            --shared;
        }
        cnt[c] = b;
        --freq[a];
        ++freq[b];
        if (b == 0) --distinct;
    }

    // ------------------------------------------------------------ 主流程

    std::vector<Answer> run() {
        // 迭代版 Sack：frame = (u, keep, stage)
        //   stage 0 = 進入（安排子節點的處理順序）
        //   stage 1 = 子節點都處理完，輪到自己
        struct Frame {
            int u;
            bool keep;
            int stage;
        };
        std::vector<Frame> st;
        st.push_back({1, true, 0});
        while (!st.empty()) {
            Frame fr = st.back();
            st.pop_back();
            int u = fr.u;
            if (fr.stage == 0) {
                // 堆疊後進先出：push 順序要反過來，
                // pop 出來的順序才是「輕兒子 -> 重兒子 -> 自己收尾」
                st.push_back({u, fr.keep, 1});
                int h = heavy[u];
                if (h) st.push_back({h, true, 0});
                for (int i = static_cast<int>(light[u].size()) - 1; i >= 0; --i) {
                    st.push_back({light[u][i], false, 0});
                }
            } else {
                // (c) 把所有輕兒子的子樹加回來
                for (int c : light[u]) {
                    int start = tin[c];
                    for (int i = start; i < start + sz[c]; ++i) {
                        addColor(color[order[i]]);
                    }
                }
                // 再加上 u 自己 —— 此刻 bag 恰好是 u 的子樹
                addColor(color[u]);
                // 回答 u 的詢問
                auto it = queries.find(u);
                if (it != queries.end()) {
                    for (const auto& pr : it->second) {
                        int f = pr.first;
                        int qi = pr.second;
                        answers[qi] = Answer{distinct, freq[f], shared};
                    }
                }
                // (d) 若不需要保留，把 u 的整個子樹移出 bag
                if (!fr.keep) {
                    int start = tin[u];
                    for (int i = start; i < start + sz[u]; ++i) {
                        removeColor(color[order[i]]);
                    }
                }
            }
        }
        return answers;
    }
};

// 暴力參考實作：對每個詢問單獨遍歷該節點的子樹
std::vector<Answer> solveBrute(int n, const std::vector<int>& color,
                               const std::vector<int>& parent,
                               const std::vector<std::pair<int, int>>& qlist) {
    std::vector<std::vector<int>> children(n + 1);
    for (int u = 2; u <= n; ++u) {
        children[parent[u]].push_back(u);
    }
    std::vector<int> total(n + 1, 0);
    for (int u = 1; u <= n; ++u) {
        total[color[u]] += 1;
    }
    int palette = 0;
    for (int c = 1; c <= n; ++c) {
        if (total[c] > 0) ++palette;
    }

    std::vector<Answer> out;
    for (const auto& pr : qlist) {
        int u = pr.first;
        int f = pr.second;
        std::map<int, int> cnt;
        std::vector<int> st;
        st.push_back(u);
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            cnt[color[x]] += 1;
            for (int c : children[x]) {
                st.push_back(c);
            }
        }
        Answer a;
        a.distinct = static_cast<int>(cnt.size());
        if (f == 0) {
            a.eq = palette - a.distinct;
        } else {
            for (const auto& kv : cnt) {
                if (kv.second == f) ++a.eq;
            }
        }
        for (const auto& kv : cnt) {
            if (kv.second < total[kv.first]) ++a.shared;
        }
        out.push_back(a);
    }
    return out;
}

// 簡易 LCG，兩種語言共用同一組常數以確保亂數序列一致
struct Lcg {
    unsigned long long s;
    explicit Lcg(unsigned long long seed) : s(seed) {}
    unsigned long long next() {
        s = (s * 6364136223846793005ULL + 1442695040888963407ULL) &
            0xFFFFFFFFFFFFFFFFULL;
        return s;
    }
    int range(int lo, int hi) {  // [lo, hi]
        return lo + static_cast<int>(
                          next() % static_cast<unsigned long long>(hi - lo + 1));
    }
};

void fail(const char* msg) {
    std::cerr << msg << "\n";
    std::exit(1);
}

void runTests() {
    // ---- 固定用例 ----
    // 樹：1-2, 1-3, 2-4, 3-5, 3-6（根 = 1）
    // 顏色：c1=1, c2=2, c3=1, c4=3, c5=1, c6=1
    {
        int n = 6;
        std::vector<int> color = {0, 1, 2, 1, 3, 1, 1};
        std::vector<int> parent = {0, 0, 1, 1, 2, 3, 3};
        std::vector<std::pair<int, int>> qlist = {
            {1, 0}, {1, 1}, {1, 2}, {1, 4}, {1, 3}, {2, 1}, {3, 1},
            {4, 1}, {5, 1}, {2, 0}, {3, 0}, {1, 6}, {6, 1}};
        std::map<int, std::vector<std::pair<int, int>>> queries;
        for (int i = 0; i < static_cast<int>(qlist.size()); ++i) {
            queries[qlist[i].first].push_back({qlist[i].second, i});
        }
        Solver s;
        s.build(n, color, parent, queries);
        std::vector<Answer> got = s.run();
        std::vector<Answer> want = solveBrute(n, color, parent, qlist);
        if (got != want) fail("case A brute mismatch");
        // 手算期望值
        const Answer e0{3, 0, 0};  // (u=1, f=0)
        const Answer e1{3, 2, 0};  // (u=1, f=1)
        const Answer e2{3, 0, 0};  // (u=1, f=2)
        const Answer e3{3, 1, 0};  // (u=1, f=4)
        const Answer e4{3, 0, 0};  // (u=1, f=3)
        const Answer e5{2, 2, 0};  // (u=2, f=1)
        const Answer e6{1, 0, 1};  // (u=3, f=1)
        const Answer e7{1, 1, 0};  // (u=4, f=1)
        const Answer e8{1, 1, 1};  // (u=5, f=1)
        const Answer e9{2, 1, 0};  // (u=2, f=0)
        const Answer e10{1, 2, 1}; // (u=3, f=0)
        const Answer e11{3, 0, 0}; // (u=1, f=6)
        const Answer e12{1, 1, 1}; // (u=6, f=1)
        const Answer expect[] = {e0, e1, e2, e3, e4, e5, e6,
                                 e7, e8, e9, e10, e11, e12};
        for (int i = 0; i < 13; ++i) {
            if (!(got[i] == expect[i])) fail("case A expected mismatch");
        }
    }
    // 單節點樹
    {
        std::vector<int> color = {0, 1};
        std::vector<int> parent = {0, 0};
        std::vector<std::pair<int, int>> qlist = {{1, 0}, {1, 1}, {1, 2}};
        std::map<int, std::vector<std::pair<int, int>>> queries;
        for (int i = 0; i < 3; ++i) queries[qlist[i].first].push_back({qlist[i].second, i});
        Solver s;
        s.build(1, color, parent, queries);
        if (s.run() != solveBrute(1, color, parent, qlist)) {
            fail("case B failed");
        }
    }
    // 所有節點同色
    {
        int n = 5;
        std::vector<int> color = {0, 5, 5, 5, 5, 5};
        std::vector<int> parent = {0, 0, 1, 2, 3, 4};
        std::vector<std::pair<int, int>> qlist;
        for (int u = 1; u <= 5; ++u) {
            qlist.push_back({u, 0});
            qlist.push_back({u, 1});
            qlist.push_back({u, 5});
        }
        std::map<int, std::vector<std::pair<int, int>>> queries;
        for (int i = 0; i < static_cast<int>(qlist.size()); ++i) {
            queries[qlist[i].first].push_back({qlist[i].second, i});
        }
        Solver s;
        s.build(n, color, parent, queries);
        if (s.run() != solveBrute(n, color, parent, qlist)) fail("case C failed");
    }

    // ---- 隨機對拍 ----
    {
        Lcg rng(20261030ULL);
        for (int iter = 0; iter < 200; ++iter) {
            int n = rng.range(1, 12);
            std::vector<int> color(n + 1, 0);
            for (int u = 1; u <= n; ++u) color[u] = rng.range(1, n);
            std::vector<int> parent(n + 1, 0);
            for (int u = 2; u <= n; ++u) parent[u] = rng.range(1, u - 1);
            int q = rng.range(1, 15);
            std::vector<std::pair<int, int>> qlist;
            for (int i = 0; i < q; ++i) {
                qlist.push_back({rng.range(1, n), rng.range(0, n)});
            }
            std::map<int, std::vector<std::pair<int, int>>> queries;
            for (int i = 0; i < q; ++i) {
                queries[qlist[i].first].push_back({qlist[i].second, i});
            }
            Solver s;
            s.build(n, color, parent, queries);
            if (s.run() != solveBrute(n, color, parent, qlist)) {
                fail("random mismatch");
            }
        }
    }
    // 稍大規模：隨機樹 + 每個節點都問
    {
        Lcg rng(13572468ULL);
        for (int iter = 0; iter < 12; ++iter) {
            int n = rng.range(50, 250);
            std::vector<int> color(n + 1, 0);
            for (int u = 1; u <= n; ++u) color[u] = rng.range(1, n);
            std::vector<int> parent(n + 1, 0);
            for (int u = 2; u <= n; ++u) parent[u] = rng.range(1, u - 1);
            std::vector<std::pair<int, int>> qlist;
            for (int u = 1; u <= n; ++u) qlist.push_back({u, rng.range(0, n)});
            std::map<int, std::vector<std::pair<int, int>>> queries;
            for (int i = 0; i < static_cast<int>(qlist.size()); ++i) {
                queries[qlist[i].first].push_back({qlist[i].second, i});
            }
            Solver s;
            s.build(n, color, parent, queries);
            if (s.run() != solveBrute(n, color, parent, qlist)) {
                fail("mid random mismatch");
            }
        }
    }
    // 深鏈：驗證迭代實作不爆堆
    {
        int n = 5000;
        std::vector<int> color(n + 1, 0);
        for (int u = 1; u <= n; ++u) color[u] = (u % 3) + 1;
        std::vector<int> parent(n + 1, 0);
        for (int u = 2; u <= n; ++u) parent[u] = u - 1;
        std::vector<std::pair<int, int>> qlist;
        for (int u = 1; u <= n; ++u) qlist.push_back({u, 1});
        std::map<int, std::vector<std::pair<int, int>>> queries;
        for (int i = 0; i < static_cast<int>(qlist.size()); ++i) {
            queries[qlist[i].first].push_back({qlist[i].second, i});
        }
        Solver s;
        s.build(n, color, parent, queries);
        if (s.run() != solveBrute(n, color, parent, qlist)) {
            fail("deep chain failed");
        }
    }
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& tok) {
    if (tok.empty()) return std::nullopt;
    std::size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i == tok.size()) return std::nullopt;
    for (std::size_t j = i; j < tok.size(); ++j) {
        if (tok[j] < '0' || tok[j] > '9') return std::nullopt;
    }
    // 用 __int128 累加避免溢出；超過 long long 視為非法
    __int128 val = 0;
    bool neg = (tok[0] == '-');
    for (std::size_t j = i; j < tok.size(); ++j) {
        val = val * 10 + (tok[j] - '0');
        if (val > (static_cast<__int128>(1) << 100)) return std::nullopt;
    }
    if (neg) val = -val;
    if (val > static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    if (val < -static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    return static_cast<long long>(val);
}

void runIo(const std::string& raw) {
    std::vector<std::string> toks;
    {
        std::size_t i = 0;
        while (i < raw.size()) {
            while (i < raw.size() && static_cast<unsigned char>(raw[i]) <= ' ') ++i;
            if (i >= raw.size()) break;
            std::size_t j = i;
            while (j < raw.size() && static_cast<unsigned char>(raw[j]) > ' ') ++j;
            toks.push_back(raw.substr(i, j - i));
            i = j;
        }
    }
    std::size_t pos = 0;
    auto nxt = [&toks, &pos]() -> const std::string* {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    auto bad = []() { std::cout << "error\n"; };

    const std::string* pt = nxt();
    std::optional<long long> nOpt = pt ? parseInt(*pt) : std::nullopt;
    if (!nOpt || *nOpt < 1 || *nOpt > NMAX) {
        bad();
        return;
    }
    int n = static_cast<int>(*nOpt);
    std::vector<int> color(n + 1, 0);
    for (int u = 1; u <= n; ++u) {
        pt = nxt();
        std::optional<long long> v = pt ? parseInt(*pt) : std::nullopt;
        if (!v || *v < 1 || *v > n) {
            bad();
            return;
        }
        color[u] = static_cast<int>(*v);
    }
    std::vector<int> parent(n + 1, 0);
    for (int u = 2; u <= n; ++u) {
        pt = nxt();
        std::optional<long long> p = pt ? parseInt(*pt) : std::nullopt;
        if (!p || *p < 1 || *p >= u) {
            bad();
            return;
        }
        parent[u] = static_cast<int>(*p);
    }
    pt = nxt();
    std::optional<long long> qOpt = pt ? parseInt(*pt) : std::nullopt;
    if (!qOpt || *qOpt < 1 || *qOpt > QMAX) {
        bad();
        return;
    }
    int q = static_cast<int>(*qOpt);

    std::vector<std::pair<int, int>> qlist;
    std::map<int, std::vector<std::pair<int, int>>> queries;
    for (int i = 0; i < q; ++i) {
        pt = nxt();
        std::optional<long long> uOpt = std::nullopt;
        if (pt) uOpt = parseInt(*pt);
        const std::string* ft = nxt();
        std::optional<long long> fOpt = std::nullopt;
        if (ft) fOpt = parseInt(*ft);
        if (!uOpt || !fOpt || *uOpt < 1 || *uOpt > n || *fOpt < 0 || *fOpt > n) {
            bad();
            return;
        }
        int u = static_cast<int>(*uOpt);
        int f = static_cast<int>(*fOpt);
        qlist.push_back({u, f});
        queries[u].push_back({f, i});
    }

    Solver solver;
    solver.build(n, color, parent, queries);
    std::vector<Answer> answers = solver.run();
    for (const Answer& a : answers) {
        std::cout << a.distinct << " " << a.eq << " " << a.shared << "\n";
    }
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (static_cast<unsigned char>(raw[i]) > ' ') {
            blank = false;
            break;
        }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
