// 虛擬樹（Virtual Tree）：把 O(n) 的原樹壓縮成 O(k) 的小樹
//
// 題意：
//     給定一棵有 n 個節點的根樹（根為 1），節點 u 有權重 w[u]。
//     有 m 次詢問，每次詢問給出 k 個「關鍵節點」（key nodes）。
//     對每次詢問，要建立這 k 個關鍵節點的**虛擬樹**並回答：
//         (1) 虛擬樹共有幾個節點
//         (2) 虛擬樹所有邊的長度總和（邊長 = 原樹上兩端點的距離）
//         (3) 對每個關鍵節點 u（依節點編號遞增輸出）：
//             cnt = 位於 u 的虛擬樹子樹中的關鍵節點個數
//             sum = 位於 u 的虛擬樹子樹中的關鍵節點的權重總和
//
//     限制：1 <= n <= 100000，0 <= m <= 100000，0 <= k <= n，
//           |w[u]| <= 10^9。
//
// 思路：
//     ### 1. 為什麼需要虛擬樹
//     原樹有 n 個節點，但每次詢問只關心 k 個關鍵節點（k << n）。
//     若每次都在原樹上做 DFS，複雜度是 O(nm)，完全不可行。
//     虛擬樹的想法是：**只保留 k 個關鍵節點，以及它們彼此的 LCA**，
//     把中間沒有關鍵節點的部分壓縮成一條邊。壓縮後的樹最多 2k-1 個節點。
//
//     ### 2. 虛擬樹的節點集合
//     令 V = {關鍵節點} ∪ {LCA(相鄰的兩個關鍵節點)}。
//     這裡「相鄰」指的是把關鍵節點依 DFS 序排序後相鄰的那一對。
//     定理：這 k-1 個 LCA 恰好涵蓋了所有需要的內部節點 ——
//     因為一個節點若是兩個關鍵節點的 LCA，那麼在 DFS 序上它們之間
//     必然存在一對「相鄰且分居兩側」的關鍵節點。
//     所以 |V| <= 2k-1。
//
//     ### 3. 建構虛擬樹
//     把 V 依 DFS 序排序，用一個堆疊線性建樹：
//     對每個節點 v，把堆疊頂端「不是 v 的祖先」的節點彈掉，
//     剩下的堆疊頂就是 v 在虛擬樹上的父親。
//     因為 V 已依 DFS 序排序，任何時刻堆疊中的節點都是一條祖先鏈。
//     邊長 = depth[child] - depth[parent]。
//
//     ### 4. 樹上統計
//     虛擬樹建好後，cnt / sum 就是普通的樹上 DP：
//     先把關鍵節點的自身權重放進去，再由子往父累加。
//     因為虛擬樹的節點依 DFS 序排序時「父一定排在子前面」，
//     倒序掃描一次即可完成累加。
//
//     ### 5. 為什麼 cnt / sum 跟原樹子樹的答案一樣
//     虛擬樹保留了關鍵節點之間的祖先關係：若 u 是 v 的祖先
//     （u、v 都是關鍵節點），則 u 在虛擬樹上依然是 v 的祖先。
//     所以「u 的虛擬樹子樹中的關鍵節點」與「u 的原樹子樹中的關鍵節點」
//     是同一個集合。
//
// 複雜度：
//     預處理（LCA 倍增表 + DFS 序）  時間 O(n log n)，空間 O(n log n)
//     單次詢問                      O(k log k + k log n)
//     總計                          O((n + sum(k)) log n)
//
// 應用場景：
//     - 樹上批量路徑查詢（Virtual Tree / 虛擬樹優化）：把 O(n) 的一次遍歷
//       壓縮成 O(k)，是「多次少量查詢」的標準加速手段
//     - 靜態 LCA 問題（ZJOI / SNOI 的「接近的祖先」類題目）
//     - 集合運算型子樹查詢：合併子集合、統計子樹內不同元素的個數
//     - 動態規劃換根 / 樹形 DP 中的「壓縮樹」技巧
//     - 莫隊 + 樹上啟發式合併的基礎設施

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr int MMAX = 100000;
constexpr long long WMAX = 1000000000LL;
constexpr long long PARSE_LIMIT = 1000000000000000000LL;

struct Tree {
    int n = 0;
    std::vector<int> parent;
    std::vector<int> depth;
    std::vector<int> tin;
    std::vector<int> tout;
    std::vector<int> order;                      // DFS 先序
    std::vector<std::vector<int>> children;
    std::vector<std::vector<int>> up;           // up[j][v] = 2^j 級祖先
    std::vector<long long> weight;

    void build(int n_, const std::vector<int>& parent_,
               const std::vector<long long>& weight_) {
        n = n_;
        parent = parent_;
        weight = weight_;
        children.assign(n + 1, std::vector<int>());
        for (int u = 2; u <= n; ++u) {
            children[parent[u]].push_back(u);
        }
        depth.assign(n + 1, 0);
        tin.assign(n + 1, 0);
        tout.assign(n + 1, 0);
        order.clear();
        order.reserve(n);
        // 迭代 DFS：避免深鏈造成堆疊溢出
        std::vector<int> st;
        st.push_back(1);
        while (!st.empty()) {
            int u = st.back();
            st.pop_back();
            tin[u] = static_cast<int>(order.size());
            order.push_back(u);
            for (int c : children[u]) {
                depth[c] = depth[u] + 1;
                st.push_back(c);
            }
        }
        // 子樹大小 -> tout（DFS 序中子樹是一段連續區間）
        std::vector<int> size(n + 1, 1);
        for (int i = n - 1; i >= 0; --i) {
            int u = order[i];
            for (int c : children[u]) {
                size[u] += size[c];
            }
            tout[u] = tin[u] + size[u];
        }
        // 倍增表
        int log = 1;
        while ((1 << log) <= n) ++log;
        up.assign(log, std::vector<int>(n + 1, 0));
        up[0] = parent;
        for (int j = 1; j < log; ++j) {
            for (int v = 1; v <= n; ++v) {
                up[j][v] = up[j - 1][up[j - 1][v]];
            }
        }
    }

    bool isAncestor(int u, int v) const {
        return tin[u] <= tin[v] && tin[v] < tout[u];
    }

    int lca(int a, int b) const {
        if (depth[a] < depth[b]) std::swap(a, b);
        int diff = depth[a] - depth[b];
        int j = 0;
        while (diff) {
            if (diff & 1) a = up[j][a];
            diff >>= 1;
            ++j;
        }
        if (a == b) return a;
        for (int i = static_cast<int>(up.size()) - 1; i >= 0; --i) {
            if (up[i][a] != up[i][b]) {
                a = up[i][a];
                b = up[i][b];
            }
        }
        return parent[a];
    }
};

// 虛擬樹建構結果
struct VirtualTree {
    std::vector<int> nodes;    // 依 DFS 序遞增
    std::vector<int> vparent;  // 在 nodes 中的索引，根為 -1
};

VirtualTree buildVirtual(const Tree& tree, const std::vector<int>& keys) {
    VirtualTree vt;
    if (keys.empty()) return vt;
    std::vector<int> skeys = keys;
    std::sort(skeys.begin(), skeys.end(),
              [&tree](int a, int b) { return tree.tin[a] < tree.tin[b]; });
    std::vector<int> nodes = skeys;
    for (std::size_t i = 0; i + 1 < skeys.size(); ++i) {
        nodes.push_back(tree.lca(skeys[i], skeys[i + 1]));
    }
    std::sort(nodes.begin(), nodes.end(),
              [&tree](int a, int b) { return tree.tin[a] < tree.tin[b]; });
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    vt.nodes = nodes;
    vt.vparent.assign(nodes.size(), -1);
    std::vector<int> st;  // 存 nodes 的索引
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        while (!st.empty() && !tree.isAncestor(nodes[st.back()], nodes[i])) {
            st.pop_back();
        }
        if (!st.empty()) vt.vparent[i] = st.back();
        st.push_back(static_cast<int>(i));
    }
    return vt;
}

struct QueryResult {
    int nodeCount = 0;
    long long edgeTotal = 0;
    std::vector<std::pair<long long, long long>> stats;  // (cnt, sum)
};

QueryResult solveQuery(const Tree& tree, const std::vector<int>& keys) {
    QueryResult res;
    VirtualTree vt = buildVirtual(tree, keys);
    std::size_t k = vt.nodes.size();
    if (k == 0) return res;
    std::set<int> keyset(keys.begin(), keys.end());
    std::vector<long long> cnt(k, 0), wsum(k, 0);
    for (std::size_t i = 0; i < k; ++i) {
        if (keyset.count(vt.nodes[i])) {
            cnt[i] = 1;
            wsum[i] = tree.weight[vt.nodes[i]];
        }
    }
    long long edgeTotal = 0;
    // 倒序掃描：子節點的 DFS 序必定大於父節點，故子先於父被處理
    for (std::size_t i = k; i-- > 0;) {
        int p = vt.vparent[i];
        if (p >= 0) {
            cnt[p] += cnt[i];
            wsum[p] += wsum[i];
            edgeTotal += tree.depth[vt.nodes[i]] - tree.depth[vt.nodes[vt.vparent[i]]];
        }
    }
    std::map<int, std::size_t> pos;
    for (std::size_t i = 0; i < k; ++i) pos[vt.nodes[i]] = i;
    std::vector<int> sortedKeys = keys;
    std::sort(sortedKeys.begin(), sortedKeys.end());
    for (int u : sortedKeys) {
        std::size_t idx = pos[u];
        res.stats.emplace_back(cnt[idx], wsum[idx]);
    }
    res.nodeCount = static_cast<int>(k);
    res.edgeTotal = edgeTotal;
    return res;
}

// 暴力參考實作：O(k^2) 求虛擬樹節點集合與每個關鍵節點的統計
QueryResult bruteQuery(const Tree& tree, const std::vector<int>& keys) {
    QueryResult res;
    if (keys.empty()) return res;
    std::set<int> vset(keys.begin(), keys.end());
    for (std::size_t i = 0; i < keys.size(); ++i) {
        for (std::size_t j = 0; j < keys.size(); ++j) {
            vset.insert(tree.lca(keys[i], keys[j]));
        }
    }
    long long edgeTotal = 0;
    for (int v : vset) {
        int bestAncestor = -1;
        for (int u : vset) {
            if (u != v && tree.isAncestor(u, v)) {
                if (bestAncestor < 0 || tree.depth[u] > tree.depth[bestAncestor]) {
                    bestAncestor = u;
                }
            }
        }
        if (bestAncestor >= 0) {
            edgeTotal += tree.depth[v] - tree.depth[bestAncestor];
        }
    }
    std::vector<int> sortedKeys = keys;
    std::sort(sortedKeys.begin(), sortedKeys.end());
    for (int u : sortedKeys) {
        long long c = 0, s = 0;
        for (int v : keys) {
            if (tree.isAncestor(u, v)) {
                ++c;
                s += tree.weight[v];
            }
        }
        res.stats.emplace_back(c, s);
    }
    res.nodeCount = static_cast<int>(vset.size());
    res.edgeTotal = edgeTotal;
    return res;
}

bool sameResult(const QueryResult& a, const QueryResult& b) {
    return a.nodeCount == b.nodeCount && a.edgeTotal == b.edgeTotal &&
           a.stats == b.stats;
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
        return lo + static_cast<int>(next() % static_cast<unsigned long long>(hi - lo + 1));
    }
};

void runTests() {
    // ---- 固定用例 ----
    // 樹：1-2, 1-3, 2-4, 2-5, 3-6（根 = 1）
    {
        std::vector<int> parent = {0, 0, 1, 1, 2, 2, 3};
        std::vector<long long> weight = {0, 10, 20, 30, 40, 50, 60};
        Tree tree;
        tree.build(6, parent, weight);

        // 單一關鍵節點：虛擬樹只有它自己
        {
            std::vector<int> keys = {4};
            QueryResult r = solveQuery(tree, keys);
            if (!(r.nodeCount == 1 && r.edgeTotal == 0 &&
                  r.stats.size() == 1 && r.stats[0].first == 1 &&
                  r.stats[0].second == 40)) {
                std::cerr << "case 1 failed\n";
                std::exit(1);
            }
            if (!sameResult(r, bruteQuery(tree, keys))) {
                std::cerr << "case 1 brute mismatch\n";
                std::exit(1);
            }
        }
        // 兩個兄弟節點 4、5：LCA = 2，虛擬樹是 2-4、2-5，邊長 1+1 = 2
        {
            std::vector<int> keys = {4, 5};
            QueryResult r = solveQuery(tree, keys);
            if (!(r.nodeCount == 3 && r.edgeTotal == 2 && r.stats.size() == 2 &&
                  r.stats[0].first == 1 && r.stats[0].second == 40 &&
                  r.stats[1].first == 1 && r.stats[1].second == 50)) {
                std::cerr << "case 2 failed\n";
                std::exit(1);
            }
            if (!sameResult(r, bruteQuery(tree, keys))) {
                std::cerr << "case 2 brute mismatch\n";
                std::exit(1);
            }
        }
        // 祖先與後代：3 與 6，LCA = 3；3 的子樹含 {3, 6}
        {
            std::vector<int> keys = {3, 6};
            QueryResult r = solveQuery(tree, keys);
            if (!(r.nodeCount == 2 && r.edgeTotal == 1 && r.stats.size() == 2 &&
                  r.stats[0].first == 2 && r.stats[0].second == 90 &&
                  r.stats[1].first == 1 && r.stats[1].second == 60)) {
                std::cerr << "case 3 failed\n";
                std::exit(1);
            }
            if (!sameResult(r, bruteQuery(tree, keys))) {
                std::cerr << "case 3 brute mismatch\n";
                std::exit(1);
            }
        }
        // 跨分支：4 與 6，LCA = 1，邊長 2+2 = 4
        {
            std::vector<int> keys = {4, 6};
            QueryResult r = solveQuery(tree, keys);
            if (!(r.nodeCount == 3 && r.edgeTotal == 4 && r.stats.size() == 2 &&
                  r.stats[0].first == 1 && r.stats[0].second == 40 &&
                  r.stats[1].first == 1 && r.stats[1].second == 60)) {
                std::cerr << "case 4 failed\n";
                std::exit(1);
            }
            if (!sameResult(r, bruteQuery(tree, keys))) {
                std::cerr << "case 4 brute mismatch\n";
                std::exit(1);
            }
        }
        // 四個關鍵節點：V = {1, 2, 4, 5, 6}
        {
            std::vector<int> keys = {1, 4, 5, 6};
            QueryResult r = solveQuery(tree, keys);
            if (!(r.nodeCount == 5 && r.edgeTotal == 5 && r.stats.size() == 4 &&
                  r.stats[0].first == 4 && r.stats[0].second == 160)) {
                std::cerr << "case 5 failed\n";
                std::exit(1);
            }
            if (!sameResult(r, bruteQuery(tree, keys))) {
                std::cerr << "case 5 brute mismatch\n";
                std::exit(1);
            }
        }
        // 空詢問
        {
            std::vector<int> keys;
            if (!sameResult(solveQuery(tree, keys), bruteQuery(tree, keys))) {
                std::cerr << "case 6 failed\n";
                std::exit(1);
            }
        }
    }
    // 負權重：鏈 1-2-3，權重 -5 / 7 / -1
    {
        std::vector<int> parent = {0, 0, 1, 2};
        std::vector<long long> weight = {0, -5, 7, -1};
        Tree tree;
        tree.build(3, parent, weight);
        std::vector<int> keys = {1, 3};
        QueryResult r = solveQuery(tree, keys);
        if (!(r.nodeCount == 2 && r.edgeTotal == 2 && r.stats.size() == 2 &&
              r.stats[0].first == 2 && r.stats[0].second == -6 &&
              r.stats[1].first == 1 && r.stats[1].second == -1)) {
            std::cerr << "case 7 failed\n";
            std::exit(1);
        }
        if (!sameResult(r, bruteQuery(tree, keys))) {
            std::cerr << "case 7 brute mismatch\n";
                std::exit(1);
        }
    }
    // 單節點樹
    {
        std::vector<int> parent = {0, 0};
        std::vector<long long> weight = {0, 42};
        Tree tree;
        tree.build(1, parent, weight);
        std::vector<int> keys = {1};
        QueryResult r = solveQuery(tree, keys);
        if (!(r.nodeCount == 1 && r.edgeTotal == 0 && r.stats.size() == 1 &&
              r.stats[0].first == 1 && r.stats[0].second == 42)) {
            std::cerr << "case 8 failed\n";
            std::exit(1);
        }
    }

    // ---- 隨機對拍 ----
    {
        Lcg rng(20261030ULL);
        for (int iter = 0; iter < 200; ++iter) {
            int n = rng.range(1, 14);
            std::vector<int> parent(n + 1, 0);
            for (int u = 2; u <= n; ++u) parent[u] = rng.range(1, u - 1);
            std::vector<long long> weight(n + 1, 0);
            for (int u = 1; u <= n; ++u) weight[u] = rng.range(-20, 20);
            Tree tree;
            tree.build(n, parent, weight);
            std::vector<int> keys;
            for (int u = 1; u <= n; ++u) keys.push_back(u);
            // 隨機取子集：用 Fisher-Yates 打亂後取前 k 個
            for (int i = n - 1; i > 0; --i) {
                int j = rng.range(0, i);
                std::swap(keys[i], keys[j]);
            }
            int k = rng.range(0, n);
            keys.resize(k);
            QueryResult a = solveQuery(tree, keys);
            QueryResult b = bruteQuery(tree, keys);
            if (!sameResult(a, b)) {
                std::cerr << "random mismatch at iter " << iter << "\n";
                std::exit(1);
            }
        }
    }
    // 較大規模
    {
        Lcg rng(99887766ULL);
        for (int iter = 0; iter < 20; ++iter) {
            int n = rng.range(50, 400);
            std::vector<int> parent(n + 1, 0);
            for (int u = 2; u <= n; ++u) parent[u] = rng.range(1, u - 1);
            std::vector<long long> weight(n + 1, 0);
            for (int u = 1; u <= n; ++u) weight[u] = rng.range(-1000, 1000);
            Tree tree;
            tree.build(n, parent, weight);
            std::vector<int> keys;
            for (int u = 1; u <= n; ++u) keys.push_back(u);
            for (int i = n - 1; i > 0; --i) {
                int j = rng.range(0, i);
                std::swap(keys[i], keys[j]);
            }
            keys.resize(rng.range(1, n));
            QueryResult a = solveQuery(tree, keys);
            QueryResult b = bruteQuery(tree, keys);
            if (!sameResult(a, b)) {
                std::cerr << "big random mismatch at iter " << iter << "\n";
                std::exit(1);
            }
        }
    }
    // 深鏈：驗證迭代 DFS 不爆堆
    {
        int n = 5000;
        std::vector<int> parent(n + 1, 0);
        for (int u = 2; u <= n; ++u) parent[u] = u - 1;
        std::vector<long long> weight(n + 1, 0);
        for (int u = 1; u <= n; ++u) weight[u] = 1;
        Tree tree;
        tree.build(n, parent, weight);
        std::vector<int> keys = {1, n / 2, n};
        QueryResult r = solveQuery(tree, keys);
        if (!(r.nodeCount == 3 && r.edgeTotal == n - 1)) {
            std::cerr << "deep chain failed\n";
            std::exit(1);
        }
        if (!(r.stats.size() == 3 && r.stats[0].first == 3 && r.stats[0].second == 3 &&
              r.stats[2].first == 1 && r.stats[2].second == 1)) {
            std::cerr << "deep chain stats failed\n";
            std::exit(1);
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
    // 以空白切分 token
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
    std::vector<long long> weight(n + 1, 0);
    for (int u = 1; u <= n; ++u) {
        pt = nxt();
        std::optional<long long> v = pt ? parseInt(*pt) : std::nullopt;
        if (!v || *v < -WMAX || *v > WMAX) {
            bad();
            return;
        }
        weight[u] = *v;
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
    std::optional<long long> mOpt = pt ? parseInt(*pt) : std::nullopt;
    if (!mOpt || *mOpt < 0 || *mOpt > MMAX) {
        bad();
        return;
    }
    int m = static_cast<int>(*mOpt);

    Tree tree;
    tree.build(n, parent, weight);

    for (int q = 0; q < m; ++q) {
        pt = nxt();
        std::optional<long long> kOpt = pt ? parseInt(*pt) : std::nullopt;
        if (!kOpt || *kOpt < 0 || *kOpt > n) {
            bad();
            return;
        }
        int k = static_cast<int>(*kOpt);
        std::vector<int> keys;
        std::set<int> seen;
        bool ok = true;
        for (int i = 0; i < k; ++i) {
            pt = nxt();
            std::optional<long long> x = pt ? parseInt(*pt) : std::nullopt;
            if (!x || *x < 1 || *x > n || seen.count(static_cast<int>(*x))) {
                ok = false;
                break;
            }
            seen.insert(static_cast<int>(*x));
            keys.push_back(static_cast<int>(*x));
        }
        if (!ok) {
            bad();
            return;
        }
        QueryResult r = solveQuery(tree, keys);
        std::string line = std::to_string(r.nodeCount) + " " + std::to_string(r.edgeTotal);
        for (const auto& st : r.stats) {
            line += " " + std::to_string(st.first);
            line += " " + std::to_string(st.second);
        }
        std::cout << line << "\n";
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
