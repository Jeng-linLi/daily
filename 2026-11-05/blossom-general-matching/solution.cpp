// 帶花樹 Blossom：一般圖最大匹配（Edmonds 花收縮）
//
// 題意：
//     給定 n 個點、m 條邊的無向簡單圖（可能有奇環，不限二分圖），
//     求最大匹配：最多的兩兩不共點的邊集合。
//     輸出匹配數，以及每個點的配對點（未匹配點輸出 0）。
//
//     輸入：第一行 n m；接下來 m 行每行兩個整數 u v（1-indexed）。
//     輸出：第一行一個整數 k；第二行 n 個整數，第 i 個是點 i 的配對點。
//
// 思路：
//     二分圖的匈牙利算法靠交錯路徑增廣，但一般圖中交錯路徑可能經過
//     奇環，交替染色失效。Edmonds 1965 提出把「花」——從樹根出發的
//     交錯路徑 + 一條閉合奇環的交錯邊——整朵收縮成一個新節點：
//     花上任意一點都能以正確的奇偶性延伸增廣路，因此收縮後繼續 BFS
//     不丟解。找到增廣路後再沿 parent 鏈展開花、翻轉匹配邊。
//
//     實作要點：
//     - base[i]：i 所在花的代表點（並查集，帶路徑壓縮）
//     - p[i]：交錯森林中 i 的父親（i 是外點時有意義）
//     - 每輪從未匹配點 BFS：外點把鄰居分類——未拜訪 → 找增廣路或
//       把其配偶入隊（內點）；拜訪過的外點且 base 不同 → 收縮新花
//
// 複雜度：
//     每個未匹配點最多增廣一次，單次找增廣路 O(n^2)～O(n^3)，
//     總時間 O(n^3)，空間 O(n + m)。
//
// 應用場景（Applications）：
//     - 舞伴配對、宿舍分配、導師-學生互選（圖不是二分圖）
//     - 網路設計中的邊不相交配對規劃
//
// 內建測試：三角形、雙三角形、Petersen 圖等固定用例 + 300 組
//     隨機圖（n <= 10）與「枚舉邊子集」暴力對拍。

#include <algorithm>
#include <cstdint>
#include <deque>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

const int NMAX = 300;
const long long VAL_LIM = 1000000000000000000LL;

struct Blossom {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> match;  // match[u] = 與 u 配對的點（-1 未匹配）
    std::vector<int> p;      // 交錯森林父親（外點）
    std::vector<int> base;   // 所在花的代表（並查集根）
    std::vector<bool> used;

    Blossom(int n_, const std::vector<std::vector<int>>& adj_) : n(n_), adj(adj_) {
        match.assign(n, -1);
        p.assign(n, 0);
        base.assign(n, 0);
        used.assign(n, false);
    }

    // 並查集（帶路徑壓縮）
    int findBase(int x) {
        while (base[x] != x) {
            base[x] = base[base[x]];
            x = base[x];
        }
        return x;
    }

    // 求兩個外點在交錯森林中的最近公共花代表
    int lca(int a, int b) {
        std::vector<bool> mark(n, false);
        while (true) {
            a = findBase(a);
            mark[a] = true;
            if (match[a] == -1) break;
            a = p[match[a]];
        }
        while (true) {
            b = findBase(b);
            if (mark[b]) return b;
            b = p[match[b]];
        }
    }

    // 把 v 到花代表 b 的路徑全標成花，child 是閉合邊的另一端
    void markPath(int v, int b, int child, std::vector<bool>& blossom) {
        while (findBase(v) != b) {
            blossom[findBase(v)] = true;
            blossom[findBase(match[v])] = true;
            p[v] = child;
            child = match[v];
            v = p[match[v]];
        }
    }

    // 從 root（未匹配點）BFS 找一條增廣路；找到並翻轉後回 true
    bool findPath(int root) {
        used.assign(n, false);
        p.assign(n, -1);
        for (int i = 0; i < n; ++i) base[i] = i;
        used[root] = true;
        std::deque<int> q;
        q.push_back(root);
        while (!q.empty()) {
            int v = q.front();
            q.pop_front();
            for (int to : adj[v]) {
                // 同一朵花內、或 v-to 本身是匹配邊：跳過
                if (findBase(v) == findBase(to) || match[v] == to) continue;
                // 情況 1：to 是已被拜訪過的外點（或樹根）→ 找到花
                if (to == root || (match[to] != -1 && p[match[to]] != -1)) {
                    int curbase = lca(v, to);
                    std::vector<bool> blossom(n, false);
                    markPath(v, curbase, to, blossom);
                    markPath(to, curbase, v, blossom);
                    for (int i = 0; i < n; ++i) {
                        if (blossom[findBase(i)]) {
                            base[i] = curbase;
                            if (!used[i]) {
                                used[i] = true;
                                q.push_back(i);
                            }
                        }
                    }
                }
                // 情況 2：to 還沒在交錯森林裡
                else if (p[to] == -1) {
                    p[to] = v;
                    if (match[to] == -1) {
                        // 找到增廣路：沿 p 鏈倒序翻轉匹配邊
                        int u = to;
                        while (u != -1) {
                            int pv = p[u];
                            int ppv = match[pv];
                            match[u] = pv;
                            match[pv] = u;
                            u = ppv;
                        }
                        return true;
                    } else {
                        // to 已匹配 → 其配偶作為內點延伸
                        used[match[to]] = true;
                        q.push_back(match[to]);
                    }
                }
            }
        }
        return false;
    }

    int solve() {
        int res = 0;
        for (int v = 0; v < n; ++v) {
            if (match[v] == -1 && findPath(v)) ++res;
        }
        return res;
    }
};

// 一般圖最大匹配：回傳 (匹配數, match 陣列；未匹配為 0)
std::pair<int, std::vector<int>> generalMatching(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> adj(n);
    for (auto [a, b] : edges) {
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    Blossom bl(n, adj);
    int res = bl.solve();
    std::vector<int> out(n);
    for (int i = 0; i < n; ++i) out[i] = (bl.match[i] == -1) ? 0 : bl.match[i] + 1;
    return {res, out};
}

// --------------------------------------------------------------- 暴力對拍

// 枚舉邊子集（m <= 15），取合法匹配的最大值
long long bruteMaxMatching(int n, const std::vector<std::pair<int, int>>& edges) {
    int m = (int)edges.size();
    long long best = 0;
    for (int mask = 0; mask < (1 << m); ++mask) {
        std::vector<int> deg(n, 0);
        int cnt = 0;
        bool ok = true;
        for (int i = 0; i < m; ++i) {
            if (!((mask >> i) & 1)) continue;
            int a = edges[i].first, b = edges[i].second;
            if (++deg[a] > 1 || ++deg[b] > 1) {
                ok = false;
                break;
            }
            ++cnt;
        }
        if (ok && cnt > best) best = cnt;
    }
    return best;
}

// --------------------------------------------------------------- 內建測試

// LCG 確定性隨機數（與 Python 版完全一致）
static uint32_t g_rng = 987654321u;
int rnd(int mod) {
    g_rng = (uint32_t)((g_rng * 1103515245ULL + 12345ULL) % 2147483648ULL);
    return (int)(g_rng % (uint32_t)mod);
}

void fail(const std::string& msg) {
    std::cerr << msg << std::endl;
    std::exit(1);
}

void runTests() {
    // 固定用例 1：三角形（奇環）→ 最大匹配 1
    {
        auto [res, match] = generalMatching(3, {{0, 1}, {1, 2}, {2, 0}});
        if (res != 1) fail("case1 count");
        std::vector<int> nz;
        for (int x : match)
            if (x) nz.push_back(x);
        if (nz.size() != 2 || match[nz[0] - 1] != nz[1]) fail("case1 match");
    }

    // 固定用例 2：兩個三角形由一條邊相連 → 匹配數 3
    {
        auto [res, match] = generalMatching(6, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {2, 3}});
        if (res != 3) fail("case2");
    }

    // 固定用例 3：Petersen 圖（10 點 15 邊）→ 有完美匹配，答案 5
    {
        std::vector<std::pair<int, int>> pet = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
                                                {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9},
                                                {5, 7}, {7, 9}, {9, 6}, {6, 8}, {8, 5}};
        auto [res, match] = generalMatching(10, pet);
        if (res != 5) fail("case3");
    }

    // 固定用例 4：無邊
    {
        auto [res, match] = generalMatching(4, {});
        if (res != 0 || match != std::vector<int>({0, 0, 0, 0})) fail("case4");
    }

    // 固定用例 5：奇環 + 尾巴（經典 blossom 觸發形態）
    {
        auto [res, match] = generalMatching(5, {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}});
        if (res != 2) fail("case5");
    }

    // 隨機對拍（n <= 10, m <= 15）
    for (int it = 0; it < 300; ++it) {
        int n = rnd(10) + 1;
        // 枚舉所有點對後 Fisher-Yates 洗牌取前 m條：
        // 直接 while 循環抽點在 LCG 低比特周期短時會湊不滿 m 條而死循環
        std::vector<std::pair<int, int>> allEdges;
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b) allEdges.push_back({a, b});
        for (int i = (int)allEdges.size() - 1; i > 0; --i) {
            int j = rnd(i + 1);
            std::swap(allEdges[i], allEdges[j]);
        }
        int m = rnd((int)std::min<size_t>(allEdges.size(), 15) + 1);
        std::vector<std::pair<int, int>> edges(allEdges.begin(), allEdges.begin() + m);
        auto [got, match] = generalMatching(n, edges);
        long long want = bruteMaxMatching(n, edges);
        if ((long long)got != want) fail("random mismatch at it=" + std::to_string(it));
        // 自洽：match 陣列必須互相配對且無共用點
        for (int a = 0; a < n; ++a) {
            if (match[a] != 0 && match[match[a] - 1] != a + 1) fail("random match broken");
        }
    }
}

// --------------------------------------------------------------- IO 模式

std::optional<long long> parseLL(const std::string& tok) {
    if (tok.empty()) return std::nullopt;
    size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i == tok.size()) return std::nullopt;
    for (size_t k = i; k < tok.size(); ++k) {
        if (tok[k] < '0' || tok[k] > '9') return std::nullopt;
    }
    // 手動轉換避免 strtoll 對超長數字的行為差異
    long long val = 0;
    for (size_t k = i; k < tok.size(); ++k) {
        val = val * 10 + (tok[k] - '0');
        if (val > VAL_LIM) return std::nullopt;
    }
    if (tok[0] == '-') val = -val;
    return val;
}

std::vector<std::string> splitTokens(const std::string& raw) {
    std::vector<std::string> toks;
    std::string cur;
    for (char c : raw) {
        if (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
            if (!cur.empty()) toks.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) toks.push_back(cur);
    return toks;
}

void runIO(const std::string& raw) {
    std::vector<std::string> toks = splitTokens(raw);
    size_t pos = 0;
    auto nxt = [&]() -> std::string { return pos < toks.size() ? toks[pos++] : std::string(); };
    auto bad = []() { std::cout << "error\n"; };

    auto nv = parseLL(nxt());
    if (!nv || *nv < 0 || *nv > NMAX) {
        bad();
        return;
    }
    int n = (int)*nv;
    auto mv = parseLL(nxt());
    if (!mv || *mv < 0 || *mv > (long long)n * NMAX) {
        bad();
        return;
    }
    int m = (int)*mv;
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < m; ++i) {
        auto a = parseLL(nxt());
        auto b = parseLL(nxt());
        if (!a || !b || *a < 1 || *a > n || *b < 1 || *b > n || *a == *b) {
            bad();
            return;
        }
        edges.push_back({(int)*a - 1, (int)*b - 1});
    }
    auto [res, match] = generalMatching(n, edges);
    std::cout << res << "\n";
    for (int i = 0; i < n; ++i) {
        if (i) std::cout << ' ';
        std::cout << match[i];
    }
    std::cout << "\n";
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw, line;
    while (std::getline(std::cin, line)) {
        raw += line;
        raw += '\n';
    }
    bool blank = true;
    for (char c : raw) {
        if (c != ' ' && c != '\n' && c != '\r' && c != '\t') {
            blank = false;
            break;
        }
    }
    if (blank) {
        runTests();
        std::cout << "all tests passed" << std::endl;
    } else {
        runIO(raw);
    }
    return 0;
}
