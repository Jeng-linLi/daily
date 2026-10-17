// A* 與 IDA*（啟發式搜索）
//
// 題意：
//     三小題，一次輸入全部做完，每一題都同時給出「啟發式搜索」與「無資訊的最優對照組」並對拍。
//
//     1. 網格最短路：R×C 格子，'#' 是障礙，四方向移動、每步代價 1。
//        A*（曼哈頓距離啟發）vs BFS，比較「最短路長度」與「展開節點數」。
//     2. 帶權圖最短路：n 點 m 條有向邊（邊權非負），外加每個點的啟發值 h[v]
//        （h 為到終點的估計距離）。A* vs Dijkstra，同樣比較距離與展開節點數。
//     3. 八碼難題（8-puzzle）：3×3 棋盤，給定初始局面，用 IDA* 求最少步數；
//        與「從目標局面 BFS 得到的精確距離」對照（BFS 只做到深度 12，查不到就標記未知）。
//
// 思路：
//     ### A*：把「已經花的代價 g」和「預估還要花的代價 h」加起來排序
//     優先佇列依照 f = g + h 由小到大取出節點。
//       - h ≡ 0 時，A* 就退化成 Dijkstra（本題的對照組）。
//       - h 若**可採納**（h(v) ≤ 真實距離 d(v,t)）且**一致**（h(u) ≤ w(u,v) + h(v)），
//         A* 保證第一次取出終點就是最短路。網格四方向 + 曼哈頓距離正好滿足這兩個條件。
//       - A* 的價值在於**展開的節點更少**：h 越準，搜索越像「直接朝終點走」。
//     為了讓 Python 與 C++ 的展開數完全一致，佇列的排序鍵固定為 (f, h, 節點編號)，
//     鄰居的枚舉順序也固定為 上 / 下 / 左 / 右。
//
//     ### IDA*：迭代加深 + 記憶「這輪超標多少」當下一輪門檻
//     A* 要存 open/closed 表，八碼難題的狀態太多會爆炸。IDA* 用 DFS 換空間：
//       - 設門檻 bound，DFS 只展開 f = g + h ≤ bound 的節點；
//       - 這輪找不到就記錄「最小的超標值」當下一輪的 bound，從頭再搜一次；
//       - 因為每輪門檻只增不減，且 h 可採納，**第一次找到的解一定是最優解**。
//     空間複雜度只有 O(最長路徑)，代價是重複搜索（但每輪的分支都被 h 砍掉大半）。
//     本題的 h 用「曼哈頓距離 + 線性衝突」：
//       曼哈頓距離給出每塊磁磚至少要走幾步；
//       同一行（列）裡有兩塊磁磚的目標位置也在這一行（列）但相對順序相反時，
//       它們必定要「繞出去」一次，至少多 2 步。
//
//     ### 八碼難題的可解性
//     3×3 棋盤寬度為奇數，所以「逆序數（不含空格）為偶數」↔ 可解。
//     不可解時直接回報 -1，不必展開任何節點。
//
// 應用場景：
//     - A*：遊戲 NPC 尋路（網格 / 導航網）、地圖導航（真實路網 + 直線距離或 landmark 啟發）、
//       機械手臂與機器人路徑規劃、拼圖 / 魔方的狀態空間搜索
//     - 權衡：h 越強（越接近真實距離）展開越少，但計算 h 本身也要時間；h = 0 是 Dijkstra；
//       h 過強且不可採納則會得到次優解
//     - IDA*：狀態空間大、記憶體吃緊時的首選（15-puzzle、分支限界），
//       也是「迭代加深 + 啟發剪枝」這個通用框架的典型代表
//
// 複雜度：
//     網格 A* / BFS        最壞 O(R·C) 時間 / O(R·C) 空間
//     帶權圖 A* / Dijkstra  O((n + m) log n) 時間 / O(n) 空間
//     IDA*                 時間 O(b^d)（b 分支數、d 解答深度），空間 O(d)

#include <algorithm>
#include <array>
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const int MAXR = 20;                 // IO 模式下網格最大邊長
const int MAXN = 40;                 // 帶權圖最大點數
const int MAXM = 200;                // 最大邊數
const int WMAX = 1000;               // 邊權上限
const long long HMAX = 1000000;      // 啟發值上限
const int IDA_LIMIT = 500000;       // IDA* 節點預算（超標就放棄，回報 -2）
const int ORACLE_DEPTH = 12;         // 對照用 BFS 只做到這個深度
const int INF = 1000000000;

using Board = std::array<int, 9>;
const Board GOAL = {{1, 2, 3, 4, 5, 6, 7, 8, 0}};

// ------------------------------------------------------- 八碼難題：前置表

std::array<std::vector<int>, 9> buildMoveTable() {
    std::array<std::vector<int>, 9> tbl;
    const int dr[4] = {-1, 1, 0, 0};      // 上 / 下 / 左 / 右，順序固定
    const int dc[4] = {0, 0, -1, 1};
    for (int p = 0; p < 9; ++p) {
        int r = p / 3, c = p % 3;
        for (int k = 0; k < 4; ++k) {
            int nr = r + dr[k], nc = c + dc[k];
            if (0 <= nr && nr < 3 && 0 <= nc && nc < 3) tbl[p].push_back(nr * 3 + nc);
        }
    }
    return tbl;
}
const std::array<std::vector<int>, 9> MOVE_TABLE = buildMoveTable();

int manhattan(const Board& st) {
    int total = 0;
    for (int i = 0; i < 9; ++i) {
        int tile = st[i];
        if (tile == 0) continue;
        int gr = (tile - 1) / 3, gc = (tile - 1) % 3;
        int r = i / 3, c = i % 3;
        total += std::abs(gr - r) + std::abs(gc - c);
    }
    return total;
}

int linearConflict(const Board& st) {
    int extra = 0;
    for (int line = 0; line < 3; ++line) {
        std::vector<int> rows;
        for (int c = 0; c < 3; ++c) {
            int t = st[line * 3 + c];
            if (t != 0 && (t - 1) / 3 == line) rows.push_back(t);
        }
        for (size_t i = 0; i < rows.size(); ++i) {
            for (size_t j = i + 1; j < rows.size(); ++j) {
                if ((rows[i] - 1) % 3 > (rows[j] - 1) % 3) extra += 2;
            }
        }
        std::vector<int> cols;
        for (int r = 0; r < 3; ++r) {
            int t = st[r * 3 + line];
            if (t != 0 && (t - 1) % 3 == line) cols.push_back(t);
        }
        for (size_t i = 0; i < cols.size(); ++i) {
            for (size_t j = i + 1; j < cols.size(); ++j) {
                if ((cols[i] - 1) / 3 > (cols[j] - 1) / 3) extra += 2;
            }
        }
    }
    return extra;
}

int heur(const Board& st) { return manhattan(st) + linearConflict(st); }

bool solvable(const Board& st) {
    std::vector<int> tiles;
    for (int t : st) {
        if (t != 0) tiles.push_back(t);
    }
    int inv = 0;
    for (size_t i = 0; i < tiles.size(); ++i) {
        for (size_t j = i + 1; j < tiles.size(); ++j) {
            if (tiles[i] > tiles[j]) ++inv;
        }
    }
    return inv % 2 == 0;
}

// ------------------------------------------------------- 八碼難題：IDA*

struct IdaState {
    long long counter = 0;
    bool aborted = false;
};

struct DfsRet {
    int res;        // ≥ 0 找到解；-1 這輪沒找到
    int nxt;        // 下一個門檻（INF 表示不用再試）
};

DfsRet idaDfs(const Board& st, int g, int bound, int prevBlank, int limit, IdaState& io) {
    ++io.counter;
    if (io.counter > limit) {
        io.aborted = true;
        return {-1, INF};
    }
    int f = g + heur(st);
    if (f > bound) return {-1, f};
    if (st == GOAL) return {g, INF};
    int blank = -1;
    for (int i = 0; i < 9; ++i) {
        if (st[i] == 0) { blank = i; break; }
    }
    int nxtBound = INF;
    for (int nb : MOVE_TABLE[blank]) {
        if (nb == prevBlank) continue;
        Board child = st;
        std::swap(child[blank], child[nb]);
        DfsRet r = idaDfs(child, g + 1, bound, blank, limit, io);
        if (io.aborted) return {-1, INF};
        if (r.res >= 0) return {r.res, INF};
        if (r.nxt < nxtBound) nxtBound = r.nxt;
    }
    return {-1, nxtBound};
}

std::pair<int, long long> idaStar(const Board& start, int limit = IDA_LIMIT) {
    if (!solvable(start)) return {-1, 0};
    if (start == GOAL) return {0, 1};
    IdaState io;
    int bound = heur(start);
    while (bound <= 40) {
        DfsRet r = idaDfs(start, 0, bound, -1, limit, io);
        if (io.aborted) return {-2, io.counter};
        if (r.res >= 0) return {r.res, io.counter};
        if (r.nxt == INF) return {-1, io.counter};
        bound = r.nxt;
    }
    return {-2, io.counter};
}

// ------------------------------------------------------- 八碼難題：BFS 對照組

long long encode(const Board& st) {
    long long key = 0, mul = 1;
    for (int i = 0; i < 9; ++i) {
        key += static_cast<long long>(st[i]) * mul;
        mul *= 9;
    }
    return key;
}

int oracleDist(const Board& st) {
    static std::unordered_map<long long, int> cache;
    static bool built = false;
    if (!built) {
        built = true;
        cache[encode(GOAL)] = 0;
        std::vector<Board> frontier(1, GOAL);
        for (int d = 1; d <= ORACLE_DEPTH; ++d) {
            std::vector<Board> nxt;
            for (const Board& cur : frontier) {
                int blank = -1;
                for (int i = 0; i < 9; ++i) {
                    if (cur[i] == 0) { blank = i; break; }
                }
                for (int nb : MOVE_TABLE[blank]) {
                    Board child = cur;
                    std::swap(child[blank], child[nb]);
                    long long k = encode(child);
                    if (!cache.count(k)) {
                        cache[k] = d;
                        nxt.push_back(child);
                    }
                }
            }
            frontier.swap(nxt);
            if (frontier.empty()) break;
        }
    }
    auto it = cache.find(encode(st));
    return it == cache.end() ? -1 : it->second;
}

// ------------------------------------------------------- 網格：A* 與 BFS

std::pair<int, int> astarGrid(const std::vector<std::string>& grid, int sr, int sc, int tr, int tc) {
    int R = static_cast<int>(grid.size());
    int C = R ? static_cast<int>(grid[0].size()) : 0;
    if (R == 0 || C == 0) return {-1, 0};
    auto h = [&](int r, int c) { return std::abs(r - tr) + std::abs(c - tc); };
    const int INF_D = 1 << 29;
    std::vector<std::vector<int>> dist(R, std::vector<int>(C, INF_D));
    std::vector<std::vector<char>> closed(R, std::vector<char>(C, 0));
    dist[sr][sc] = 0;
    using Item = std::tuple<int, int, int>;          // (f, h, 編號)
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    pq.emplace(h(sr, sc), h(sr, sc), sr * C + sc);
    int expanded = 0;
    const int dr[4] = {-1, 1, 0, 0};
    const int dc[4] = {0, 0, -1, 1};
    while (!pq.empty()) {
        Item cur = pq.top();
        pq.pop();
        int code = std::get<2>(cur);
        int r = code / C, c = code % C;
        if (closed[r][c]) continue;
        closed[r][c] = 1;
        ++expanded;
        if (r == tr && c == tc) break;
        for (int k = 0; k < 4; ++k) {
            int nr = r + dr[k], nc = c + dc[k];
            if (0 <= nr && nr < R && 0 <= nc && nc < C && grid[nr][nc] != '#') {
                if (!closed[nr][nc] && dist[r][c] + 1 < dist[nr][nc]) {
                    dist[nr][nc] = dist[r][c] + 1;
                    int hv = h(nr, nc);
                    pq.emplace(dist[nr][nc] + hv, hv, nr * C + nc);
                }
            }
        }
    }
    int d = dist[tr][tc];
    return {d >= INF_D ? -1 : d, expanded};
}

std::pair<int, int> bfsGrid(const std::vector<std::string>& grid, int sr, int sc, int tr, int tc) {
    int R = static_cast<int>(grid.size());
    int C = R ? static_cast<int>(grid[0].size()) : 0;
    if (R == 0 || C == 0) return {-1, 0};
    const int INF_D = 1 << 29;
    std::vector<std::vector<int>> dist(R, std::vector<int>(C, INF_D));
    dist[sr][sc] = 0;
    std::deque<std::pair<int, int>> q;
    q.emplace_back(sr, sc);
    int expanded = 0;
    const int dr[4] = {-1, 1, 0, 0};
    const int dc[4] = {0, 0, -1, 1};
    while (!q.empty()) {
        int r = q.front().first, c = q.front().second;
        q.pop_front();
        if (dist[r][c] == INF_D) continue;
        ++expanded;
        if (r == tr && c == tc) break;
        for (int k = 0; k < 4; ++k) {
            int nr = r + dr[k], nc = c + dc[k];
            if (0 <= nr && nr < R && 0 <= nc && nc < C && grid[nr][nc] != '#' && dist[nr][nc] == INF_D) {
                dist[nr][nc] = dist[r][c] + 1;
                q.emplace_back(nr, nc);
            }
        }
    }
    int d = dist[tr][tc];
    return {d >= INF_D ? -1 : d, expanded};
}

// ------------------------------------------------------- 帶權圖：A* 與 Dijkstra

using Adj = std::vector<std::vector<std::pair<int, int>>>;

std::pair<int, int> dijkstra(int n, const Adj& adj, int s, int t) {
    const int INF_D = 1 << 29;
    std::vector<int> dist(n, INF_D);
    std::vector<char> closed(n, 0);
    dist[s] = 0;
    using Item = std::pair<int, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    pq.emplace(0, s);
    int expanded = 0;
    while (!pq.empty()) {
        Item cur = pq.top();
        pq.pop();
        int d = cur.first, v = cur.second;
        if (closed[v]) continue;
        closed[v] = 1;
        ++expanded;
        if (v == t) break;
        for (auto& e : adj[v]) {
            int u = e.first, w = e.second;
            if (!closed[u] && d + w < dist[u]) {
                dist[u] = d + w;
                pq.emplace(dist[u], u);
            }
        }
    }
    return {dist[t] >= INF_D ? -1 : dist[t], expanded};
}

std::pair<int, int> astarWeighted(int n, const Adj& adj, const std::vector<int>& h, int s, int t) {
    const int INF_D = 1 << 29;
    std::vector<int> dist(n, INF_D);
    std::vector<char> closed(n, 0);
    dist[s] = 0;
    using Item = std::tuple<int, int, int>;          // (f, h, 節點)
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    pq.emplace(h[s], h[s], s);
    int expanded = 0;
    while (!pq.empty()) {
        Item cur = pq.top();
        pq.pop();
        int v = std::get<2>(cur);
        if (closed[v]) continue;
        closed[v] = 1;
        ++expanded;
        if (v == t) break;
        for (auto& e : adj[v]) {
            int u = e.first, w = e.second;
            if (!closed[u] && dist[v] + w < dist[u]) {
                dist[u] = dist[v] + w;
                pq.emplace(dist[u] + h[u], h[u], u);
            }
        }
    }
    return {dist[t] >= INF_D ? -1 : dist[t], expanded};
}

// ------------------------------------------------------- 輸出

std::string run(const std::vector<std::string>& grid, int sr, int sc, int tr, int tc,
                int n, const std::vector<std::array<int, 3>>& edges, int gs, int gt,
                const std::vector<int>& h, const Board& board) {
    Adj adj(n);
    for (const auto& e : edges) adj[e[0]].emplace_back(e[1], e[2]);

    std::pair<int, int> a = astarGrid(grid, sr, sc, tr, tc);
    std::pair<int, int> b = bfsGrid(grid, sr, sc, tr, tc);
    std::pair<int, int> w = {-1, 0}, d = {-1, 0};
    if (n > 0) {
        w = astarWeighted(n, adj, h, gs, gt);
        d = dijkstra(n, adj, gs, gt);
    }
    std::pair<int, long long> ida = idaStar(board);
    int od = oracleDist(board);

    int ok = 1;
    if (a.first != b.first) ok = 0;
    if (w.first != d.first) ok = 0;                 // h 可採納時兩者必須相等
    if (od >= 0 && ida.first != od) ok = 0;

    std::ostringstream oss;
    oss << a.first << " " << a.second << " " << b.first << " " << b.second << "\n";
    oss << w.first << " " << w.second << " " << d.first << " " << d.second << "\n";
    oss << ida.first << " " << ida.second << " " << od << "\n";
    oss << ok << "\n";
    return oss.str();
}

// ------------------------------------------------------- 內置測試

std::vector<std::string> genMaze(std::mt19937& rng, int R, int C) {
    std::vector<std::string> g;
    for (int r = 0; r < R; ++r) {
        std::string row;
        for (int c = 0; c < C; ++c) row.push_back(rng() % 100 < 28 ? '#' : '.');
        g.push_back(row);
    }
    return g;
}

struct GraphCase {
    int n;
    std::vector<std::array<int, 3>> edges;
    int s, t;
    std::vector<int> h;
};

GraphCase genGraph(std::mt19937& rng) {
    GraphCase out;
    out.n = 1 + static_cast<int>(rng() % 12);
    int m = static_cast<int>(rng() % 31);
    for (int i = 0; i < m; ++i) {
        int u = static_cast<int>(rng() % out.n);
        int v = static_cast<int>(rng() % out.n);
        int w = static_cast<int>(rng() % 10);
        out.edges.push_back({{u, v, w}});
    }
    out.s = static_cast<int>(rng() % out.n);
    out.t = static_cast<int>(rng() % out.n);
    const int INF_D = 1 << 29;
    std::vector<std::vector<std::pair<int, int>>> rev(out.n);
    for (auto& e : out.edges) rev[e[1]].emplace_back(e[0], e[2]);
    std::vector<int> dist(out.n, INF_D);
    dist[out.t] = 0;
    using Item = std::pair<int, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    pq.emplace(0, out.t);
    while (!pq.empty()) {
        Item cur = pq.top();
        pq.pop();
        int dd = cur.first, v = cur.second;
        if (dd > dist[v]) continue;
        for (auto& e : rev[v]) {
            if (dd + e.second < dist[e.first]) {
                dist[e.first] = dd + e.second;
                pq.emplace(dist[e.first], e.first);
            }
        }
    }
    static const int NUMS[4] = {0, 1, 3, 1};
    static const int DENS[4] = {1, 2, 4, 1};
    int pick = static_cast<int>(rng() % 4);
    for (int i = 0; i < out.n; ++i) {
        if (dist[i] >= INF_D) {
            out.h.push_back(0);
        } else {
            out.h.push_back(dist[i] * NUMS[pick] / DENS[pick]);
        }
    }
    return out;
}

void runTests() {
    // ---- 網格固定用例 ----
    {
        std::vector<std::string> g = {"...", "...", "..."};
        auto a = astarGrid(g, 0, 0, 2, 2);
        auto b = bfsGrid(g, 0, 0, 2, 2);
        assert(a.first == 4);
        assert(a.first == b.first);
        assert(a.second <= b.second);
    }
    {
        std::vector<std::string> g = {"....", ".##.", "....", "...."};
        auto a = astarGrid(g, 0, 0, 3, 3);
        auto b = bfsGrid(g, 0, 0, 3, 3);
        assert(a.first == b.first && a.first > 0);
        assert(a.second <= b.second);
    }
    {
        std::vector<std::string> g = {".....", ".###.", ".#.#.", ".###.", "....."};
        assert(astarGrid(g, 0, 0, 2, 2).first == -1);
        assert(bfsGrid(g, 0, 0, 2, 2).first == -1);
    }
    assert(astarGrid({""}, 0, 0, 0, 0).first == -1);
    assert(astarGrid({".#"}, 0, 0, 0, 1).first == -1);
    assert(astarGrid({}, 0, 0, 0, 0).first == -1);

    // ---- 帶權圖固定用例 ----
    {
        int n = 5;
        std::vector<std::array<int, 3>> edges = {
            {{0, 1, 4}}, {{0, 2, 1}}, {{2, 1, 1}}, {{1, 3, 2}}, {{2, 3, 5}}, {{3, 4, 3}}};
        Adj adj(n);
        for (auto& e : edges) adj[e[0]].emplace_back(e[1], e[2]);
        std::vector<int> h0(n, 0);
        assert(dijkstra(n, adj, 0, 4) == astarWeighted(n, adj, h0, 0, 4));
        assert(dijkstra(n, adj, 0, 4).first == 7);
        std::vector<int> htrue = {7, 5, 6, 3, 0};
        assert(astarWeighted(n, adj, htrue, 0, 4).first == 7);
        assert(astarWeighted(n, adj, htrue, 0, 4).second <= dijkstra(n, adj, 0, 4).second);
        assert(dijkstra(n, adj, 4, 0).first == -1);
    }

    // ---- 八碼難題固定用例 ----
    assert(manhattan(GOAL) == 0 && heur(GOAL) == 0);
    assert(idaStar(GOAL) == std::make_pair(0, 1LL));
    {
        Board one = {{1, 2, 3, 4, 5, 0, 7, 8, 6}};      // 空格上移一步
        assert(oracleDist(one) == 1);
        assert(idaStar(one).first == 1);
        assert(solvable(one));
    }
    {
        Board bad = {{2, 1, 3, 4, 5, 6, 7, 8, 0}};      // 逆序數 1 → 不可解
        assert(!solvable(bad));
        assert(idaStar(bad) == std::make_pair(-1, 0LL));
        assert(oracleDist(bad) == -1);
    }
    {
        Board conflict = {{3, 1, 2, 4, 5, 6, 7, 8, 0}};
        assert(manhattan(conflict) == 4);
        assert(linearConflict(conflict) == 4);
        assert(heur(conflict) == manhattan(conflict) + linearConflict(conflict));
    }

    // ---- 隨機對拍 ----
    {
        std::mt19937 rng(20261017);
        for (int it = 0; it < 300; ++it) {
            int R = 1 + static_cast<int>(rng() % 8);
            int C = 1 + static_cast<int>(rng() % 8);
            std::vector<std::string> grid = genMaze(rng, R, C);
            int sr = static_cast<int>(rng() % R), sc = static_cast<int>(rng() % C);
            int tr = static_cast<int>(rng() % R), tc = static_cast<int>(rng() % C);
            GraphCase gc = genGraph(rng);
            Board board = GOAL;
            int pos = 8;
            int walks = static_cast<int>(rng() % 10);
            for (int i = 0; i < walks; ++i) {
                const std::vector<int>& cand = MOVE_TABLE[pos];
                int np = cand[rng() % cand.size()];
                std::swap(board[pos], board[np]);
                pos = np;
            }
            std::string out = run(grid, sr, sc, tr, tc, gc.n, gc.edges, gc.s, gc.t, gc.h, board);
            assert(out.size() >= 2 && out[out.size() - 2] == '1');
        }
    }

    // ---- 深一點的八碼局面（只驗證 IDA* 與 BFS 對照組一致）----
    {
        int depths[2] = {10, 12};
        for (int di = 0; di < 2; ++di) {
            std::mt19937 rng2(static_cast<unsigned>(depths[di]));
            for (int it = 0; it < 3; ++it) {
                Board board = GOAL;
                int pos = 8;
                for (int i = 0; i < depths[di]; ++i) {
                    const std::vector<int>& cand = MOVE_TABLE[pos];
                    int np = cand[rng2() % cand.size()];
                    std::swap(board[pos], board[np]);
                    pos = np;
                }
                assert(idaStar(board).first == oracleDist(board));
            }
        }
    }
}

// ------------------------------------------------------- IO

bool isInt(const std::string& t) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i) {
        if (t[i] < '0' || t[i] > '9') return false;
    }
    return true;
}

long long clampLL(const std::string& t, long long lo, long long hi, long long def) {
    if (!isInt(t)) return def;
    errno = 0;
    char* end = nullptr;
    long long v = std::strtoll(t.c_str(), &end, 10);
    if (errno == ERANGE || end != t.c_str() + t.size()) return def;
    if (v < lo) v = lo;
    if (v > hi) v = hi;
    return v;
}

void runIo(const std::string& raw) {
    std::istringstream iss(raw);
    std::vector<std::string> toks;
    std::string t;
    while (iss >> t) toks.push_back(t);
    size_t idx = 0;
    auto nextInt = [&](long long def) -> long long {
        if (idx >= toks.size()) return def;
        std::string s = toks[idx++];
        return clampLL(s, -1000000000000000LL, 1000000000000000LL, def);
    };
    auto nextStr = [&]() -> std::string {
        if (idx >= toks.size()) return "";
        return toks[idx++];
    };

    int R = static_cast<int>(nextInt(0));
    if (R < 0) R = 0;
    if (R > MAXR) R = MAXR;
    int C = static_cast<int>(nextInt(0));
    if (C < 0) C = 0;
    if (C > MAXR) C = MAXR;
    std::vector<std::string> grid;
    for (int r = 0; r < R; ++r) {
        std::string row = nextStr();
        std::string fixed;
        for (int i = 0; i < C; ++i) {
            char ch = i < static_cast<int>(row.size()) ? row[i] : '.';   // 不足補 '.'
            if (ch != '.' && ch != '#') ch = '#';                        // 非法字元視為障礙
            fixed.push_back(ch);
        }
        grid.push_back(fixed);
    }
    if (R == 0) C = 0;

    auto cell = [&](int hi) -> int {
        long long v = nextInt(0);
        if (hi <= 0) return 0;
        if (v < 0) v = 0;
        if (v > hi - 1) v = hi - 1;
        return static_cast<int>(v);
    };
    int sr = cell(R), sc = cell(C), tr = cell(R), tc = cell(C);

    int n = static_cast<int>(nextInt(0));
    if (n < 0) n = 0;
    if (n > MAXN) n = MAXN;
    int m = static_cast<int>(nextInt(0));
    if (m < 0) m = 0;
    if (m > MAXM) m = MAXM;
    std::vector<std::array<int, 3>> edges;
    for (int i = 0; i < m; ++i) {
        long long u = nextInt(0), v = nextInt(0), w = nextInt(0);
        if (w < 0) w = 0;
        if (w > WMAX) w = WMAX;
        if (u >= 0 && u < n && v >= 0 && v < n) {
            edges.push_back({{static_cast<int>(u), static_cast<int>(v), static_cast<int>(w)}});
        }
    }
    int gs = cell(n), gt = cell(n);
    std::vector<int> h(n, 0);
    for (int i = 0; i < n; ++i) {
        long long hv = nextInt(0);
        if (hv < 0) hv = 0;
        if (hv > HMAX) hv = HMAX;
        h[i] = static_cast<int>(hv);
    }

    std::string bt = nextStr();
    Board board = GOAL;
    if (bt.size() == 9) {
        bool digits = true;
        for (char ch : bt) {
            if (ch < '0' || ch > '9') digits = false;
        }
        if (digits) {
            Board cand;
            std::vector<int> seen(9, 0);
            for (int i = 0; i < 9; ++i) {
                cand[i] = bt[i] - '0';
                seen[cand[i]] = 1;
            }
            bool perm = true;
            for (int i = 0; i < 9; ++i) {
                if (!seen[i]) perm = false;
            }
            if (perm) board = cand;
        }
    }
    std::cout << run(grid, sr, sc, tr, tc, n, edges, gs, gt, h, board);
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
    bool blank = true;
    for (char c : raw) {
        unsigned char u = static_cast<unsigned char>(c);
        if (u > ' ') { blank = false; break; }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
