// 帶權併查集：關係推斷（奇偶性 / 食物鏈）與擴展域建模
//
// 題意：
//     給定 n 個元素，元素之間的關係是一個模 d 的「循環群」：
//     d = 2 時 0 = 同類、1 = 不同類；d = 3 時 0 = 同類、1 = x 吃 y、2 = x 被 y 吃。
//     關係可傳遞：已知 rel(a,b) 與 rel(b,c) 就能推出
//         rel(a,c) = (rel(a,b) + rel(b,c)) mod d
//     依序給出 m 條陳述 (x, y, r)，表示「rel(x,y) = r」。
//     與已知資訊矛盾的陳述稱為**假話**；假話不會被記錄，不影響後續推斷。
//     統計假話的數量。
//
//     限制：d ∈ {2, 3}，1 <= n, m <= 100000，1 <= x, y <= n，0 <= r < d。
//
// 思路：
//     ### 1. 關係的代數結構
//     把每個元素看成一個「種類值」k[i] ∈ Z_d，定義
//         rel(x, y) = (k[x] - k[y]) mod d
//     於是 rel(a,c) = rel(a,b) + rel(b,c) 自動成立 —— 傳遞性免費得到。
//     一條陳述 (x, y, r) 本質上是對兩個元素的種類差給出一個約束；
//     矛盾 ⟺ 兩個約束給出了不同的差。
//
//     ### 2. 帶權併查集（主實作）
//     併查集的每個節點附一個權重 rel[x] = rel(x, parent[x])。
//     查找時做路徑壓縮，順手把權重累加成「相對於根」的
//         rel(x, root) = rel(x, parent) + rel(parent, root)
//     合併兩棵樹時，根據要維護的 rel(x, y) = r 解出新邊的權重：
//         rel(ry, rx) = (rel[x] - rel[y] - r) mod d
//     查詢時若 x, y 已連通，直接比較 (rel[x] - rel[y]) mod d 與 r。
//     兩次操作都接近 O(α(n))。
//
//     ### 3. 擴展域併查集（第二實作，測試對拍用）
//     另一個經典建模：把每個元素複製成 d 個「域節點」(i, t)
//     表示「元素 i 的種類是 t」。陳述 (x, y, r) 等價於
//         對所有 u：union(x + u·n, y + ((u - r) mod d)·n)
//     假話的判定也變得很直接：rel(x,y) = t 已知
//         ⟺ find(x + t·n) == find(y)
//     所以只要存在 t ≠ r 使得 find(x + t·n) == find(y)，就是假話。
//     x == y 且 r ≠ 0 的假話也會被正確抓到（此時 t = 0 就會命中）。
//     域節點之間永不矛盾，正是「種類唯一」的編碼。
//
//     ### 4. 三套實作互相對拍
//     測試裡同時跑 帶權併查集 / 擴展域併查集 / 樸素傳遞閉包（O(n^3)），
//     隨機數據三方結果必須一致 —— 關係取模的邊界（負數、x==y、
//     跨合併才暴露的矛盾）極易寫錯，多一套實作就多一道保險。
//
// 複雜度（α 為阿克曼函數的反函數）：
//     帶權併查集        每次操作接近 O(α(n))
//     擴展域併查集      每次操作 O(d·α(n))，此處 d <= 3
//     樸素閉包（僅測試）O(n^3) 預處理 + O(n^2) 每條陳述
//
// 應用場景：
//     - 食物鏈推斷（NOI2001 經典題）：根據觀察記錄判定哪些記錄是假的
//     - 奇偶性 / 二分圖判定：d = 2 就是「同色 / 異色」，離線判斷圖是否
//       二分圖、找出破壞二分性的矛盾邊
//     - 種類合一（unification）：Prolog / 類型推導裡「變數相等或不等」
//       的約束求解，關係群從 Z_d 推廣到更一般的置換群
//     - 版本相容性：模組 A 與 B 相容、B 與 C 不相容……傳遞推斷衝突組合
//     - 密碼學份額：Shamir 份額之間「同組 / 不同組」關係的批量校驗

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr int MMAX = 100000;

void fail(const std::string& msg) {
    std::cerr << "test failed: " << msg << std::endl;
    std::exit(1);
}

using Stmt = std::tuple<int, int, int>;

// ---------------------------------------------------------------- 帶權併查集

class WeightedDSU {
public:
    WeightedDSU(int n, int d) : n_(n), d_(d) {
        parent_.resize(n + 1);
        size_.assign(n + 1, 1);
        rel_.assign(n + 1, 0);
        for (int i = 0; i <= n; ++i) {
            parent_[i] = i;
        }
    }

    int find(int x) {
        // 第一遍：找根
        int root = x;
        while (parent_[root] != root) {
            root = parent_[root];
        }
        // 第二遍：自頂向下壓縮，權重累加為相對於根
        std::vector<int> path;
        int cur = x;
        while (parent_[cur] != root) {
            path.push_back(cur);
            cur = parent_[cur];
        }
        for (size_t i = path.size(); i-- > 0;) {
            int node = path[i];
            int p = parent_[node];
            rel_[node] = (rel_[node] + rel_[p]) % d_;
            parent_[node] = root;
        }
        return root;
    }

    void unite(int x, int y, int r) {
        // 記錄一條確定為真的陳述 rel(x,y) = r
        int rx = find(x);
        int ry = find(y);
        if (rx == ry) {
            return;
        }
        int dx = rel_[x];  // rel(x, rx)
        int dy = rel_[y];  // rel(y, ry)
        if (size_[rx] < size_[ry]) {
            // 把 rx 掛到 ry 下面：rel(rx, ry) = (r - dx + dy) mod d
            parent_[rx] = ry;
            rel_[rx] = ((r - dx + dy) % d_ + d_) % d_;
            size_[ry] += size_[rx];
        } else {
            // 把 ry 掛到 rx 下面：rel(ry, rx) = (dx - dy - r) mod d
            parent_[ry] = rx;
            rel_[ry] = ((dx - dy - r) % d_ + d_) % d_;
            size_[rx] += size_[ry];
        }
    }

    std::optional<int> relation(int x, int y) {
        // x, y 已連通時回傳 rel(x,y)，否則回 nullopt
        int rx = find(x);
        int ry = find(y);
        if (rx != ry) {
            return std::nullopt;
        }
        return ((rel_[x] - rel_[y]) % d_ + d_) % d_;
    }

private:
    int n_, d_;
    std::vector<int> parent_, size_, rel_;
};

// ---------------------------------------------------------------- 擴展域併查集

class ExtendedDSU {
public:
    ExtendedDSU(int n, int d) : n_(n), d_(d) {
        parent_.resize(n * d + 1);
        size_.assign(n * d + 1, 1);
        for (int i = 0; i <= n * d; ++i) {
            parent_[i] = i;
        }
    }

    int find(int x) {
        int root = x;
        while (parent_[root] != root) {
            root = parent_[root];
        }
        while (parent_[x] != root) {
            int nxt = parent_[x];
            parent_[x] = root;
            x = nxt;
        }
        return root;
    }

    void unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) {
            return;
        }
        if (size_[ra] < size_[rb]) {
            std::swap(ra, rb);
        }
        parent_[rb] = ra;
        size_[ra] += size_[rb];
    }

    bool addStatement(int x, int y, int r) {
        // 回傳 true 表示為真並已記錄，false 表示是假話
        for (int t = 0; t < d_; ++t) {
            if (t != r && find(x + t * n_) == find(y)) {
                return false;  // 已知關係是 t 而非 r
            }
        }
        for (int u = 0; u < d_; ++u) {
            unite(x + u * n_, y + ((u - r) % d_ + d_) % d_ * n_);
        }
        return true;
    }

private:
    int n_, d_;
    std::vector<int> parent_, size_;
};

// ---------------------------------------------------------------- 樸素閉包（測試用）

bool closeClosure(std::vector<std::vector<std::optional<int>>>& R, int n, int d) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (int a = 1; a <= n; ++a) {
            for (int b = 1; b <= n; ++b) {
                if (!R[a][b]) {
                    continue;
                }
                for (int c = 1; c <= n; ++c) {
                    if (!R[b][c]) {
                        continue;
                    }
                    int v = (*R[a][b] + *R[b][c]) % d;
                    if (!R[a][c]) {
                        R[a][c] = v;
                        changed = true;
                    } else if (*R[a][c] != v) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

int solveBrute(int n, int d, const std::vector<Stmt>& stmts) {
    std::vector<std::vector<std::optional<int>>> R(
        n + 1, std::vector<std::optional<int>>(n + 1, std::nullopt));
    for (int i = 1; i <= n; ++i) {
        R[i][i] = 0;
    }
    int falseCnt = 0;
    for (const Stmt& s : stmts) {
        int x = std::get<0>(s), y = std::get<1>(s), r = std::get<2>(s);
        if (R[x][y]) {
            if (*R[x][y] != r) {
                ++falseCnt;
            }
            continue;
        }
        R[x][y] = r;
        R[y][x] = (d - r) % d;
        if (!closeClosure(R, n, d)) {
            ++falseCnt;
            R[x][y] = std::nullopt;
            R[y][x] = std::nullopt;
        }
    }
    return falseCnt;
}

// ---------------------------------------------------------------- 三套 solve

int solveWeighted(int n, int d, const std::vector<Stmt>& stmts) {
    WeightedDSU dsu(n, d);
    int falseCnt = 0;
    for (const Stmt& s : stmts) {
        int x = std::get<0>(s), y = std::get<1>(s), r = std::get<2>(s);
        std::optional<int> known = dsu.relation(x, y);
        if (!known) {
            dsu.unite(x, y, r);
        } else if (*known != r) {
            ++falseCnt;
        }
    }
    return falseCnt;
}

int solveExtended(int n, int d, const std::vector<Stmt>& stmts) {
    ExtendedDSU dsu(n, d);
    int falseCnt = 0;
    for (const Stmt& s : stmts) {
        int x = std::get<0>(s), y = std::get<1>(s), r = std::get<2>(s);
        if (!dsu.addStatement(x, y, r)) {
            ++falseCnt;
        }
    }
    return falseCnt;
}

int solve(int n, int d, const std::vector<Stmt>& stmts) {
    return solveWeighted(n, d, stmts);
}

// ---------------------------------------------------------------- 測試

void runTests() {
    auto eq = [](int got, int want, const std::string& tag) {
        if (got != want) {
            fail(tag + ": got " + std::to_string(got) + " want " +
                 std::to_string(want));
        }
    };

    // ---- 固定用例（手工驗算）----
    // d=3 食物鏈：第 3、5 條是假話
    std::vector<Stmt> stmts = {{1, 2, 1}, {2, 3, 1}, {3, 1, 2},
                               {1, 1, 0}, {1, 1, 1}, {1, 3, 2}};
    eq(solve(3, 3, stmts), 2, "fixed food d3");
    eq(solveExtended(3, 3, stmts), 2, "fixed food d3 ext");
    eq(solveBrute(3, 3, stmts), 2, "fixed food d3 brute");
    // d=2 奇偶
    std::vector<Stmt> stmts2 = {{1, 2, 1}, {2, 3, 1}, {1, 3, 0},
                                {1, 4, 1}, {4, 2, 1}, {3, 4, 1}};
    eq(solve(4, 2, stmts2), 1, "fixed parity d2");
    eq(solveExtended(4, 2, stmts2), 1, "fixed parity d2 ext");
    eq(solveBrute(4, 2, stmts2), 1, "fixed parity d2 brute");
    // d=2 奇環
    std::vector<Stmt> stmts3 = {{1, 2, 1}, {2, 3, 1}, {1, 3, 1}};
    eq(solve(3, 2, stmts3), 1, "fixed odd cycle");
    eq(solveExtended(3, 2, stmts3), 1, "fixed odd cycle ext");
    eq(solveBrute(3, 2, stmts3), 1, "fixed odd cycle brute");
    // x == x 的各種情形
    eq(solve(1, 2, {{1, 1, 0}, {1, 1, 1}}), 1, "fixed self d2");
    eq(solve(1, 3, {{1, 1, 0}, {1, 1, 1}, {1, 1, 2}}), 2, "fixed self d3");
    eq(solveExtended(1, 3, {{1, 1, 0}, {1, 1, 1}, {1, 1, 2}}), 2,
       "fixed self d3 ext");
    // 空語句
    eq(solve(5, 3, {}), 0, "fixed empty");
    // 跨合併才暴露的矛盾
    std::vector<Stmt> stmts4 = {{1, 2, 1}, {3, 4, 1}, {2, 3, 0}, {1, 4, 0}};
    eq(solve(4, 3, stmts4), 1, "fixed late merge");
    eq(solveExtended(4, 3, stmts4), 1, "fixed late merge ext");
    eq(solveBrute(4, 3, stmts4), 1, "fixed late merge brute");

    // ---- 隨機對拍：三套實作結果必須一致 ----
    std::mt19937 rng(20261026);
    for (int tc = 0; tc < 300; ++tc) {
        int d = (rng() % 2 == 0) ? 2 : 3;
        int n = static_cast<int>(rng() % 10) + 1;
        int m = static_cast<int>(rng() % 40) + 1;
        std::vector<Stmt> stmts;
        for (int i = 0; i < m; ++i) {
            int x = static_cast<int>(rng() % n) + 1;
            int y = static_cast<int>(rng() % n) + 1;
            int r = static_cast<int>(rng() % d);
            stmts.emplace_back(x, y, r);
        }
        int a = solveWeighted(n, d, stmts);
        int b = solveExtended(n, d, stmts);
        int c = solveBrute(n, d, stmts);
        if (a != b || b != c) {
            fail("random triple mismatch: d=" + std::to_string(d) +
                 " n=" + std::to_string(n) + " a=" + std::to_string(a) +
                 " b=" + std::to_string(b) + " c=" + std::to_string(c));
        }
    }

    // ---- 稍大規模：帶權 vs 擴展域 ----
    std::mt19937 big(13572468);
    int n = 3000, m = 6000;
    std::vector<Stmt> stmtsBig;
    for (int i = 0; i < m; ++i) {
        int x = static_cast<int>(big() % n) + 1;
        int y = static_cast<int>(big() % n) + 1;
        int r = static_cast<int>(big() % 3);
        stmtsBig.emplace_back(x, y, r);
    }
    eq(solveWeighted(n, 3, stmtsBig), solveExtended(n, 3, stmtsBig),
       "random large");
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& s) {
    // 與 Python 版一致：可選正負號 + 純數字，空串/非數字/越界皆 nullopt
    if (s.empty()) {
        return std::nullopt;
    }
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') {
        i = 1;
    }
    if (i == s.size()) {
        return std::nullopt;
    }
    for (size_t j = i; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') {
            return std::nullopt;
        }
    }
    long long val = 0;
    for (size_t j = i; j < s.size(); ++j) {
        val = val * 10 + (s[j] - '0');
        if (val > 1000000000000000000LL) {
            return std::nullopt;
        }
    }
    if (s[0] == '-') {
        val = -val;
    }
    return val;
}

void runIo(const std::string& raw) {
    std::istringstream in(raw);
    std::vector<std::string> toks;
    std::string t;
    while (in >> t) {
        toks.push_back(t);
    }
    size_t pos = 0;
    auto nxt = [&]() -> const std::string* {
        if (pos >= toks.size()) {
            return nullptr;
        }
        return &toks[pos++];
    };

    auto bad = []() {
        std::cout << "error\n";
    };

    const std::string* p = nxt();
    std::optional<long long> d = p ? parseInt(*p) : std::nullopt;
    if (!d || (*d != 2 && *d != 3)) {
        bad();
        return;
    }
    p = nxt();
    std::optional<long long> n = p ? parseInt(*p) : std::nullopt;
    if (!n || !(*n >= 1 && *n <= NMAX)) {
        bad();
        return;
    }
    p = nxt();
    std::optional<long long> m = p ? parseInt(*p) : std::nullopt;
    if (!m || !(*m >= 1 && *m <= MMAX)) {
        bad();
        return;
    }
    std::vector<Stmt> stmts;
    for (int i = 0; i < *m; ++i) {
        p = nxt();
        std::optional<long long> xo = p ? parseInt(*p) : std::nullopt;
        p = nxt();
        std::optional<long long> yo = p ? parseInt(*p) : std::nullopt;
        p = nxt();
        std::optional<long long> ro = p ? parseInt(*p) : std::nullopt;
        // 解析失敗與越界最終都導向同一結果（error），先用值取代 optional
        long long x = xo.value_or(-1);
        long long y = yo.value_or(-1);
        long long r = ro.value_or(-1);
        if (x < 1 || x > *n || y < 1 || y > *n || r < 0 || r >= *d) {
            bad();
            return;
        }
        stmts.emplace_back(static_cast<int>(x), static_cast<int>(y),
                           static_cast<int>(r));
    }
    std::cout << solve(static_cast<int>(*n), static_cast<int>(*d), stmts)
              << "\n";
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (size_t i = 0; i < raw.size(); ++i) {
        unsigned char u = static_cast<unsigned char>(raw[i]);
        if (u > ' ') {
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
