// 線段樹合併（Segment Tree Merge）：樹上子樹的相異值統計
//
// 題意：
//     給定一棵有 n 個節點的**根樹**（根為 1），每個節點 u 有一個整數標籤 c[u]。
//     對每個節點 u，求其子樹中出現過的**相異**（distinct）標籤數量。
//     限制：1 <= n <= 100000，|c[u]| <= 10^9。
//
// 思路：
//     ### 1. 為什麼不能用樸素做法
//     最直接的想法是「對每個節點都遍歷一次子樹，把標籤塞進 set」。
//     但一個鏈狀樹上，子樹大小總和是 O(n^2)，直接超時。
//     改進一點：把小集合併進大集合（启发式合并 / small-to-large），
//     可以做到 O(n log^2 n)，但仍然要維護 hash 集合，常數不小。
//
//     ### 2. 線段樹合併的核心想法
//     換一個視角：**自底向上，把子樹的資訊「合併」到父節點**。
//     對每個節點 u 維護一棵「值域線段樹」，代表 u 子樹裡出現過的所有標籤：
//         葉子（某個具體標籤值）= 該標籤是否出現（0/1）
//         內部節點 = 左右兩邊的出現數之和
//     於是 u 子樹的相異標籤數 = u 那棵線段樹的根節點的 cnt。
//
//     合併兩棵線段樹 merge(a, b) 的規則非常簡單：
//         若 a 為空 → 直接回傳 b
//         若 b 為空 → 直接回傳 a
//         否則（兩者都非空，這是最後一種情形也是唯一需要遞迴的情形）：
//             a 的左子樹 = merge(a.left, b.left)
//             a 的右子樹 = merge(a.right, b.right)
//             a.cnt = a.left.cnt + a.right.cnt
//     底層標籤是「出現/不出現」，所以葉子合併時取 or（1）而不是相加 ——
//     這正是「相異」的語意。
//
//     ### 3. 為什麼總複雜度是 O(n log V) 而不是 O(n^2)
//     merge 裡唯一會「消滅」節點的情形是「a、b 都非空」；
//     此時 b 的根節點被丟棄（不再被任何線段樹引用）。
//     每次 merge 的「有消耗」呼叫都會消滅一個節點，
//     而全程建立的節點總數不超過 O(n log V)（每個節點插入一個標籤，
//     最多開 log V 個節點）。
//     所以所有 merge 呼叫加總是 O(節點總數) = O(n log V)。
//     這是線段樹合併最漂亮的地方：**它不需要任何平衡技巧**。
//
//     ### 4. 實作順序（後序遍歷）
//     用顯式堆疊做 DFS 取得「父在子前」的順序，再倒過來就是「子在父前」。
//     對每個節點 u（依此順序）：
//         1. 把 u 自己的標籤插入 root[u]（此時 root[u] 已含所有子樹的貢獻）
//         2. ans[u] = cnt[root[u]]
//         3. 把 root[u] 合併進 root[parent[u]]
//     注意 ans[u] 必須在「合併進父親」之前記錄，
//     因為合併會改寫（甚至換掉）那些節點。
//
// 複雜度：
//     離散化            O(n log n)
//     建立 + 插入       O(n log V)
//     所有 merge 合計    O(n log V)（節點消耗 amortize）
//     總計              時間 O(n log n)，空間 O(n log V)
//
// 應用場景（Applications）：
//     - **專案成員技能盤點**：每個節點是一個人，子樹是團隊，
//       「團隊涵蓋了幾種技能」就是子樹的相異值數
//     - **物種普查**：把觀測點依地理層級組成樹，
//       每個區域的物種豐富度（species richness）就是這個問題
//     - **授權 / 權限盤點**：群組樹上的「有效權限種類數」
//     - **程式碼重複偵測**：以 AST 節點組成樹，統計子樹內的相異節點型別數

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr long long VAL_MAX = 1000000000LL;
constexpr long long PARSE_LIMIT = 1000000000000000000LL;

// 動態開點的線段樹節點池（平行陣列，避免大量物件配置）
class SegmentTreePool {
public:
    SegmentTreePool() {
        left_.push_back(0);
        right_.push_back(0);
        cnt_.push_back(0);
    }

    int newNode() {
        left_.push_back(0);
        right_.push_back(0);
        cnt_.push_back(0);
        return static_cast<int>(cnt_.size()) - 1;
    }

    // 在 pos 處插入標籤（存在則維持 1，即集合語意）
    int insert(int node, int pos, int lo, int hi) {
        if (node == 0) node = newNode();
        if (lo == hi) {
            cnt_[node] = 1;
            return node;
        }
        int mid = (lo + hi) >> 1;
        if (pos <= mid) {
            left_[node] = insert(left_[node], pos, lo, mid);
        } else {
            right_[node] = insert(right_[node], pos, mid + 1, hi);
        }
        cnt_[node] = cnt_[left_[node]] + cnt_[right_[node]];
        return node;
    }

    // 合併兩棵線段樹（集合聯集）
    int merge(int a, int b, int lo, int hi) {
        if (a == 0) return b;
        if (b == 0) return a;
        if (lo == hi) {
            // 同一個標籤值：出現即為 1（相異，不是計數）
            cnt_[a] = 1;
            return a;
        }
        int mid = (lo + hi) >> 1;
        left_[a] = merge(left_[a], left_[b], lo, mid);
        right_[a] = merge(right_[a], right_[b], mid + 1, hi);
        cnt_[a] = cnt_[left_[a]] + cnt_[right_[a]];
        return a;
    }

    int count(int node) const { return node ? cnt_[node] : 0; }

private:
    std::vector<int> left_;
    std::vector<int> right_;
    std::vector<int> cnt_;
};

std::vector<int> subtreeDistinctCounts(int n, const std::vector<long long>& labels,
                                        const std::vector<std::pair<int, int>>& edges) {
    // ---- 1. 標籤離散化 ----
    std::vector<long long> uniq = labels;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    int m = static_cast<int>(uniq.size());
    std::unordered_map<long long, int> rankOf;
    rankOf.reserve(static_cast<std::size_t>(m) * 2 + 1);
    for (int i = 0; i < m; ++i) rankOf[uniq[i]] = i + 1;  // 1-based 值域

    // ---- 2. 建鄰接表 ----
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n) + 1);
    for (const auto& e : edges) {
        adj[static_cast<std::size_t>(e.first)].push_back(e.second);
        adj[static_cast<std::size_t>(e.second)].push_back(e.first);
    }

    // ---- 3. 取得「父在子前」的遍歷順序（顯式堆疊，避免遞迴爆栈）----
    std::vector<int> parent(static_cast<std::size_t>(n) + 1, 0);
    std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
    std::vector<int> order;
    order.reserve(n);
    std::vector<int> stack;
    stack.push_back(1);
    seen[1] = 1;
    parent[1] = -1;
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        order.push_back(u);
        for (int v : adj[static_cast<std::size_t>(u)]) {
            if (!seen[static_cast<std::size_t>(v)]) {
                seen[static_cast<std::size_t>(v)] = 1;
                parent[static_cast<std::size_t>(v)] = u;
                stack.push_back(v);
            }
        }
    }

    // ---- 4. 自底向上：插入 + 記錄答案 + 合併到父親 ----
    SegmentTreePool pool;
    std::vector<int> root(static_cast<std::size_t>(n) + 1, 0);
    std::vector<int> ans(static_cast<std::size_t>(n) + 1, 0);
    for (int i = static_cast<int>(order.size()) - 1; i >= 0; --i) {
        int u = order[static_cast<std::size_t>(i)];
        int r = pool.insert(root[static_cast<std::size_t>(u)],
                            rankOf[labels[static_cast<std::size_t>(u - 1)]], 1, m);
        root[static_cast<std::size_t>(u)] = r;
        ans[static_cast<std::size_t>(u)] = pool.count(r);
        int p = parent[static_cast<std::size_t>(u)];
        if (p > 0) {
            root[static_cast<std::size_t>(p)] =
                pool.merge(root[static_cast<std::size_t>(p)], r, 1, m);
        }
    }
    return std::vector<int>(ans.begin() + 1, ans.end());
}

// 參考實作：對每個節點遍歷子樹收集標籤去重（O(n^2)，僅供對拍）
std::vector<int> bruteForce(int n, const std::vector<long long>& labels,
                             const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n) + 1);
    for (const auto& e : edges) {
        adj[static_cast<std::size_t>(e.first)].push_back(e.second);
        adj[static_cast<std::size_t>(e.second)].push_back(e.first);
    }
    // 先求出以 1 為根的父子關係，樸素遍歷只能往下走（不能回頭走到父親）
    std::vector<int> parent(static_cast<std::size_t>(n) + 1, 0);
    std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
    std::vector<int> stack;
    stack.push_back(1);
    seen[1] = 1;
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (int v : adj[static_cast<std::size_t>(u)]) {
            if (!seen[static_cast<std::size_t>(v)]) {
                seen[static_cast<std::size_t>(v)] = 1;
                parent[static_cast<std::size_t>(v)] = u;
                stack.push_back(v);
            }
        }
    }
    std::vector<std::vector<int>> children(static_cast<std::size_t>(n) + 1);
    for (int v = 2; v <= n; ++v) {
        children[static_cast<std::size_t>(parent[static_cast<std::size_t>(v)])].push_back(v);
    }

    std::vector<int> ans(static_cast<std::size_t>(n) + 1, 0);
    std::vector<long long> vals;
    for (int start = 1; start <= n; ++start) {
        vals.clear();
        stack.clear();
        stack.push_back(start);
        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            vals.push_back(labels[static_cast<std::size_t>(u - 1)]);
            for (int v : children[static_cast<std::size_t>(u)]) stack.push_back(v);
        }
        std::sort(vals.begin(), vals.end());
        ans[static_cast<std::size_t>(start)] =
            static_cast<int>(std::unique(vals.begin(), vals.end()) - vals.begin());
    }
    return std::vector<int>(ans.begin() + 1, ans.end());
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

void checkEqual(const std::vector<int>& got, const std::vector<int>& want,
                const char* msg) {
    if (got != want) fail(msg);
}

// 應用程式示範：把「資源池合併」抽象成同一套 merge 機制
void demoResourcePoolMerge() {
    SegmentTreePool pool;
    std::vector<std::pair<int, std::vector<int>>> sets = {
        {1, {1, 5, 9}}, {2, {5, 7}}, {3, {9, 11}}};
    std::vector<int> roots;
    for (const auto& s : sets) {
        int r = 0;
        for (int v : s.second) r = pool.insert(r, v, 1, 16);
        roots.push_back(r);
    }
    if (pool.count(roots[0]) != 3) fail("demo a failed");
    if (pool.count(roots[1]) != 2) fail("demo b failed");
    if (pool.count(roots[2]) != 2) fail("demo c failed");

    int mergedAb = pool.merge(roots[0], roots[1], 1, 16);
    if (pool.count(mergedAb) != 4) fail("demo d failed");  // {1,5,9,7}
    int mergedAbc = pool.merge(mergedAb, roots[2], 1, 16);
    if (pool.count(mergedAbc) != 5) fail("demo e failed");
    int again = pool.merge(mergedAbc, mergedAbc, 1, 16);
    if (pool.count(again) != 5) fail("demo f failed");
}

void runTests() {
    // ---- 固定用例 ----
    {
        std::vector<long long> init{7};
        std::vector<int> want{1};
        checkEqual(subtreeDistinctCounts(1, init, {}), want, "case 1 failed");
    }
    {
        std::vector<long long> init{1, 2, 3};
        std::vector<std::pair<int, int>> edges{{1, 2}, {2, 3}};
        std::vector<int> want{3, 2, 1};
        checkEqual(subtreeDistinctCounts(3, init, edges), want, "case 2 failed");
    }
    {
        std::vector<long long> init{5, 5, 5};
        std::vector<std::pair<int, int>> edges{{1, 2}, {2, 3}};
        std::vector<int> want{1, 1, 1};
        checkEqual(subtreeDistinctCounts(3, init, edges), want, "case 3 failed");
    }
    {
        std::vector<long long> init{9, 1, 2, 1};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3}, {1, 4}};
        std::vector<int> want{3, 1, 1, 1};
        checkEqual(subtreeDistinctCounts(4, init, edges), want, "case 4 failed");
    }
    {
        std::vector<long long> init{1, 2, 3, 1, 4};
        std::vector<std::pair<int, int>> edges{{1, 2}, {2, 4}, {1, 3}, {3, 5}};
        std::vector<int> want{4, 2, 2, 1, 1};
        checkEqual(subtreeDistinctCounts(5, init, edges), want, "case 5 failed");
    }
    {
        std::vector<long long> init{-1, 0, -1};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3}};
        std::vector<int> want{2, 1, 1};
        checkEqual(subtreeDistinctCounts(3, init, edges), want, "case 6 failed");
    }
    {
        std::vector<long long> init{8, 8, 8, 8};
        std::vector<std::pair<int, int>> edges{{1, 2}, {2, 3}, {3, 4}};
        std::vector<int> want{1, 1, 1, 1};
        checkEqual(subtreeDistinctCounts(4, init, edges), want, "case 7 failed");
    }
    {
        std::vector<long long> init{1, 2, 3, 4, 2, 3, 4};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3}, {2, 4},
                                               {2, 5}, {3, 6}, {3, 7}};
        std::vector<int> want{4, 2, 2, 1, 1, 1, 1};
        checkEqual(subtreeDistinctCounts(7, init, edges), want, "case 8 failed");
    }
    // 9) 應用程式示範
    demoResourcePoolMerge();

    // ---- 與參考實作對拍 ----
    {
        Lcg rng(20261102ULL);
        for (int iter = 0; iter < 200; ++iter) {
            int n = rng.range(1, 10);
            std::vector<long long> labels;
            for (int i = 0; i < n; ++i) labels.push_back(rng.range(-3, 3));
            std::vector<std::pair<int, int>> edges;
            for (int v = 2; v <= n; ++v) edges.push_back({rng.range(1, v - 1), v});
            checkEqual(subtreeDistinctCounts(n, labels, edges),
                       bruteForce(n, labels, edges), "random mismatch");
        }
    }
    {
        Lcg rng(99887766ULL);
        for (int iter = 0; iter < 150; ++iter) {
            int n = rng.range(1, 12);
            std::vector<long long> labels;
            for (int i = 0; i < n; ++i) labels.push_back(rng.range(-1000, 1000));
            std::vector<std::pair<int, int>> edges;
            for (int v = 2; v <= n; ++v) edges.push_back({rng.range(1, v - 1), v});
            checkEqual(subtreeDistinctCounts(n, labels, edges),
                       bruteForce(n, labels, edges), "mid random mismatch");
        }
    }
    {
        Lcg rng(13572468ULL);
        for (int iter = 0; iter < 150; ++iter) {
            int n = rng.range(1, 12);
            std::vector<long long> labels;
            for (int i = 0; i < n; ++i) labels.push_back(rng.range(-2, 2));
            std::vector<std::pair<int, int>> edges;
            for (int v = 2; v <= n; ++v) edges.push_back({rng.range(1, v - 1), v});
            // 打亂邊的順序（根仍然由 1 沿父方向決定）
            for (std::size_t i = edges.size(); i > 1; --i) {
                std::size_t j = static_cast<std::size_t>(rng.range(1, static_cast<int>(i)));
                std::swap(edges[i - 1], edges[j - 1]);
            }
            checkEqual(subtreeDistinctCounts(n, labels, edges),
                       bruteForce(n, labels, edges), "shuffled random mismatch");
        }
    }
    {
        // 較大規模：深鏈 + 大量重複
        Lcg rng(24681357ULL);
        int n = 400;
        std::vector<long long> labels;
        for (int i = 0; i < n; ++i) labels.push_back(rng.range(0, 5));
        std::vector<std::pair<int, int>> edges;
        for (int i = 1; i < n; ++i) edges.push_back({i, i + 1});
        checkEqual(subtreeDistinctCounts(n, labels, edges), bruteForce(n, labels, edges),
                   "large chain mismatch");
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
    std::optional<long long> nOpt = std::nullopt;
    if (pt) nOpt = parseInt(*pt);
    if (!nOpt || *nOpt < 1 || *nOpt > NMAX) {
        bad();
        return;
    }
    int n = static_cast<int>(*nOpt);

    std::vector<long long> labels;
    labels.reserve(n);
    for (int i = 0; i < n; ++i) {
        const std::string* t = nxt();
        std::optional<long long> v = std::nullopt;
        if (t) v = parseInt(*t);
        if (!v || *v < -VAL_MAX || *v > VAL_MAX) {
            bad();
            return;
        }
        labels.push_back(*v);
    }

    std::vector<std::pair<int, int>> edges;
    edges.reserve(n > 0 ? n - 1 : 0);
    for (int i = 0; i + 1 < n; ++i) {
        const std::string* tu = nxt();
        const std::string* tv = nxt();
        std::optional<long long> u = std::nullopt;
        std::optional<long long> v = std::nullopt;
        if (tu) u = parseInt(*tu);
        if (tv) v = parseInt(*tv);
        if (!u || !v || *u < 1 || *u > n || *v < 1 || *v > n) {
            bad();
            return;
        }
        edges.push_back({static_cast<int>(*u), static_cast<int>(*v)});
    }

    std::vector<int> answers = subtreeDistinctCounts(n, labels, edges);
    for (int a : answers) std::cout << a << "\n";
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
