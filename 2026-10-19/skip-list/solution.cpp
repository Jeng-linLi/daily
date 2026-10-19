// 跳表（Skip List）：以隨機層高達成期望 O(log n) 的有序多重集合
//
// 題意：
//     維護一個**有序多重集合**（允許重複元素），支援以下操作，並與「有序陣列 + 二分」
//     的參考實作逐項對拍，確保跳表的行為完全正確。
//
//     操作碼（每行 `op x`）：
//         1  insert x    插入 x（允許重複）
//         2  erase  x    刪除一個 x（不存在則忽略）
//         3  count  x    x 的出現次數
//         4  rank   x    嚴格小於 x 的元素個數
//         5  kth    k    第 k 小的元素（k 從 1 起算），越界回傳 -1
//         6  pred   x    嚴格小於 x 的最大元素，不存在回傳 -1
//         7  succ   x    嚴格大於 x 的最小元素，不存在回傳 -1
//
// 思路：
//     ### 1. 為什麼需要跳表
//     有序陣列的插入/刪除是 O(n)；平衡樹（紅黑樹、AVL）能做到 O(log n) 但旋轉規則複雜。
//     跳表用「**機率平衡**」換來簡潔：每個節點隨機一個層高 lv（P(lv ≥ k) = 2^{-(k-1)}），
//     第 i 層的指標串起一條「稀疏高速公路」，越高層跳過越多節點。
//     搜尋時從最高層往右走，走不過去就下降一層，期望步數 O(log n)。
//
//     ### 2. span（跨度）欄位：讓「排名」也能 O(log n)
//     純指標只能回答「在不在」，要回答 rank / kth 必須知道每條邊**跨過幾個節點**。
//     每個節點在第 i 層額外存 span[i]：從自己沿第 i 層指標走一步，等於在第 0 層往前跳幾個節點。
//     插入時若節點只出現在 0..lv-1 層，更高層的每條邊都被「多跨過一個節點」：
//         - i <  lv：新節點切斷舊邊，兩段跨度是 (rank[0]-rank[i]) 與 舊跨度 - (rank[0]-rank[i])
//         - i ≥  lv：舊邊仍然存在，跨度直接 +1
//     刪除則是逆向操作：update[i].span[i] += x->span[i] - 1。
//
//     ### 3. 隨機層高必須跨語言一致
//     Python 的 random 與 C++ 的 mt19937 序列不同，兩邊的樹形會長得不一樣。
//     這裡兩個版本都用**同一組 LCG 常數**自製亂數，取 64 位元的第 33 位當硬幣：
//         state = state * 6364136223846793005 + 1442695040888963407  (mod 2^64)
//     擲到 1 就停止，否則層高 +1（上限 32 層）。
//
//     ### 4. 不變量檢查
//     對每一層的每一條邊，都用第 0 層指標實際走一趟，確認 span[i] 等於真實距離。
//
// 應用場景：
//     - Redis 的 zset（有序集合）底層就是跳表（zskiplist），ZRANK / ZRANGE 靠 span 實作
//     - LevelDB / RocksDB 的 MemTable（併發友善的記憶體索引）
//     - 需要細粒度區間鎖、或想要無鎖（lock-free）實作的場合，跳表比平衡樹容易做
//     - 任何「需要有序 + 頻繁插入刪除 + 排名查詢」的中介層（排行榜、時間軸、價格檔位）
//
// 複雜度：
//     搜尋 / 插入 / 刪除              期望 O(log n) 時間
//     rank / kth / pred / succ        期望 O(log n) 時間（靠 span 欄位）
//     count                          期望 O(log n) 時間（count_leq - count_less）
//     空間                           期望 O(n)（每個節點期望 2 條邊），最壞 O(n log n)

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const int MAX_LEVEL = 32;
const unsigned long long LCG_MUL = 6364136223846793005ULL;
const unsigned long long LCG_ADD = 1442695040888963407ULL;
const unsigned long long INIT_SEED = 0x9E3779B97F4A7C15ULL;

const int QMAX = 300;               // IO 模式下的操作數上限
const long long VMAX = 1000000;     // 元素值上限（夾到 [0, VMAX]，讓 -1 可以當哨兵）

// ---------------------------------------------------------------- 跳表節點

struct Node {
    long long key;
    std::vector<Node*> forward;
    std::vector<long long> span;
    explicit Node(long long k, int level) : key(k), forward(level, nullptr), span(level, 0) {}
};

class SkipList {
public:
    explicit SkipList(unsigned long long seed = INIT_SEED)
        : state_(seed), head_(new Node(0, MAX_LEVEL)), level_(1), length_(0) {}

    long long size() const { return length_; }

    void insert(long long key) {
        std::vector<Node*> update(MAX_LEVEL, nullptr);
        std::vector<long long> rank(MAX_LEVEL, 0);
        findPreds(key, true, update, rank);

        int lv = randomLevel();
        if (lv > level_) {
            for (int i = level_; i < lv; ++i) {
                update[i] = head_;
                update[i]->span[i] = length_;
                rank[i] = 0;
            }
            level_ = lv;
        }
        Node* node = new Node(key, lv);
        for (int i = 0; i < lv; ++i) {
            node->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = node;
            node->span[i] = update[i]->span[i] - (rank[0] - rank[i]);
            update[i]->span[i] = (rank[0] - rank[i]) + 1;
        }
        for (int i = lv; i < level_; ++i) {
            update[i]->span[i] += 1;
        }
        ++length_;
    }

    bool erase(long long key) {
        std::vector<Node*> update(MAX_LEVEL, nullptr);
        std::vector<long long> rank(MAX_LEVEL, 0);
        findPreds(key, true, update, rank);

        Node* x = update[0]->forward[0];
        if (x == nullptr || x->key != key) return false;
        for (int i = 0; i < level_; ++i) {
            if (update[i]->forward[i] == x) {
                update[i]->span[i] += x->span[i] - 1;
                update[i]->forward[i] = x->forward[i];
            } else {
                update[i]->span[i] -= 1;
            }
        }
        while (level_ > 1 && head_->forward[level_ - 1] == nullptr) --level_;
        --length_;
        return true;
    }

    long long countLess(long long key) const {
        long long cnt = 0;
        Node* x = head_;
        for (int i = level_ - 1; i >= 0; --i) {
            while (true) {
                Node* nxt = x->forward[i];
                if (nxt == nullptr || nxt->key >= key) break;
                cnt += x->span[i];
                x = nxt;
            }
        }
        return cnt;
    }

    long long countLeq(long long key) const {
        long long cnt = 0;
        Node* x = head_;
        for (int i = level_ - 1; i >= 0; --i) {
            while (true) {
                Node* nxt = x->forward[i];
                if (nxt == nullptr || nxt->key > key) break;
                cnt += x->span[i];
                x = nxt;
            }
        }
        return cnt;
    }

    long long count(long long key) const { return countLeq(key) - countLess(key); }

    long long kth(long long k) const {
        if (k < 1 || k > length_) return -1;
        long long traversed = 0;
        Node* x = head_;
        for (int i = level_ - 1; i >= 0; --i) {
            while (true) {
                Node* nxt = x->forward[i];
                if (nxt == nullptr || traversed + x->span[i] > k) break;
                traversed += x->span[i];
                x = nxt;
            }
        }
        if (x == head_) return -1;
        return x->key;
    }

    long long predecessor(long long key) const {
        long long r = countLess(key);
        return r >= 1 ? kth(r) : -1;
    }

    long long successor(long long key) const { return kth(countLeq(key) + 1); }

    std::vector<long long> toList() const {
        std::vector<long long> out;
        Node* x = head_->forward[0];
        while (x != nullptr) {
            out.push_back(x->key);
            x = x->forward[0];
        }
        return out;
    }

    std::vector<long long> levels() const {
        std::vector<long long> out;
        Node* x = head_->forward[0];
        while (x != nullptr) {
            out.push_back(static_cast<long long>(x->forward.size()));
            x = x->forward[0];
        }
        return out;
    }

    bool validate() const {
        for (int i = 0; i < level_; ++i) {
            Node* x = head_;
            while (x->forward[i] != nullptr) {
                Node* target = x->forward[i];
                long long d = 0;
                Node* y = x;
                while (y != target) {
                    y = y->forward[0];
                    if (y == nullptr) return false;
                    ++d;
                    if (d > length_ + 1) return false;
                }
                if (d != x->span[i]) return false;
                x = target;
            }
        }
        return true;
    }

private:
    unsigned long long rand64() {
        state_ = state_ * LCG_MUL + LCG_ADD;   // uint64 自然溢位，等同 mod 2^64
        return state_;
    }

    int randomLevel() {
        int lv = 1;
        while (lv < MAX_LEVEL) {
            if ((rand64() >> 33) & 1ULL) break;
            ++lv;
        }
        return lv;
    }

    void findPreds(long long key, bool strict, std::vector<Node*>& update,
                   std::vector<long long>& rank) {
        Node* x = head_;
        for (int i = level_ - 1; i >= 0; --i) {
            if (i + 1 < level_) rank[i] = rank[i + 1];
            while (true) {
                Node* nxt = x->forward[i];
                if (nxt == nullptr) break;
                bool go = strict ? (nxt->key < key) : (nxt->key <= key);
                if (!go) break;
                rank[i] += x->span[i];
                x = nxt;
            }
            update[i] = x;
        }
    }

    unsigned long long state_;
    Node* head_;
    int level_;
    long long length_;
};

// ------------------------------------------------- 參考實作（有序陣列 + 二分）

class RefSet {
public:
    void insert(long long x) {
        auto it = std::lower_bound(a_.begin(), a_.end(), x);
        a_.insert(it, x);
    }

    bool erase(long long x) {
        auto it = std::lower_bound(a_.begin(), a_.end(), x);
        if (it == a_.end() || *it != x) return false;
        a_.erase(it);
        return true;
    }

    long long count(long long x) const {
        return static_cast<long long>(std::upper_bound(a_.begin(), a_.end(), x) -
                                      std::lower_bound(a_.begin(), a_.end(), x));
    }

    long long rank(long long x) const {
        return static_cast<long long>(std::lower_bound(a_.begin(), a_.end(), x) - a_.begin());
    }

    long long kth(long long k) const {
        if (k < 1 || k > static_cast<long long>(a_.size())) return -1;
        return a_[static_cast<size_t>(k - 1)];
    }

    long long predecessor(long long x) const {
        long long i = rank(x);
        return i > 0 ? a_[static_cast<size_t>(i - 1)] : -1;
    }

    long long successor(long long x) const {
        auto it = std::upper_bound(a_.begin(), a_.end(), x);
        return it == a_.end() ? -1 : *it;
    }

    long long size() const { return static_cast<long long>(a_.size()); }
    const std::vector<long long>& data() const { return a_; }

private:
    std::vector<long long> a_;
};

// ---------------------------------------------------------------- 主流程

std::string run(const std::vector<std::pair<int, long long>>& ops) {
    SkipList sl;
    RefSet ref;
    std::vector<long long> res;
    bool ok = true;

    for (const auto& op : ops) {
        int code = op.first;
        long long x = op.second;
        long long sv = 0, rv = 0;
        if (code == 1) {
            sl.insert(x);
            ref.insert(x);
        } else if (code == 2) {
            bool a = sl.erase(x);
            bool b = ref.erase(x);
            if (a != b) ok = false;
        } else if (code == 3) {
            sv = sl.count(x);
            rv = ref.count(x);
            res.push_back(sv);
            if (sv != rv) ok = false;
        } else if (code == 4) {
            sv = sl.countLess(x);
            rv = ref.rank(x);
            res.push_back(sv);
            if (sv != rv) ok = false;
        } else if (code == 5) {
            sv = sl.kth(x);
            rv = ref.kth(x);
            res.push_back(sv);
            if (sv != rv) ok = false;
        } else if (code == 6) {
            sv = sl.predecessor(x);
            rv = ref.predecessor(x);
            res.push_back(sv);
            if (sv != rv) ok = false;
        } else if (code == 7) {
            sv = sl.successor(x);
            rv = ref.successor(x);
            res.push_back(sv);
            if (sv != rv) ok = false;
        }
    }

    if (sl.size() != ref.size()) ok = false;
    if (sl.toList() != ref.data()) ok = false;
    if (!sl.validate()) ok = false;

    std::ostringstream oss;
    for (size_t i = 0; i < res.size(); ++i) {
        if (i) oss << ' ';
        oss << res[i];
    }
    oss << '\n';
    oss << sl.size();
    for (long long v : sl.toList()) oss << ' ' << v;
    oss << '\n';
    oss << (ok ? '1' : '0');
    oss << '\n';
    return oss.str();
}

// ---------------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例 1：基本插入 ----
    {
        std::vector<std::pair<int, long long>> ops = {
            {1, 5}, {1, 3}, {1, 7}, {3, 3}, {4, 5}, {5, 2}};
        std::string out = run(ops);
        if (out != "1 1 5\n3 3 5 7\n1\n") {
            std::cout << "FAIL case1:\n" << out << std::endl;
            std::exit(1);
        }
    }
    // ---- 固定用例 2：重複元素 ----
    {
        std::vector<std::pair<int, long long>> ops = {
            {1, 2}, {1, 2}, {1, 2}, {3, 2}, {4, 2}, {4, 3}};
        std::string out = run(ops);
        if (out != "3 0 3\n3 2 2 2\n1\n") {
            std::cout << "FAIL case2:\n" << out << std::endl;
            std::exit(1);
        }
    }
    // ---- 固定用例 3：前驅後繼 ----
    {
        std::vector<std::pair<int, long long>> ops = {
            {1, 10}, {1, 20}, {1, 30}, {6, 20}, {7, 20}, {6, 10}, {7, 30}};
        std::string out = run(ops);
        if (out != "10 30 -1 -1\n3 10 20 30\n1\n") {
            std::cout << "FAIL case3:\n" << out << std::endl;
            std::exit(1);
        }
    }
    // ---- 固定用例 4：刪除不存在 / 刪空 ----
    {
        std::vector<std::pair<int, long long>> ops = {
            {1, 1}, {2, 9}, {2, 1}, {3, 1}, {5, 1}};
        std::string out = run(ops);
        if (out != "0 -1\n0\n1\n") {
            std::cout << "FAIL case4:\n" << out << std::endl;
            std::exit(1);
        }
    }
    // ---- 空操作 ----
    {
        std::string out = run({});
        if (out != "\n0\n1\n") {
            std::cout << "FAIL empty:\n" << out << std::endl;
            std::exit(1);
        }
    }
    // ---- 層高序列跨語言一致性（回歸用） ----
    {
        SkipList s2;
        for (int i = 0; i < 8; ++i) s2.insert(static_cast<long long>(i));
        std::vector<long long> got = s2.levels();
        std::vector<long long> want = {3, 1, 2, 4, 1, 3, 1, 1};
        if (got != want) {
            std::cout << "FAIL levels:";
            for (long long v : got) std::cout << ' ' << v;
            std::cout << std::endl;
            std::exit(1);
        }
    }
    // ---- 隨機對拍（與 Python 用同一組種子的自製 LCG） ----
    {
        unsigned long long s = 20261019ULL;
        auto nxt = [&]() -> unsigned long long {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            return (s >> 32) & 0xFFFFFFFFULL;
        };
        for (int t = 0; t < 700; ++t) {
            long long q = static_cast<long long>(nxt() % 41ULL);
            std::vector<std::pair<int, long long>> ops;
            for (long long i = 0; i < q; ++i) {
                long long code = static_cast<long long>(nxt() % 7ULL) + 1;
                long long x = static_cast<long long>(nxt() % 13ULL);
                ops.push_back({static_cast<int>(code), x});
            }
            std::string out = run(ops);
            std::string last;
            size_t p = out.rfind('\n', out.size() - 2);
            if (p == std::string::npos) last = out;
            else last = out.substr(p + 1, out.size() - 2 - p);
            if (last != "1") {
                std::cout << "FAIL random t=" << t << "\n" << out << std::endl;
                std::exit(1);
            }
        }
    }
    // ---- 大量插入後隨機查詢 ----
    {
        unsigned long long s = 4242ULL;
        auto nxt = [&]() -> unsigned long long {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            return (s >> 32) & 0xFFFFFFFFULL;
        };
        for (int t = 0; t < 60; ++t) {
            std::vector<std::pair<int, long long>> ops;
            long long cnt = static_cast<long long>(nxt() % 121ULL);
            for (long long i = 0; i < cnt; ++i) {
                ops.push_back({1, static_cast<long long>(nxt() % 201ULL)});
            }
            for (int i = 0; i < 30; ++i) {
                ops.push_back({static_cast<int>(nxt() % 6ULL) + 2,
                               static_cast<long long>(nxt() % 211ULL)});
            }
            std::string out = run(ops);
            size_t p = out.rfind('\n', out.size() - 2);
            std::string last = (p == std::string::npos) ? out
                                                        : out.substr(p + 1, out.size() - 2 - p);
            if (last != "1") {
                std::cout << "FAIL bulk t=" << t << "\n" << out << std::endl;
                std::exit(1);
            }
        }
    }
}

// ---------------------------------------------------------------- IO

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
        return clampLL(s, -1000000000000000000LL, 1000000000000000000LL, def);
    };

    long long q = nextInt(0);
    if (q < 0) q = 0;
    if (q > QMAX) q = QMAX;
    std::vector<std::pair<int, long long>> ops;
    for (long long i = 0; i < q; ++i) {
        long long op = nextInt(0);
        long long x = nextInt(0);
        if (op < 1 || op > 7) op = 0;
        if (x < 0) x = 0;
        if (x > VMAX) x = VMAX;
        ops.push_back({static_cast<int>(op), x});
    }
    std::cout << run(ops);
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
