// 替罪羊樹（Scapegoat Tree）：不旋轉，靠「整棵子樹重建」維持平衡的權重平衡樹
//
// 題意：
//     維護一個**有序集合**（元素唯一），支援插入 / 刪除 / 查找 / 排名 / 第 k 小 /
//     前驅 / 後繼，並與「有序陣列 + 二分」的參考實作逐項對拍。
//
//     操作碼（每行 `op x`）：
//         1  insert x    插入 x（已存在則忽略）
//         2  erase  x    刪除 x（不存在則忽略）
//         3  count  x    x 是否在集合中（1 / 0）
//         4  rank   x    嚴格小於 x 的元素個數
//         5  kth    k    第 k 小的元素（k 從 1 起算），越界回傳 -1
//         6  pred   x    嚴格小於 x 的最大元素，不存在回傳 -1
//         7  succ   x    嚴格大於 x 的最小元素，不存在回傳 -1
//
// 思路：
//     ### 1. 什麼是「權重平衡」
//     紅黑樹 / AVL / Treap / Splay 都靠**旋轉**（rotation）局部調整形狀。
//     替罪羊樹完全不旋轉，改成一個更粗暴的條件：
//
//         對每個節點 u，其左右子樹大小都必須 ≤ α · size(u)   （本題 α = 7/10）
//
//     一旦某個節點違反這個條件，就把它的**整棵子樹**拉直（中序遍歷成有序陣列），
//     再從中點遞迴建成一棵完美平衡樹。這個動作叫 **rebuild（重建）**。
//
//     ### 2. 為什麼這樣還能是 O(log n)
//     - α 權重平衡 ⟹ 從根往下每走一層，子樹大小至少縮為 α 倍
//       ⟹ 樹高 ≤ log_{1/α}(n) = ln n / ln(10/7) ≈ 2.81 · log₂ n。
//     - 重建的攤還代價：一個節點被重建，代表它底下剛累積了 Θ(size) 次插入，
//       於是把 O(size) 的重建成本分攤到這 Θ(size) 次插入上，每次插入攤還 O(log n)。
//
//     ### 3. 插入後要重建「哪一個」節點（本題的關鍵設計）
//     插入後從根沿著 key 走一趟，找**最淺**（離根最近）的那個失衡節點 v 重建。
//     重建完之後，subtree(v) 內部全部完美平衡，而 v 的祖先們大小完全沒變
//     （重建不改變節點數），所以「下一次失衡的位置」只能嚴格地往根的方向移動。
//     於是把這個動作寫成迴圈：
//
//         while 存在失衡節點 v（取最淺的）:  rebuild(v)
//
//     每一輪 v 的深度嚴格遞減，最多走 depth 輪就結束；一旦 v 就是根，
//     重建後整棵樹完美平衡，必然收斂。這比「找最深的 / 只重建一次」更穩健——
//     因為插入會讓祖先的 size 與子樹 size 同時 +1，而 α < 1，
//     `child+1 ≤ α·(size+1)` 是有可能被打破的（1 > α），只修一次不保證全域平衡。
//
//     ### 4. 刪除：靠「縮水太多就整棵重建」
//     刪除走一般 BST 刪除（兩個小孩時用**中序前驅**頂替），沿途把 size -1。
//     維護 max_size（歷史最大節點數）：當 `size · 10 < max_size · 7` 時，
//     表示節點已經掉到歷史高峰的 70% 以下，整棵樹重建並把 max_size 設回 size。
//     這樣攤還起來每次刪除也是 O(log n)。
//
//     ### 5. 兩語言一致性的兩個陷阱
//     - 平衡條件用**整數**比較 `max(l, r) · 10 > size · 7`，不要用浮點 α。
//     - 建平衡樹取中點用 `mid = (l + r) / 2`（下標非負，C++ 的 truncate 與
//       Python 的 floor 結果相同），不要寫 `(l + r) / 2` 配上負座標。
//
// 應用場景：
//     - 需要「有序集合」但又想避開旋轉邏輯的場合（程式碼量比紅黑樹少一個數量級）
//     - 函數式 / 持久化結構：重建天然產生新節點，不需要可變指標
//     - 資料庫索引中「批次重建比逐次旋轉更快」的場景（LSM 的 read-only segment 合併）
//     - 競賽中 Treap / Splay 的替代：特別是需要「按子樹大小攤平」的題
//     - 不希望旋轉造成迭代器 / 外部指標失效的嵌入式與即時系統
//
// 複雜度（n = 元素個數，α = 0.7）：
//     查找 / rank / kth / pred / succ   O(log n) 最壞
//     插入                             攤還 O(log n)，單次最壞 O(n)（重建）
//     刪除                             攤還 O(log n)，單次最壞 O(n)
//     樹高                             ≤ log_{1/α}(n) + O(1) ≈ 2.81 log₂ n
//     空間                             O(n)

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

const int QMAX = 300;               // IO 模式下的操作數上限
const long long VMAX = 1000000;     // 元素值上限（夾到 [0, VMAX]，讓 -1 可以當「不存在」哨兵）

const long long ALPHA_NUM = 7;      // α = 7/10
const long long ALPHA_DEN = 10;

// ---------------------------------------------------------------- 節點與小工具

struct Node {
    long long key;
    Node* left;
    Node* right;
    long long size;
    explicit Node(long long k) : key(k), left(nullptr), right(nullptr), size(1) {}
};

inline long long sz(Node* t) { return t == nullptr ? 0 : t->size; }

inline void pull(Node* t) { t->size = 1 + sz(t->left) + sz(t->right); }

void collect(Node* t, std::vector<long long>& out) {
    // 中序遍歷，把子樹拉直成有序陣列
    if (t == nullptr) return;
    collect(t->left, out);
    out.push_back(t->key);
    collect(t->right, out);
}

Node* buildTree(const std::vector<long long>& keys, long long l, long long r) {
    // 把有序陣列 keys[l:r] 建成一棵高度平衡的樹
    if (l >= r) return nullptr;
    long long mid = (l + r) / 2;
    Node* t = new Node(keys[static_cast<size_t>(mid)]);
    t->left = buildTree(keys, l, mid);
    t->right = buildTree(keys, mid + 1, r);
    pull(t);
    return t;
}

int heightOf(Node* t) {
    if (t == nullptr) return 0;
    return 1 + std::max(heightOf(t->left), heightOf(t->right));
}

long long trueSize(Node* t) {
    // 實際走一趟算出來的子樹大小，用來驗證 size 欄位沒有算錯
    if (t == nullptr) return 0;
    return 1 + trueSize(t->left) + trueSize(t->right);
}

bool checkSizes(Node* t) {
    if (t == nullptr) return true;
    if (t->size != 1 + sz(t->left) + sz(t->right)) return false;
    return checkSizes(t->left) && checkSizes(t->right);
}

inline bool unbalanced(Node* t) {
    // α 權重平衡條件（α = 7/10）。整數運算，兩語言完全一致。
    return std::max(sz(t->left), sz(t->right)) * ALPHA_DEN > t->size * ALPHA_NUM;
}

// ---------------------------------------------------------------- 替罪羊樹

class ScapegoatTree {
public:
    ScapegoatTree() : root_(nullptr), maxSize_(0) {}

    long long size() const { return sz(root_); }
    int height() const { return heightOf(root_); }

    std::vector<long long> inorder() const {
        std::vector<long long> out;
        collect(root_, out);
        return out;
    }

    // 給測試用的自我檢查：size 欄位必須等於真實子樹大小
    bool sizeInvariant() const { return checkSizes(root_); }
    long long trueSizeAll() const { return trueSize(root_); }

    bool contains(long long key) const {
        Node* t = root_;
        while (t != nullptr) {
            if (key == t->key) return true;
            t = (key < t->key) ? t->left : t->right;
        }
        return false;
    }

    long long countLess(long long key) const {
        long long cnt = 0;
        Node* t = root_;
        while (t != nullptr) {
            if (key <= t->key) {
                t = t->left;
            } else {
                cnt += 1 + sz(t->left);
                t = t->right;
            }
        }
        return cnt;
    }

    long long kth(long long k) const {
        // 第 k 小（k 從 1 起算）
        if (k < 1 || k > sz(root_)) return -1;
        Node* t = root_;
        while (t != nullptr) {
            long long ls = sz(t->left);
            if (k <= ls) {
                t = t->left;
            } else if (k == ls + 1) {
                return t->key;
            } else {
                k -= ls + 1;
                t = t->right;
            }
        }
        return -1;
    }

    long long pred(long long key) const {
        long long best = -1;
        Node* t = root_;
        while (t != nullptr) {
            if (t->key < key) {
                best = t->key;
                t = t->right;
            } else {
                t = t->left;
            }
        }
        return best;
    }

    long long succ(long long key) const {
        long long best = -1;
        Node* t = root_;
        while (t != nullptr) {
            if (t->key > key) {
                best = t->key;
                t = t->left;
            } else {
                t = t->right;
            }
        }
        return best;
    }

    bool insert(long long key) {
        if (root_ == nullptr) {
            root_ = new Node(key);
            maxSize_ = 1;
            return true;
        }
        std::vector<Node*> path;
        Node* t = root_;
        while (t != nullptr) {
            if (key == t->key) return false;
            path.push_back(t);
            t = (key < t->key) ? t->left : t->right;
        }
        Node* parent = path.back();
        if (key < parent->key) {
            parent->left = new Node(key);
        } else {
            parent->right = new Node(key);
        }
        for (Node* u : path) u->size += 1;
        if (root_->size > maxSize_) maxSize_ = root_->size;
        fix(key);
        return true;
    }

    bool erase(long long key) {
        std::vector<Node*> path;
        Node* t = root_;
        while (t != nullptr && t->key != key) {
            path.push_back(t);
            t = (key < t->key) ? t->left : t->right;
        }
        if (t == nullptr) return false;
        if (t->left != nullptr && t->right != nullptr) {
            // 兩個小孩：用中序前驅（左子樹最大值）頂替
            std::vector<Node*> ppath;
            ppath.push_back(t);
            Node* p = t->left;
            while (p->right != nullptr) {
                ppath.push_back(p);
                p = p->right;
            }
            t->key = p->key;
            Node* parent = ppath.back();
            if (parent->left == p) {
                parent->left = p->left;
            } else {
                parent->right = p->left;
            }
            // 注意：t 的祖先（path）也要 -1，ppath 只包含 t 往下到 p 的父節點
            for (Node* u : path) u->size -= 1;
            for (Node* u : ppath) u->size -= 1;
        } else {
            Node* child = (t->left != nullptr) ? t->left : t->right;
            if (!path.empty()) {
                Node* parent = path.back();
                if (parent->left == t) {
                    parent->left = child;
                } else {
                    parent->right = child;
                }
            } else {
                root_ = child;
            }
            for (Node* u : path) u->size -= 1;
        }
        if (root_ != nullptr && root_->size * ALPHA_DEN < maxSize_ * ALPHA_NUM) {
            std::vector<long long> keys;
            collect(root_, keys);
            root_ = buildTree(keys, 0, static_cast<long long>(keys.size()));
            maxSize_ = root_->size;
        } else if (root_ == nullptr) {
            maxSize_ = 0;
        }
        return true;
    }

private:
    // 插入後反覆重建最淺的失衡節點
    void fix(long long key) {
        while (true) {
            std::vector<Node*> path;
            Node* t = root_;
            while (t != nullptr) {
                path.push_back(t);
                if (key == t->key) break;
                t = (key < t->key) ? t->left : t->right;
            }
            long long bad = -1;
            for (size_t i = 0; i < path.size(); ++i) {
                if (unbalanced(path[i])) {
                    bad = static_cast<long long>(i);
                    break;
                }
            }
            if (bad < 0) return;
            Node* v = path[static_cast<size_t>(bad)];
            std::vector<long long> keys;
            collect(v, keys);
            Node* nb = buildTree(keys, 0, static_cast<long long>(keys.size()));
            if (bad == 0) {
                root_ = nb;
            } else {
                Node* p = path[static_cast<size_t>(bad - 1)];
                if (p->left == v) {
                    p->left = nb;
                } else {
                    p->right = nb;
                }
            }
        }
    }

    Node* root_;
    long long maxSize_;
};

// ---------------------------------------------------------------- 參考實作對拍

std::string run(const std::vector<std::pair<int, long long>>& ops) {
    ScapegoatTree tree;
    std::vector<long long> ref;
    std::vector<std::string> qres;

    for (const auto& pr : ops) {
        int op = pr.first;
        long long x = pr.second;
        if (op == 1) {
            bool a = tree.insert(x);
            bool b = false;
            if (std::find(ref.begin(), ref.end(), x) == ref.end()) {
                ref.insert(std::lower_bound(ref.begin(), ref.end(), x), x);
                b = true;
            }
            if (a != b) {
                std::cout << "FATAL insert " << x << "\n";
                std::abort();
            }
        } else if (op == 2) {
            bool a = tree.erase(x);
            bool b = false;
            auto it = std::lower_bound(ref.begin(), ref.end(), x);
            if (it != ref.end() && *it == x) {
                ref.erase(it);
                b = true;
            }
            if (a != b) {
                std::cout << "FATAL erase " << x << "\n";
                std::abort();
            }
        } else if (op == 3) {
            qres.push_back(tree.contains(x) ? "1" : "0");
            auto it = std::lower_bound(ref.begin(), ref.end(), x);
            bool want = (it != ref.end() && *it == x);
            if (tree.contains(x) != want) {
                std::cout << "FATAL contains " << x << "\n";
                std::abort();
            }
        } else if (op == 4) {
            qres.push_back(std::to_string(tree.countLess(x)));
            auto it = std::lower_bound(ref.begin(), ref.end(), x);
            long long want = static_cast<long long>(it - ref.begin());
            if (tree.countLess(x) != want) {
                std::cout << "FATAL rank " << x << "\n";
                std::abort();
            }
        } else if (op == 5) {
            qres.push_back(std::to_string(tree.kth(x)));
            long long want = (x >= 1 && x <= static_cast<long long>(ref.size()))
                                 ? ref[static_cast<size_t>(x - 1)]
                                 : -1;
            if (tree.kth(x) != want) {
                std::cout << "FATAL kth " << x << " got " << tree.kth(x) << " want " << want << "\n";
                std::abort();
            }
        } else if (op == 6) {
            qres.push_back(std::to_string(tree.pred(x)));
            auto it = std::lower_bound(ref.begin(), ref.end(), x);
            long long want = (it != ref.begin()) ? *(--it) : -1;
            if (tree.pred(x) != want) {
                std::cout << "FATAL pred " << x << "\n";
                std::abort();
            }
        } else if (op == 7) {
            qres.push_back(std::to_string(tree.succ(x)));
            auto it = std::upper_bound(ref.begin(), ref.end(), x);   // 嚴格大於：用右插入點
            long long want = (it != ref.end()) ? *it : -1;
            if (tree.succ(x) != want) {
                std::cout << "FATAL succ " << x << "\n";
                std::abort();
            }
        }
    }

    // size 欄位的一致性檢查：中序內容對但 size 錯，kth / rank 仍會給錯答案
    if (!tree.sizeInvariant() || tree.trueSizeAll() != tree.size()) {
        std::cout << "FATAL size invariant\n";
        std::abort();
    }

    std::ostringstream oss;
    for (size_t i = 0; i < qres.size(); ++i) {
        if (i) oss << " ";
        oss << qres[i];
    }
    oss << "\n";

    std::vector<long long> io = tree.inorder();
    oss << tree.size() << " " << tree.height();
    for (long long v : io) oss << " " << v;
    oss << "\n";

    int ok = (io == ref) ? 1 : 0;
    oss << ok << "\n";
    return oss.str();
}

// ---------------------------------------------------------------- 自製 LCG

unsigned long long lcgState = 0x9E3779B97F4A7C15ULL;

unsigned long long lcgNext() {
    lcgState = lcgState * 6364136223846793005ULL + 1442695040888963407ULL;
    return lcgState;
}

long long lcgInt(long long lo, long long hi) {
    unsigned long long r = lcgNext() >> 11;          // 53 位元，避免乘法溢位疑慮
    long long span = hi - lo + 1;
    return lo + static_cast<long long>(r % static_cast<unsigned long long>(span));
}

// 取出輸出第三行（ok 旗標）
int flagOf(const std::string& out) {
    size_t p = out.rfind('\n', out.size() - 2);
    std::string s = out.substr(p + 1);
    return std::atoi(s.c_str());
}

int bitLen(long long v) {
    int b = 0;
    while (v > 0) {
        ++b;
        v >>= 1;
    }
    return b;
}

// ---------------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例：基本插入與查詢 ----
    {
        std::vector<std::pair<int, long long>> ops = {
            {1, 5}, {1, 3}, {1, 8}, {1, 1}, {1, 5}, {3, 3}, {3, 4},
            {4, 5}, {5, 2}, {6, 5}, {7, 5}};
        std::string out = run(ops);
        if (out != "1 0 2 3 3 8\n4 3 1 3 5 8\n1\n") {
            std::cout << "FATAL fixed1: [" << out << "]\n";
            std::abort();
        }
    }

    // ---- 空集合 ----
    {
        std::vector<std::pair<int, long long>> ops = {
            {5, 1}, {6, 0}, {7, 0}, {3, 0}, {4, 0}, {2, 0}};
        std::string out = run(ops);
        if (out != "-1 -1 -1 0 0\n0 0\n1\n") {
            std::cout << "FATAL fixed2: [" << out << "]\n";
            std::abort();
        }
    }

    // ---- 樹高：連續插入 1..2000（已排序輸入，對 BST 是最壞情況）----
    {
        ScapegoatTree tree;
        for (long long i = 1; i <= 2000; ++i) tree.insert(i);
        if (tree.size() != 2000) {
            std::cout << "FATAL sorted size\n";
            std::abort();
        }
        if (tree.height() > 30) {
            std::cout << "FATAL sorted height " << tree.height() << "\n";
            std::abort();
        }
        std::vector<long long> io = tree.inorder();
        for (long long i = 0; i < 2000; ++i) {
            if (io[static_cast<size_t>(i)] != i + 1) {
                std::cout << "FATAL sorted order\n";
                std::abort();
            }
        }
    }

    // ---- 樹高：隨機插入 2000 個 ----
    {
        lcgState = 1234ULL;
        ScapegoatTree tree;
        std::vector<long long> vals;
        for (int i = 0; i < 2000; ++i) vals.push_back(lcgInt(0, 1000000));
        for (long long v : vals) tree.insert(v);
        std::vector<long long> want = vals;
        std::sort(want.begin(), want.end());
        want.erase(std::unique(want.begin(), want.end()), want.end());
        if (tree.inorder() != want) {
            std::cout << "FATAL random order\n";
            std::abort();
        }
        if (tree.height() > 30) {
            std::cout << "FATAL random height " << tree.height() << "\n";
            std::abort();
        }
    }

    // ---- 隨機對拍 ----
    {
        lcgState = 20261020ULL;
        for (int t = 0; t < 700; ++t) {
            int q = static_cast<int>(lcgInt(0, 40));
            std::vector<std::pair<int, long long>> ops;
            for (int i = 0; i < q; ++i) {
                ops.push_back({static_cast<int>(lcgInt(1, 7)), lcgInt(0, 12)});
            }
            std::string out = run(ops);
            if (flagOf(out) != 1) {
                std::cout << "FATAL fuzz " << t << ": [" << out << "]\n";
                std::abort();
            }
        }
    }

    // ---- 大量插入後隨機查詢（順便檢查樹高）----
    {
        for (int trial = 0; trial < 60; ++trial) {
            lcgState = 9000ULL + static_cast<unsigned long long>(trial);
            std::vector<std::pair<int, long long>> ops;
            int m = static_cast<int>(lcgInt(0, 150));
            for (int i = 0; i < m; ++i) ops.push_back({1, lcgInt(0, 200)});
            for (int i = 0; i < 30; ++i) ops.push_back({static_cast<int>(lcgInt(2, 7)), lcgInt(0, 210)});
            std::string out = run(ops);
            if (flagOf(out) != 1) {
                std::cout << "FATAL big " << trial << ": [" << out << "]\n";
                std::abort();
            }
            // 解析第二行的 size 與 height
            size_t nl1 = out.find('\n');
            size_t nl2 = out.find('\n', nl1 + 1);
            std::string line2 = out.substr(nl1 + 1, nl2 - nl1 - 1);
            std::istringstream iss(line2);
            long long n = 0, h = 0;
            iss >> n >> h;
            if (n >= 2 && h > 3 * bitLen(n) + 6) {
                std::cout << "FATAL height " << trial << " n=" << n << " h=" << h << "\n";
                std::abort();
            }
        }
    }

    // ---- 刪除壓力測試：插入 800 個再全部刪掉 ----
    {
        lcgState = 555ULL;
        std::vector<std::pair<int, long long>> ops;
        std::vector<long long> vals;
        for (int i = 0; i < 800; ++i) vals.push_back(lcgInt(0, 500));
        for (long long v : vals) ops.push_back({1, v});
        for (long long v : vals) ops.push_back({2, v});
        std::string out = run(ops);
        size_t nl1 = out.find('\n');
        size_t nl2 = out.find('\n', nl1 + 1);
        std::string line2 = out.substr(nl1 + 1, nl2 - nl1 - 1);
        if (line2 != "0 0") {
            std::cout << "FATAL delall: [" << line2 << "]\n";
            std::abort();
        }
    }
}

// ---------------------------------------------------------------- IO

bool isInt(const std::string& s) {
    if (s.empty()) return false;
    size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
    if (i >= s.size()) return false;
    for (size_t j = i; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') return false;
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
