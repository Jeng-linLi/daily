// 二元搜尋樹（BST）：插入 / 搜尋 / 刪除 / 中序遍歷 / 樹高 / 合法性檢查
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：BST 的核心不變量是「左子樹所有鍵值 < x.val <= 右子樹所有鍵值」，
//   有了它，搜尋就退化成二分，插入 / 搜尋 / 刪除都只沿一條根到葉的路徑走，複雜度 O(h)。
//   不做平衡的樸素 BST 在有序插入時會退化成鏈（h = n），這是樹類題目最常見的性能坑。
//
//   刪除分三種情形：葉子直接摘；只有一棵子樹就讓子樹頂替；
//   左右子樹都在時找**中序後繼**（右子樹最小節點）把值搬上來，再刪掉那個後繼。
//   後繼至多只有右孩子，所以第二趟刪除必定落在前兩種簡單情形。
//   這裏用一個棧上的哨兵節點當「根的父節點」，連刪根都不用特判。
//   遍歷 / 求高 / 計數同樣寫成迭代版，行爲與 Python 版完全一致。
//
// 輸入（空白分隔）：n m / a1..an（插入）/ b1..bm（刪除）
// 輸出：插入後中序、插入後樹高、m 個搜尋結果（1/0）、刪除後中序、刪除後樹高
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using ll = long long;

struct Node {
    ll val;
    Node* left;
    Node* right;
    explicit Node(ll v) : val(v), left(nullptr), right(nullptr) {}
};

// 插入 val（重複值走右子樹）。時間 O(h)，空間 O(1)
Node* insert(Node* root, ll val) {
    if (root == nullptr) return new Node(val);
    Node* cur = root;
    while (true) {
        if (val < cur->val) {                 // 嚴格小於走左邊
            if (cur->left == nullptr) {
                cur->left = new Node(val);
                return root;
            }
            cur = cur->left;
        } else {                              // 大於或等於都走右邊（重複鍵值策略）
            if (cur->right == nullptr) {
                cur->right = new Node(val);
                return root;
            }
            cur = cur->right;
        }
    }
}

// 判斷 val 是否存在。時間 O(h)，空間 O(1)
bool search(Node* root, ll val) {
    Node* cur = root;
    while (cur != nullptr) {
        if (val == cur->val) return true;
        cur = (val < cur->val) ? cur->left : cur->right;
    }
    return false;
}

// 把 parent 的某個孩子 old 換成 new（用指針身份比較，不用值比較）
static void replaceChild(Node* parent, Node* oldNode, Node* newNode) {
    if (parent->left == oldNode) parent->left = newNode;
    else parent->right = newNode;
}

// 刪除一個值爲 val 的節點（不存在則原樣返回）。時間 O(h)，空間 O(1)
Node* deleteNode(Node* root, ll val) {
    Node dummy(0);
    dummy.right = root;                        // 哨兵：當「根的父節點」
    Node* parent = &dummy;
    Node* cur = root;
    while (cur != nullptr && cur->val != val) {
        parent = cur;
        cur = (val < cur->val) ? cur->left : cur->right;
    }
    if (cur == nullptr) return dummy.right;    // 沒找到

    if (cur->left == nullptr || cur->right == nullptr) {   // 情形 1 / 2
        Node* child = (cur->left != nullptr) ? cur->left : cur->right;
        replaceChild(parent, cur, child);
        delete cur;
    } else {                                               // 情形 3：找中序後繼
        Node* succParent = cur;
        Node* succ = cur->right;
        while (succ->left != nullptr) {
            succParent = succ;
            succ = succ->left;
        }
        cur->val = succ->val;                  // 值搬上來，結構不動
        replaceChild(succParent, succ, succ->right);
        delete succ;
    }
    return dummy.right;
}

// 中序遍歷（迭代版，顯式棧）。時間 O(n)，空間 O(h)
vector<ll> inorder(Node* root) {
    vector<ll> out;
    vector<Node*> st;
    Node* cur = root;
    while (cur != nullptr || !st.empty()) {
        while (cur != nullptr) {               // 一路向左到底
            st.push_back(cur);
            cur = cur->left;
        }
        cur = st.back();
        st.pop_back();
        out.push_back(cur->val);
        cur = cur->right;
    }
    return out;
}

// 樹高：空樹 0，單節點 1。迭代 DFS
int height(Node* root) {
    if (root == nullptr) return 0;
    int best = 0;
    vector<pair<Node*, int>> st;
    st.push_back({root, 1});
    while (!st.empty()) {
        auto [node, depth] = st.back();
        st.pop_back();
        if (depth > best) best = depth;
        if (node->left != nullptr) st.push_back({node->left, depth + 1});
        if (node->right != nullptr) st.push_back({node->right, depth + 1});
    }
    return best;
}

// 節點個數（迭代版）
int treeSize(Node* root) {
    int cnt = 0;
    vector<Node*> st;
    if (root != nullptr) st.push_back(root);
    while (!st.empty()) {
        Node* node = st.back();
        st.pop_back();
        ++cnt;
        if (node->left != nullptr) st.push_back(node->left);
        if (node->right != nullptr) st.push_back(node->right);
    }
    return cnt;
}

// 檢查 BST 不變量：左子樹嚴格小於，右子樹允許等於。區間邊界是否可取用旗標區分
bool isValid(Node* root) {
    if (root == nullptr) return true;
    const ll NEG = numeric_limits<ll>::min();
    const ll POS = numeric_limits<ll>::max();
    struct Item {
        Node* node;
        ll lo, hi;
        bool loStrict, hiStrict;
    };
    vector<Item> st;
    st.push_back({root, NEG, POS, false, false});
    while (!st.empty()) {
        Item it = st.back();
        st.pop_back();
        ll v = it.node->val;
        if (v < it.lo || (it.loStrict && v == it.lo)) return false;
        if (v > it.hi || (it.hiStrict && v == it.hi)) return false;
        if (it.node->left != nullptr)          // 左子樹：上界爲當前值，必須嚴格小於
            st.push_back({it.node->left, it.lo, v, it.loStrict, true});
        if (it.node->right != nullptr)         // 右子樹：下界爲當前值，可取等號
            st.push_back({it.node->right, v, it.hi, false, it.hiStrict});
    }
    return true;
}

// 依序插入整個序列
Node* build(const vector<ll>& values) {
    Node* root = nullptr;
    for (ll v : values) root = insert(root, v);
    return root;
}

static void freeTree(Node* root) {
    if (root == nullptr) return;
    freeTree(root->left);
    freeTree(root->right);
    delete root;
}

// ---------------- 獨立對照實現（用於對拍） ----------------

// 另一種求高寫法：層序 BFS 數層數
int heightBfs(Node* root) {
    if (root == nullptr) return 0;
    vector<Node*> level{root};
    int h = 0;
    while (!level.empty()) {
        ++h;
        vector<Node*> nxt;
        for (Node* node : level) {
            if (node->left != nullptr) nxt.push_back(node->left);
            if (node->right != nullptr) nxt.push_back(node->right);
        }
        level.swap(nxt);
    }
    return h;
}

// 遞迴版中序遍歷，用於交叉驗證迭代版
void inorderRec(Node* root, vector<ll>& out) {
    if (root == nullptr) return;
    inorderRec(root->left, out);
    out.push_back(root->val);
    inorderRec(root->right, out);
}

// ---------------- IO ----------------

template <typename T>
static string joinInts(const vector<T>& v) {
    ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}

static void runIo(const string& data) {
    istringstream iss(data);
    ll n = 0, m = 0;
    if (!(iss >> n)) n = 0;
    if (!(iss >> m)) m = 0;
    vector<ll> inserts, deletes;
    for (ll i = 0; i < n; ++i) {
        ll x = 0;
        if (!(iss >> x)) x = 0;
        inserts.push_back(x);
    }
    for (ll i = 0; i < m; ++i) {
        ll x = 0;
        if (!(iss >> x)) x = 0;
        deletes.push_back(x);
    }

    Node* root = build(inserts);
    cout << joinInts(inorder(root)) << '\n';      // 插入後的中序
    cout << height(root) << '\n';                 // 插入後的樹高

    ostringstream oss;
    for (size_t i = 0; i < deletes.size(); ++i) {
        if (i) oss << ' ';
        oss << (search(root, deletes[i]) ? 1 : 0);
    }
    cout << oss.str() << '\n';                    // 逐個搜尋

    for (ll v : deletes) root = deleteNode(root, v);
    cout << joinInts(inorder(root)) << '\n';      // 刪除後的中序
    cout << height(root) << '\n';                 // 刪除後的樹高
    freeTree(root);
}

static void runTests() {
    // README 示例：插入 [5, 3, 7, 3, 6]，刪除 [3, 7]
    Node* root = build({5, 3, 7, 3, 6});
    assert(inorder(root) == vector<ll>({3, 3, 5, 6, 7}));
    assert(height(root) == 3);
    assert(treeSize(root) == 5);
    assert(isValid(root));
    assert(search(root, 3) && search(root, 7) && !search(root, 9));
    root = deleteNode(root, 3);
    assert(inorder(root) == vector<ll>({3, 5, 6, 7}));   // 重複鍵值只刪掉一個
    assert(isValid(root));
    root = deleteNode(root, 7);                          // 兩個孩子：用中序後繼頂替
    assert(inorder(root) == vector<ll>({3, 5, 6}));
    assert(isValid(root));
    assert(height(root) == 2);
    freeTree(root);

    // 空樹
    assert(inorder(nullptr).empty());
    assert(height(nullptr) == 0);
    assert(treeSize(nullptr) == 0);
    assert(isValid(nullptr));
    assert(deleteNode(nullptr, 1) == nullptr);

    // 單節點
    Node* one = build({42});
    assert(inorder(one) == vector<ll>({42}));
    assert(height(one) == 1);
    assert(deleteNode(one, 42) == nullptr);

    // 有序插入 → 退化成鏈（樸素 BST 的軟肋）
    Node* chain = build({1, 2, 3, 4, 5});
    assert(inorder(chain) == vector<ll>({1, 2, 3, 4, 5}));
    assert(height(chain) == 5);
    assert(isValid(chain));
    freeTree(chain);
    Node* rnd = build({5, 3, 8, 1, 4});
    assert(height(rnd) == 3);
    freeTree(rnd);

    // 刪根 / 刪葉子 / 刪不存在的鍵
    Node* t = build({4, 2, 6, 1, 3, 5, 7});
    assert(height(t) == 3);
    t = deleteNode(t, 4);                                // 刪根：中序後繼是 5
    assert(inorder(t) == vector<ll>({1, 2, 3, 5, 6, 7}));
    assert(isValid(t));
    t = deleteNode(t, 1);                                // 刪葉子
    assert(inorder(t) == vector<ll>({2, 3, 5, 6, 7}));
    t = deleteNode(t, 100);                              // 鍵不存在
    assert(inorder(t) == vector<ll>({2, 3, 5, 6, 7}));
    assert(isValid(t));
    freeTree(t);

    // 全部刪空
    Node* e = build({2, 1, 3});
    for (ll v : {2LL, 1LL, 3LL}) {
        e = deleteNode(e, v);
        assert(isValid(e));
    }
    assert(e == nullptr);

    // 重複鍵值：全相同
    Node* dup = build({7, 7, 7});
    assert(inorder(dup) == vector<ll>({7, 7, 7}));
    assert(height(dup) == 3);
    assert(isValid(dup));
    dup = deleteNode(dup, 7);
    assert(inorder(dup) == vector<ll>({7, 7}));
    assert(isValid(dup));
    freeTree(dup);

    mt19937 rng(20260929);

    // 隨機對拍：用有序 vector 模擬多重集合，每一步都與樹對比
    for (int iter = 0; iter < 300; ++iter) {
        int n = rng() % 61;
        vector<ll> inserts;
        for (int i = 0; i < n; ++i) inserts.push_back((ll)(rng() % 17) - 8);
        Node* rt = build(inserts);
        vector<ll> model = inserts;
        sort(model.begin(), model.end());
        assert(inorder(rt) == model);
        vector<ll> rec;
        inorderRec(rt, rec);
        assert(rec == model);
        assert(isValid(rt));
        assert(treeSize(rt) == (int)model.size());
        assert(height(rt) == heightBfs(rt));
        if (n > 0) assert((int)ceil(log2(n + 1)) <= height(rt) && height(rt) <= n);

        int m = rng() % 21;
        vector<ll> deletes;
        for (int i = 0; i < m; ++i) deletes.push_back((ll)(rng() % 17) - 8);
        for (ll v : deletes) {
            bool inModel = find(model.begin(), model.end(), v) != model.end();
            assert(search(rt, v) == inModel);
            rt = deleteNode(rt, v);
            if (inModel) model.erase(find(model.begin(), model.end(), v));  // 只刪一個副本
            assert(inorder(rt) == model);
            assert(isValid(rt));
            assert(treeSize(rt) == (int)model.size());
            assert(height(rt) == heightBfs(rt));
        }
        freeTree(rt);
    }

    cout << "all tests passed" << '\n';
}

int main() {
    // 無 stdin（或 stdin 爲空）時跑內置測試
    string data, line;
    bool hasInput = false;
    while (getline(cin, line)) {
        data += line;
        data += '\n';
        if (!line.empty()) hasInput = true;
    }
    if (hasInput) runIo(data);
    else runTests();
    return 0;
}
