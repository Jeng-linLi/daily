// -*- coding: utf-8 -*-
// 配對堆 Pairing Heap（two-pass combine / 支援 decrease-key）
//
// 題意：維護 10 個（編號 0..9）可合併的小根堆，支援六種操作：
//   1 h x     插入值 x
//   2 h       彈出最小值並輸出；空堆輸出 -1
//   3 h       詢問最小值（不彈出）；空堆輸出 -1
//   4 a b     把堆 b 併入堆 a（a == b 視為無操作）
//   5 h       輸出堆 h 的元素個數
//   6 h x y   把堆 h 中的值 x 改為 y（decrease-key，無輸出）；
//             x 不在堆 h、y 已存在或 x == y 時視為無操作。
//
// 節點用 child（第一個兒子）+ sibling（下一個兄弟）表示。
//   - meld：值大的根作為值小根的第一個兒子，O(1)。
//     par 指向父親：掛新兒子時舊兒子的 par 仍是父親，不動。
//   - delete-min：兒子先兩兩配對（第一遍），再從右往左依序合併（第二遍）。
//   - decrease-key：從父親的兒子鏈摘下 → 改值 → meld 回根。
//     根節點要先把兒子鏈 two-pass 重組、改值後把自己 meld 回去，
//     否則原地改大 root.val 會破壞堆序。
//
// IO 模式（stdin 有內容時）：q 行操作，輸出每個 2/3/5 操作的結果。
// 值域 |x| ≤ 10^9；同一時刻全域值互異（decrease-key 靠值定位節點）。
// 無 stdin 時跑內建斷言（與有序陣列參考實作對拍 + 堆序不變量檢查），
// 輸出 all tests passed。
#include <cstdio>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;
typedef long long ll;
typedef unsigned long long ull;

static ull g_seed;

static inline unsigned rnd(unsigned n) {
    // 與 Python 版相同的 LCG（mod 2^31），保證兩語言測試序列一致
    g_seed = (g_seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
    return (unsigned)(g_seed % n);
}

struct Node {
    ll val;
    Node *child, *sib, *par;
    Node(ll v) : val(v), child(nullptr), sib(nullptr), par(nullptr) {}
};

static Node* meld(Node* a, Node* b) {
    if (!a) return b;
    if (!b) return a;
    if (b->val < a->val) swap(a, b);
    // b 掛到 a 的兒子鏈頭；a 的舊兒子 par 仍是 a，不動
    b->sib = a->child;
    b->par = a;
    a->child = b;
    return a;
}

// 把兒子鏈摘成 list，並清掉 par/sib 指針
static vector<Node*> collect(Node* head) {
    vector<Node*> out;
    for (Node* c = head; c;) {
        Node* nxt = c->sib;
        c->sib = nullptr;
        c->par = nullptr;
        out.push_back(c);
        c = nxt;
    }
    return out;
}

// 兩遍合併：第一遍兩兩配對，第二遍從右往左依序合併
static Node* twoPassCombine(const vector<Node*>& children) {
    if (children.empty()) return nullptr;
    vector<Node*> pairs;
    for (size_t i = 0; i + 1 < children.size(); i += 2)
        pairs.push_back(meld(children[i], children[i + 1]));
    if (children.size() % 2 == 1)
        pairs.push_back(children.back());
    Node* cur = pairs.back();
    for (int i = (int)pairs.size() - 2; i >= 0; --i)
        cur = meld(pairs[i], cur);
    return cur;
}

struct PairingHeap {
    Node* root;
    int size;
    map<ll, Node*> pos;  // 值 -> 節點（同時刻全域值互異）

    PairingHeap() : root(nullptr), size(0) {}

    void push(ll x) {
        Node* node = new Node(x);
        root = meld(root, node);
        pos[x] = node;
        ++size;
    }

    ll pop() {
        if (!root) return -1;
        Node* node = root;
        ll v = node->val;
        root = twoPassCombine(collect(node->child));
        pos.erase(v);
        --size;
        delete node;
        return v;
    }

    ll top() const { return root ? root->val : -1; }

    void decreaseKey(ll x, ll y) {
        if (x == y || pos.find(x) == pos.end() || pos.count(y)) return;
        Node* node = pos[x];
        pos.erase(x);
        if (node != root) {
            // 從父親的兒子鏈上摘下 node（兄弟的 par 仍指向父親，無需改動）
            Node* par = node->par;
            if (par->child == node) {
                par->child = node->sib;
            } else {
                Node* prev = par->child;
                while (prev->sib != node) prev = prev->sib;
                prev->sib = node->sib;
            }
        }
        // 兒子鏈摘下重組：調大時兒子 ≥ 舊值卻可能 < 新值，
        // 子樹內部堆序被破壞，必須重建；調小時重組無害。
        // 摘完清掉 child/sib/par，防止 meld 掛到陳舊鏈。
        vector<Node*> children = collect(node->child);
        node->child = nullptr;
        node->sib = nullptr;
        node->par = nullptr;
        node->val = y;
        Node* tree = meld(twoPassCombine(children), node);
        if (node == root) root = tree;
        else root = meld(root, tree);
        pos[y] = node;
    }

    void meldInto(PairingHeap& other) {
        root = meld(root, other.root);
        size += other.size;
        for (map<ll, Node*>::iterator it = other.pos.begin();
             it != other.pos.end(); ++it)
            pos[it->first] = it->second;  // b 的條目覆蓋，與 Python dict.update 一致
        other.root = nullptr;
        other.size = 0;
        other.pos.clear();
    }
};

// 實走檢查：堆序（父 ≤ 子）+ 節點計數 + par 一致；回傳節點數
static int checkInvariants(Node* root) {
    if (!root) return 0;
    int cnt = 0;
    vector<pair<Node*, ll>> st;
    st.push_back(make_pair(root, root->val));  // 根的下界取自身值，恆過
    while (!st.empty()) {
        Node* u = st.back().first;
        ll lo = st.back().second;
        st.pop_back();
        ++cnt;
        if (lo > u->val) {
            fprintf(stderr, "heap order violated\n");
            exit(1);
        }
        for (Node* c = u->child; c; c = c->sib) {
            if (c->par != u) {
                fprintf(stderr, "parent pointer broken\n");
                exit(1);
            }
            st.push_back(make_pair(c, u->val));
        }
    }
    return cnt;
}

struct RefHeap {
    vector<ll> a;  // 有序陣列
    void push(ll x) { a.insert(lower_bound(a.begin(), a.end(), x), x); }
    ll pop() {
        if (a.empty()) return -1;
        ll v = a[0];
        a.erase(a.begin());
        return v;
    }
    ll top() const { return a.empty() ? -1 : a[0]; }
    void meldInto(RefHeap& other) {
        a.insert(a.end(), other.a.begin(), other.a.end());
        sort(a.begin(), a.end());
        other.a.clear();
    }
};

static void runOps(const vector<vector<ll>>& ops, vector<PairingHeap>& hl,
                   vector<RefHeap>& hr, vector<ll>& outL, vector<ll>& outR) {
    for (size_t k = 0; k < ops.size(); ++k) {
        int c = (int)ops[k][0];
        if (c == 1) {
            int h = (int)ops[k][1];
            ll x = ops[k][2];
            hl[h].push(x);
            hr[h].push(x);
        } else if (c == 2) {
            int h = (int)ops[k][1];
            outL.push_back(hl[h].pop());
            outR.push_back(hr[h].pop());
        } else if (c == 3) {
            int h = (int)ops[k][1];
            outL.push_back(hl[h].top());
            outR.push_back(hr[h].top());
        } else if (c == 4) {
            int a = (int)ops[k][1], b = (int)ops[k][2];
            if (a != b) {
                hl[a].meldInto(hl[b]);
                hr[a].meldInto(hr[b]);
            }
        } else if (c == 5) {
            int h = (int)ops[k][1];
            outL.push_back(hl[h].size);
            outR.push_back((ll)hr[h].a.size());
        } else {
            int h = (int)ops[k][1];
            ll x = ops[k][2], y = ops[k][3];
            hl[h].decreaseKey(x, y);
            // 參考實作：值互異，直接 remove + insort
            vector<ll>& arr = hr[h].a;
            vector<ll>::iterator it = lower_bound(arr.begin(), arr.end(), x);
            if (x != y && it != arr.end() && *it == x &&
                !binary_search(arr.begin(), arr.end(), y)) {
                arr.erase(it);
                arr.insert(lower_bound(arr.begin(), arr.end(), y), y);
            }
        }
    }
}

static void die(const char* msg, int tc) {
    fprintf(stderr, "%s (case %d)\n", msg, tc);
    exit(1);
}

static void runTests() {
    g_seed = 20261109ULL;

    // 1) 空堆行為
    {
        PairingHeap h;
        if (!(h.pop() == -1 && h.top() == -1 && h.size == 0))
            die("empty heap behavior", -1);
        if (checkInvariants(h.root) != 0) die("empty invariant", -1);
    }

    // 2) 固定情境（全域值互異；期望值由實跑取得）
    {
        vector<PairingHeap> hl(10);
        vector<RefHeap> hr(10);
        vector<vector<ll>> ops = {
            {1, 0, 5}, {1, 0, 2}, {1, 0, 8}, {3, 0}, {2, 0}, {5, 0},
            {6, 0, 8, 1}, {3, 0}, {1, 1, 7}, {4, 0, 1}, {5, 0}, {5, 1},
            {2, 0}, {2, 0}, {2, 0}, {2, 0}};
        vector<ll> outL, outR;
        runOps(ops, hl, hr, outL, outR);
        if (outL != outR) die("fixed scenario mismatch vs reference", -1);
        // 實跑取值：top=2, pop=2, size=2；decrease 8→1 後 top=1；
        // 合併後 size=3, h1=0；pop 依序 1, 5, 7, 最後空堆 -1
        vector<ll> want = {2, 2, 2, 1, 3, 0, 1, 5, 7, -1};
        if (outL != want) die("fixed scenario expected values", -1);
        for (int i = 0; i < 10; ++i)
            if (checkInvariants(hl[i].root) != hl[i].size)
                die("fixed invariant size", -1);
    }

    // 3) 自製 LCG 隨機對拍（全域值互異）+ 不變量檢查
    for (int tc = 0; tc < 120; ++tc) {
        vector<PairingHeap> hl(10);
        vector<RefHeap> hr(10);
        vector<vector<ll>> ops;
        int q = 60 + rnd(120);
        map<ll, int> owner;  // 全域值 -> 所在堆（生成 decrease-key 用）
        for (int t = 0; t < q; ++t) {
            int c = rnd(100);
            if (c < 40) {
                ll v = (ll)rnd(2000001) - 1000000;
                while (owner.count(v)) v = (ll)rnd(2000001) - 1000000;
                int h = rnd(10);
                ops.push_back({1, h, v});
                owner[v] = h;
            } else if (c < 62) {
                int h = rnd(10);
                ops.push_back({2, h});
            } else if (c < 77) {
                ops.push_back({3, rnd(10)});
            } else if (c < 87) {
                int a = rnd(10), b = rnd(10);
                ops.push_back({4, a, b});
            } else if (c < 92) {
                ops.push_back({5, rnd(10)});
            } else {
                int h = rnd(10);
                // 找堆 h 裡已存在的值（沿 owner 掃描一份快照）
                vector<ll> vals;
                for (map<ll, int>::iterator it = owner.begin();
                     it != owner.end(); ++it)
                    if (it->second == h) vals.push_back(it->first);
                if (vals.empty()) {
                    ll v = (ll)rnd(2000001) - 1000000;
                    while (owner.count(v)) v = (ll)rnd(2000001) - 1000000;
                    ops.push_back({1, h, v});
                    owner[v] = h;
                } else {
                    ll x = vals[rnd((unsigned)vals.size())];
                    ll y = (ll)rnd(2000001) - 1000000;
                    while (owner.count(y)) y = (ll)rnd(2000001) - 1000000;
                    ops.push_back({6, h, x, y});
                    owner[x] = -1;  // 移除舊值標記
                    owner[y] = h;
                }
            }
        }
        vector<ll> outL, outR;
        runOps(ops, hl, hr, outL, outR);
        if (outL != outR) die("random scenario mismatch", tc);
        for (int i = 0; i < 10; ++i)
            if (checkInvariants(hl[i].root) != hl[i].size)
                die("invariant size mismatch", tc);
    }
}

int main() {
    long long first;
    if (scanf("%lld", &first) != 1) {
        runTests();
        printf("all tests passed\n");
        return 0;
    }
    int n = (int)first;
    int q = n;
    vector<PairingHeap> heaps(10);
    bool firstOut = true;
    for (int t = 0; t < q; ++t) {
        int c;
        if (scanf("%d", &c) != 1) break;
        if (c == 1) {
            int h;
            long long x;
            if (scanf("%d %lld", &h, &x) != 2) break;
            if (0 <= h && h <= 9) heaps[h].push(x);
        } else if (c == 2) {
            int h;
            if (scanf("%d", &h) != 1) break;
            ll v = (0 <= h && h <= 9) ? heaps[h].pop() : -1;
            printf(firstOut ? "%lld" : "\n%lld", v);
            firstOut = false;
        } else if (c == 3) {
            int h;
            if (scanf("%d", &h) != 1) break;
            ll v = (0 <= h && h <= 9) ? heaps[h].top() : -1;
            printf(firstOut ? "%lld" : "\n%lld", v);
            firstOut = false;
        } else if (c == 4) {
            int a, b;
            if (scanf("%d %d", &a, &b) != 2) break;
            if (0 <= a && a <= 9 && 0 <= b && b <= 9 && a != b)
                heaps[a].meldInto(heaps[b]);
        } else if (c == 5) {
            int h;
            if (scanf("%d", &h) != 1) break;
            ll v = (0 <= h && h <= 9) ? heaps[h].size : 0;
            printf(firstOut ? "%lld" : "\n%lld", v);
            firstOut = false;
        } else if (c == 6) {
            int h;
            long long x, y;
            if (scanf("%d %lld %lld", &h, &x, &y) != 3) break;
            if (0 <= h && h <= 9) heaps[h].decreaseKey(x, y);
        }
    }
    if (!firstOut) printf("\n");
    return 0;
}
