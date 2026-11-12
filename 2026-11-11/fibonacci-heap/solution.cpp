// 費波那契堆（Fibonacci heap）：可合併堆，支援攤還 O(1) 的插入與 decrease-key。
//
// 命令（IO 模式，從 stdin 讀入）：
//   push k      -> 插入鍵值 k，輸出分配到的 id（從 0 遞增）
//   popmin      -> 彈出最小鍵值並輸出；空堆輸出 -1
//   min         -> 輸出當前最小鍵值（不彈出）；空堆輸出 -1
//   dec id nk   -> 把 id 的鍵值降為 nk（僅當 nk <= 現值且 id 仍存活，否則忽略）
//   size        -> 輸出當前元素個數
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * push O(1)：掛進根串列，必要時更新 minRoot。
//   * dec 攤還 O(1)：破壞堆序就 cut（剪回根串列），父節點失去第二個
//     孩子時級聯 cut。
//   * popmin 攤還 O(log n)：孩子升為根後 consolidate，依 degree 合併
//     同度數的樹，保證每個度數至多一棵（第 k 小的樹至少 F(k) 個節點）。
//   * popmin 只輸出鍵值，但 dec 之後 id 的死活依賴「平手時彈哪個節點」。
//     因此堆序採字典序 (key, nid)：父節點嚴格小於孩子，全域最小節點
//     必為根，平手彈出 id 最小者——與 Python 版 / 參考實作完全一致。
//   * 應用：Dijkstra / Prim 的 O(E + V log V) 理論最優實作。

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <climits>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

using namespace std;

typedef long long ll;

struct Node {
    ll key;
    int nid;
    int deg;
    bool mark;
    bool alive;
    Node* parent;
    vector<Node*> children;
};

static Node* newNode(ll key, int nid) {
    Node* n = new Node();
    n->key = key;
    n->nid = nid;
    n->deg = 0;
    n->mark = false;
    n->alive = true;
    n->parent = NULL;
    return n;
}

// 字典序 (key, nid) 比較：堆序的基本序，平手時取 id 小者
static bool lexLess(Node* a, Node* b) {
    return a->key < b->key || (a->key == b->key && a->nid < b->nid);
}

struct FibonacciHeap {
    vector<Node*> roots;
    Node* minRoot;
    vector<Node*> nodes;   // id -> Node（含已彈出的死節點）
    int cnt;

    FibonacciHeap() : minRoot(NULL), cnt(0) {}

    // 把根 y 掛到根 x 下面（(y.key,y.nid) >= (x.key,x.nid)）
    void link(Node* y, Node* x) {
        y->parent = x;
        x->children.push_back(y);
        x->deg += 1;
        y->mark = false;
    }

    int push(ll key) {
        Node* node = newNode(key, cnt);
        cnt++;
        nodes.push_back(node);
        roots.push_back(node);
        if (minRoot == NULL || lexLess(node, minRoot)) minRoot = node;
        return node->nid;
    }

    void consolidate() {
        map<int, Node*> byDeg;
        for (size_t t = 0; t < roots.size(); ++t) {
            Node* x = roots[t];
            int d = x->deg;
            while (byDeg.count(d)) {
                Node* y = byDeg[d];
                byDeg.erase(d);
                if (lexLess(y, x)) {
                    Node* tmp = x;
                    x = y;
                    y = tmp;
                }
                link(y, x);
                d += 1;
            }
            byDeg[d] = x;
        }
        roots.clear();
        minRoot = NULL;
        for (map<int, Node*>::iterator it = byDeg.begin(); it != byDeg.end(); ++it) {
            roots.push_back(it->second);
            if (minRoot == NULL || lexLess(it->second, minRoot)) minRoot = it->second;
        }
    }

    ll popmin() {
        Node* z = minRoot;
        if (z == NULL) return -1;
        z->alive = false;
        // 最小根的孩子全部升為根
        vector<Node*> newRoots;
        for (size_t t = 0; t < roots.size(); ++t) {
            if (roots[t] != z) newRoots.push_back(roots[t]);
        }
        for (size_t t = 0; t < z->children.size(); ++t) {
            z->children[t]->parent = NULL;
            z->children[t]->mark = false;  // 根不該帶標記
            newRoots.push_back(z->children[t]);
        }
        z->children.clear();
        roots = newRoots;
        consolidate();
        return z->key;
    }

    void cut(Node* x) {
        Node* p = x->parent;
        for (size_t t = 0; t < p->children.size(); ++t) {
            if (p->children[t] == x) {
                p->children.erase(p->children.begin() + t);
                break;
            }
        }
        p->deg -= 1;
        x->parent = NULL;
        x->mark = false;
        roots.push_back(x);
    }

    void cascadingCut(Node* y) {
        Node* z = y->parent;
        if (z == NULL) return;
        if (!y->mark) {
            y->mark = true;
        } else {
            cut(y);
            cascadingCut(z);
        }
    }

    void dec(int nid, ll nk) {
        if (nid < 0 || nid >= (int)nodes.size()) return;
        Node* node = nodes[nid];
        if (!node->alive || nk > node->key) return;
        node->key = nk;
        Node* p = node->parent;
        // 堆序是字典序 (key, nid)：與父節點「平手但 id 更小」也算破壞堆序
        if (p != NULL && lexLess(node, p)) {
            cut(node);
            cascadingCut(p);
        }
        if (minRoot == NULL || lexLess(node, minRoot)) minRoot = node;
    }

    ll peekMin() const { return minRoot == NULL ? -1 : minRoot->key; }

    int size() const {
        int s = 0;
        for (size_t t = 0; t < nodes.size(); ++t) {
            if (nodes[t]->alive) s++;
        }
        return s;
    }
};

// 線性掃描參考實作（僅供內建測試對拍）
struct RefHeap {
    map<int, ll> alive;
    int cnt;
    RefHeap() : cnt(0) {}

    int push(ll key) {
        alive[cnt] = key;
        cnt++;
        return cnt - 1;
    }
    ll popmin() {
        if (alive.empty()) return -1;
        int bestId = -1;
        ll bestV = 0;
        for (map<int, ll>::iterator it = alive.begin(); it != alive.end(); ++it) {
            if (bestId < 0 || it->second < bestV) {
                bestId = it->first;
                bestV = it->second;
            }
        }
        alive.erase(bestId);
        return bestV;
    }
    ll peekMin() const {
        if (alive.empty()) return -1;
        ll best = LLONG_MAX;
        for (map<int, ll>::const_iterator it = alive.begin(); it != alive.end(); ++it) {
            if (it->second < best) best = it->second;
        }
        return best;
    }
    void dec(int nid, ll nk) {
        map<int, ll>::iterator it = alive.find(nid);
        if (it != alive.end() && nk <= it->second) it->second = nk;
    }
    int size() const { return (int)alive.size(); }
};

static bool tryParse(const string& s, ll& v) {
    if (s.empty()) return false;
    size_t p = 0;
    if (s[0] == '+' || s[0] == '-') p = 1;
    if (p >= s.size()) return false;
    for (size_t t = p; t < s.size(); ++t) {
        if (s[t] < '0' || s[t] > '9') return false;
    }
    v = strtoll(s.c_str(), NULL, 10);
    return true;
}

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    string data;
    char buf[65536];
    size_t got;
    while ((got = fread(buf, 1, sizeof(buf), stdin)) > 0) data.append(buf, got);

    bool hasInput = false;
    for (char c : data) {
        if (!isspace((unsigned char)c)) {
            hasInput = true;
            break;
        }
    }

    if (hasInput) {
        vector<string> toks;
        {
            string cur;
            for (char c : data) {
                if (isspace((unsigned char)c)) {
                    if (!cur.empty()) {
                        toks.push_back(cur);
                        cur.clear();
                    }
                } else {
                    cur += c;
                }
            }
            if (!cur.empty()) toks.push_back(cur);
        }
        size_t pos = 0;
        auto nxt = [&](string& outTok) -> bool {
            if (pos >= toks.size()) return false;
            outTok = toks[pos++];
            return true;
        };
        FibonacciHeap heap;
        string cmd;
        while (nxt(cmd)) {
            if (cmd == "push") {
                string t;
                ll k;
                if (!nxt(t) || !tryParse(t, k)) break;
                printf("%d\n", heap.push(k));
            } else if (cmd == "popmin") {
                printf("%lld\n", heap.popmin());
            } else if (cmd == "min") {
                printf("%lld\n", heap.peekMin());
            } else if (cmd == "dec") {
                string t;
                ll id, k;
                if (!nxt(t) || !tryParse(t, id)) break;
                if (!nxt(t) || !tryParse(t, k)) break;
                heap.dec((int)id, k);
            } else if (cmd == "size") {
                printf("%d\n", heap.size());
            } else {
                break;
            }
        }
        fflush(stdout);
        return 0;
    }

    // ---- 內建斷言測試 ----
    {
        FibonacciHeap h;
        if (h.push(5) != 0) return 1;
        if (h.push(3) != 1) return 1;
        if (h.push(8) != 2) return 1;
        if (h.peekMin() != 3) return 1;
        if (h.size() != 3) return 1;
        if (h.popmin() != 3) return 1;
        if (h.peekMin() != 5) return 1;
        h.dec(2, 1);
        if (h.peekMin() != 1) return 1;
        if (h.popmin() != 1) return 1;
        if (h.popmin() != 5) return 1;
        if (h.popmin() != -1) return 1;
        if (h.size() != 0) return 1;
    }
    {
        FibonacciHeap h;
        if (h.popmin() != -1) return 1;
        if (h.peekMin() != -1) return 1;
        if (h.size() != 0) return 1;
        h.dec(99, 0);  // 非法 id 忽略
        if (h.size() != 0) return 1;
    }
    {
        FibonacciHeap h;
        vector<int> ids;
        for (int i = 0; i < 50; ++i) ids.push_back(h.push(i % 7));
        if (h.size() != 50) return 1;
        for (int i = 0; i < 50; i += 5) h.dec(ids[i], -i);
        if (h.peekMin() != -45) return 1;
        vector<ll> vals;
        for (int i = 0; i < 50; ++i) vals.push_back(h.popmin());
        vector<ll> expect;
        for (int i = 0; i < 50; ++i) expect.push_back(i % 5 == 0 ? -i : i % 7);
        sort(vals.begin(), vals.end());
        sort(expect.begin(), expect.end());
        if (vals != expect) return 1;
        if (h.size() != 0) return 1;
    }
    {
        // 與線性掃描參考實作對拍
        unsigned int st = 20261111u;
        auto rnd = [&](unsigned int mm) -> unsigned int {
            st = st * 1103515245u + 12345u;
            return (st >> 1) % mm;
        };
        for (int tc = 0; tc < 300; ++tc) {
            FibonacciHeap h;
            RefHeap ref;
            int ops = (int)rnd(200) + 30;
            for (int t = 0; t < ops; ++t) {
                unsigned int op = rnd(10);
                if (op < 4) {
                    ll k = (ll)rnd(20);
                    if (h.push(k) != ref.push(k)) return 1;
                } else if (op < 6) {
                    if (h.popmin() != ref.popmin()) return 1;
                } else if (op < 7) {
                    if (h.peekMin() != ref.peekMin()) return 1;
                } else if (op < 9) {
                    if (h.cnt > 0) {
                        int nid = (int)rnd((unsigned int)h.cnt);
                        ll k = (ll)rnd(20);
                        h.dec(nid, k);
                        ref.dec(nid, k);
                    }
                } else {
                    if (h.size() != ref.size()) return 1;
                }
            }
            while (ref.size() > 0) {
                if (h.popmin() != ref.popmin()) return 1;
            }
            if (h.popmin() != -1 || ref.popmin() != -1) return 1;
        }
    }

    printf("all tests passed\n");
    return 0;
}
