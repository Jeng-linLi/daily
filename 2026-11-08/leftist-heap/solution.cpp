// 左偏樹 Leftist Heap（可合併堆 Mergeable Heap）
//
// 題意：維護 10 個（編號 0..9）可合併的小根堆，支援五種操作：
//   1 h x   把值 x 插入堆 h
//   2 h     彈出堆 h 的最小值並輸出；空堆輸出 -1
//   3 h     詢問堆 h 的最小值（不彈出）；空堆輸出 -1
//   4 a b   把堆 b 併入堆 a（之後 b 為空；a == b 時視為無操作）
//   5 h     輸出堆 h 的元素個數
//
// 要點：dist(空)=-1、葉子 dist=0、dist(x)=dist(右子)+1；
//       左偏性質：左兒子 dist >= 右兒子 dist ⟹ 樹高 O(log n)；
//       merge 值大的根掛到值小的根右子樹遞歸合併，回溯時必要時交換左右。
//       push/pop/top/merge 全部 O(log n)。
//
// IO 模式（stdin 有內容時）：q 行操作，輸出每個 2/3/5 操作的結果，一行一個。
// 隨機測試保證同時刻各堆內值互異（輸出順序才能唯一）。
// 無 stdin 時跑內建斷言（與有序陣列參考實作對拍 + 不變量實走檢查）
// 並輸出 all tests passed。
#include <cstdio>
#include <cctype>
#include <cstring>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;
typedef long long ll;

struct Node {
    ll val;
    int dist;
    Node* l;
    Node* r;
    Node(ll v) : val(v), dist(0), l(NULL), r(NULL) {}
};

static Node* meld(Node* a, Node* b) {
    if (!a) return b;
    if (!b) return a;
    if (b->val < a->val) swap(a, b);
    a->r = meld(a->r, b);
    if (!a->l || (a->r && a->l->dist < a->r->dist)) swap(a->l, a->r);
    a->dist = a->r ? a->r->dist + 1 : 0;
    return a;
}

struct LeftistHeap {
    Node* root;
    int size;
    LeftistHeap() : root(NULL), size(0) {}
    void push(ll x) {
        root = meld(root, new Node(x));
        size++;
    }
    ll pop() {
        if (!root) return -1;
        ll v = root->val;
        Node* nr = meld(root->l, root->r);
        delete root;
        root = nr;
        size--;
        return v;
    }
    ll top() const { return root ? root->val : -1; }
};

static int checkInvariants(const Node* u) {
    if (!u) return 0;
    int dl = u->l ? u->l->dist : -1;  // dist(空) = -1
    int dr = u->r ? u->r->dist : -1;
    assert(!u->l || u->l->val >= u->val);  // 堆序性
    assert(!u->r || u->r->val >= u->val);
    assert(dl >= dr);                       // 左偏性質
    assert(u->dist == dr + 1);              // dist 欄位一致
    return 1 + checkInvariants(u->l) + checkInvariants(u->r);
}

// 參考實作：有序陣列（測試規模足夠）
struct RefHeap {
    vector<ll> a;
    void push(ll x) { a.insert(lower_bound(a.begin(), a.end(), x), x); }
    ll pop() {
        if (a.empty()) return -1;
        ll v = a[0];
        a.erase(a.begin());
        return v;
    }
    ll top() const { return a.empty() ? -1 : a[0]; }
};

static void runAllTests() {
    // 1) 空堆行為
    {
        LeftistHeap h;
        assert(h.pop() == -1 && h.top() == -1 && h.size == 0);
    }
    // 2) 固定情境（值互異）
    {
        vector<LeftistHeap> hl(10);
        vector<RefHeap> hr(10);
        struct Op { int c, h, x; };
        vector<Op> ops = {
            {1,0,5},{1,0,2},{1,0,8},{3,0,0},{2,0,0},{5,0,0},
            {1,1,7},{4,0,1},{5,0,0},{5,1,0},
            {2,0,0},{2,0,0},{2,0,0}};
        vector<ll> outL, outR;
        for (size_t i = 0; i < ops.size(); i++) {
            const Op& o = ops[i];
            if (o.c == 1) { hl[o.h].push(o.x); hr[o.h].push(o.x); }
            else if (o.c == 2) { outL.push_back(hl[o.h].pop()); outR.push_back(hr[o.h].pop()); }
            else if (o.c == 3) { outL.push_back(hl[o.h].top()); outR.push_back(hr[o.h].top()); }
            else if (o.c == 4) {
                if (o.h != o.x) {
                    hl[o.h].root = meld(hl[o.h].root, hl[o.x].root);
                    hl[o.h].size += hl[o.x].size;
                    hl[o.x].root = NULL; hl[o.x].size = 0;
                    // 參考版：合併兩個有序陣列
                    vector<ll> m;
                    merge(hr[o.h].a.begin(), hr[o.h].a.end(),
                          hr[o.x].a.begin(), hr[o.x].a.end(), back_inserter(m));
                    hr[o.h].a = m; hr[o.x].a.clear();
                }
            } else { outL.push_back(hl[o.h].size); outR.push_back((ll)hr[o.h].a.size()); }
        }
        assert(outL == outR);
        for (int i = 0; i < 10; i++)
            if (hl[i].root) assert(checkInvariants(hl[i].root) == hl[i].size);
        // 實跑取值：top=2、pop=2、size=2 → 合併後 size=3、h1=0 → pop 5,7,8
        assert(outL == vector<ll>({2, 2, 2, 3, 0, 5, 7, 8}));
        assert(hl[0].size == 0 && hl[0].pop() == -1);
    }
    // 3) 自製 LCG 隨機對拍 + 不變量檢查
    unsigned long long seed = 20261108ULL;
    auto rnd = [&](int n) -> int {
        seed = seed * 1103515245ULL + 12345ULL;
        seed &= 0x7fffffffULL;
        return (int)(seed % (unsigned)n);
    };
    for (int tc = 0; tc < 120; tc++) {
        vector<LeftistHeap> hl(10);
        vector<RefHeap> hr(10);
        // 全域值互異：從 1 開始線性取值再映射（值域內互異即可）
        vector<ll> pool;
        for (int i = 0; i < 400; i++) pool.push_back((ll)rnd(2000001) - 1000000);
        sort(pool.begin(), pool.end());
        pool.erase(unique(pool.begin(), pool.end()), pool.end());
        size_t pi = 0;
        int q = 60 + rnd(120);
        vector<ll> outL, outR;
        for (int i = 0; i < q; i++) {
            int c = rnd(100);
            if (c < 45) {
                int h = rnd(10);
                ll v = pool[pi++];
                hl[h].push(v); hr[h].push(v);
            } else if (c < 70) {
                int h = rnd(10);
                outL.push_back(hl[h].pop()); outR.push_back(hr[h].pop());
            } else if (c < 85) {
                int h = rnd(10);
                outL.push_back(hl[h].top()); outR.push_back(hr[h].top());
            } else if (c < 95) {
                int a = rnd(10), b = rnd(10);
                if (a != b) {
                    hl[a].root = meld(hl[a].root, hl[b].root);
                    hl[a].size += hl[b].size;
                    hl[b].root = NULL; hl[b].size = 0;
                    vector<ll> m;
                    merge(hr[a].a.begin(), hr[a].a.end(),
                          hr[b].a.begin(), hr[b].a.end(), back_inserter(m));
                    hr[a].a = m; hr[b].a.clear();
                }
            } else {
                int h = rnd(10);
                outL.push_back(hl[h].size); outR.push_back((ll)hr[h].a.size());
            }
        }
        assert(outL == outR);
        for (int i = 0; i < 10; i++)
            if (hl[i].root) assert(checkInvariants(hl[i].root) == hl[i].size);
    }
}

int main() {
    stringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    bool anyTok = false;
    for (size_t i = 0; i < input.size(); i++)
        if (!isspace((unsigned char)input[i])) { anyTok = true; break; }
    if (!anyTok) { runAllTests(); printf("all tests passed\n"); return 0; }
    istringstream it(input);
    long long q;
    if (!(it >> q)) { return 0; }
    vector<LeftistHeap> heaps(10);
    ostringstream oss;
    for (long long i = 0; i < q; i++) {
        int c;
        it >> c;
        if (c == 1) {
            int h; ll x;
            it >> h >> x;
            if (0 <= h && h <= 9) heaps[h].push(x);
        } else if (c == 2) {
            int h; it >> h;
            oss << (0 <= h && h <= 9 ? heaps[h].pop() : -1) << '\n';
        } else if (c == 3) {
            int h; it >> h;
            oss << (0 <= h && h <= 9 ? heaps[h].top() : -1) << '\n';
        } else if (c == 4) {
            int a, b;
            it >> a >> b;
            if (0 <= a && a <= 9 && 0 <= b && b <= 9 && a != b) {
                heaps[a].root = meld(heaps[a].root, heaps[b].root);
                heaps[a].size += heaps[b].size;
                heaps[b].root = NULL; heaps[b].size = 0;
            }
        } else if (c == 5) {
            int h; it >> h;
            oss << (0 <= h && h <= 9 ? heaps[h].size : 0) << '\n';
        }
    }
    fputs(oss.str().c_str(), stdout);
    return 0;
}
