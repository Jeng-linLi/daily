// 開放尋址哈希表（線性 / 二次 / 雙重探測）與 LRU Cache
//
// 題意：
//     實現一個不使用鏈地址法的哈希表（string -> int），支持 put / get / erase，
//     三種探測策略（線性、二次、雙重散列），墓碑刪除，負載因子 > 0.5 自動擴容再散列，
//     並統計探測次數。另外實現 LRU Cache（哈希表 + 雙向鏈表），get / put 均 O(1)。
//
// 思路：
//     ### 為什麼刪除要用墓碑
//     開放尋址的查找沿探測序列走到「第一個 EMPTY 槽」才宣告失敗。若刪除直接設回 EMPTY，
//     會把探測鏈截斷，後方元素永遠查不到。所以刪除標記 DELETED：
//       查找遇到 DELETED 繼續前進（不停），但記住第一個 DELETED 供插入複用；
//       插入一路掃到 EMPTY，再回頭用第一個墓碑。
//     ### 三種探測
//       線性 h+i：緩存友好，但有一次聚集（primary clustering）。
//       二次 h+i*i：消除一次聚集，仍有二次聚集；M 為質數且負載 <= 0.5 時，
//                   i^2 mod M 取到 (M+1)/2 個不同值，與 >= M/2 個空槽必有交集，保證找得到空位。
//       雙重 h+i*(1+h' mod (M-1))：步長與質數 M 互質，可走遍全表，聚集最少。
//     ### LRU Cache
//     哈希表 O(1) 定位，雙向鏈表 O(1) 維護順序：頭哨兵後是最久未使用，尾哨兵前是最近使用。
//       get 命中 -> 摘下掛尾端；put 已存在 -> 改值掛尾端；不存在且已滿 -> 淘汰頭部後第一個。
//
// 輸入格式（stdin，全部以空白分隔）：
//     mode                      linear | quadratic | double
//     cap                       初始容量（提升到 >= 2 的質數）
//     m                         哈希表操作數
//     m 行：put key value | get key | del key
//     lru_cap                   LRU 容量
//     q                         LRU 操作數
//     q 行：put key value | get key
// 輸出格式（stdout）：
//     第 1 行：哈希表最終容量
//     第 2 行：存活元素個數
//     第 3 行：累計探測次數
//     第 4 行：所有元素，按 key 字典序，k:v 以空白分隔（空表輸出空行）
//     接著每個 get 一行：`1 <value>` 或 `0 -1`
//     再一行：LRU 目前大小
//     再一行：LRU 內容，由最久未使用到最近使用（空則空行）
//     再接著每個 LRU get 一行：<value> 或 -1
// 輸入被截斷時，缺的部分按「0 個操作 / 容量 0」處理。
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// 槽位的三種狀態
enum SlotState { EMPTY = 0, OCCUPIED = 1, DELETED = 2 };

static const uint64_t MASK32 = 0xFFFFFFFFULL;

// ---------------------------------------------------------------- 雜湊函數
// 多項式雜湊（基數 131，截斷到 32 位無號）。用確定性雜湊才能和 Python 版逐字節對拍
// （只對 ASCII 保證一致；Python 是逐碼點，C++ 是逐位元組）。
static uint64_t polyHash(const string& s) {
    uint64_t h = 0;
    for (unsigned char c : s) h = (h * 131 + c) & MASK32;
    return h;
}

static bool isPrime(long long x) {
    if (x < 2) return false;
    if (x % 2 == 0) return x == 2;
    for (long long d = 3; d * d <= x; d += 2)
        if (x % d == 0) return false;
    return true;
}

static long long nextPrime(long long x) {
    if (x <= 2) return 2;
    if (x % 2 == 0) x += 1;
    while (!isPrime(x)) x += 2;
    return x;
}

// ---------------------------------------------------------------- 開放尋址哈希表
class OpenAddressingHashTable {
public:
    int cap;
    string mode;
    vector<string> keys;
    vector<long long> vals;
    vector<int> state;
    int size;     // 存活元素數（不含墓碑）
    int used;     // 被佔用槽位數（含墓碑）
    long long probes;

    OpenAddressingHashTable(int capacity = 8, string m = "linear")
        : cap((int)nextPrime(max(capacity, 1))), mode(m), size(0), used(0), probes(0) {
        if (cap < 2) cap = 2;
        if (mode != "linear" && mode != "quadratic" && mode != "double") mode = "linear";
        keys.assign((size_t)cap, string(""));
        vals.assign((size_t)cap, 0);
        state.assign((size_t)cap, EMPTY);
    }

    // 第 i 次探測的槽位下標
    int probe(const string& key, int i) const {
        long long base = (long long)(polyHash(key) % (uint64_t)cap);
        if (mode == "linear") return (int)((base + i) % cap);
        if (mode == "quadratic") return (int)((base + (long long)i * i) % cap);
        long long step = 1 + (long long)(polyHash(key) % (uint64_t)(cap - 1));
        return (int)((base + (long long)i * step) % cap);
    }

    // 返回 (槽位下標, 是否已存在)；掃到 EMPTY 才判定不存在，中途 DELETED 記錄供插入複用
    pair<int, bool> findSlot(const string& key) {
        int firstTomb = -1;
        for (int i = 0; i < cap; ++i) {
            int idx = probe(key, i);
            probes++;
            int st = state[(size_t)idx];
            if (st == EMPTY) return {(firstTomb >= 0) ? firstTomb : idx, false};
            if (st == DELETED) {
                if (firstTomb < 0) firstTomb = idx;
            } else if (keys[(size_t)idx] == key) {
                return {idx, true};
            }
        }
        return {(firstTomb >= 0) ? firstTomb : -1, false};
    }

    void rehash(int newCap) {
        vector<pair<string, long long>> pairs;
        for (int i = 0; i < cap; ++i)
            if (state[(size_t)i] == OCCUPIED) pairs.push_back({keys[(size_t)i], vals[(size_t)i]});
        int nc = (int)nextPrime(newCap);
        cap = (nc < 2) ? 2 : nc;
        keys.assign((size_t)cap, string(""));
        vals.assign((size_t)cap, 0);
        state.assign((size_t)cap, EMPTY);
        size = 0;
        used = 0;
        for (auto& kv : pairs) put(kv.first, kv.second);
    }

    void put(const string& key, long long value) {
        if ((long long)(used + 1) * 2 > cap) rehash(cap * 2);
        pair<int, bool> r = findSlot(key);
        if (r.second) {
            vals[(size_t)r.first] = value;
            return;
        }
        while (r.first < 0) {                 // 理論上不會發生，保險
            rehash(cap * 2);
            r = findSlot(key);
        }
        int idx = r.first;
        keys[(size_t)idx] = key;
        vals[(size_t)idx] = value;
        state[(size_t)idx] = OCCUPIED;
        size++;
        used++;
    }

    pair<bool, long long> get(const string& key) {
        pair<int, bool> r = findSlot(key);
        if (r.second) return {true, vals[(size_t)r.first]};
        return {false, -1};
    }

    bool erase(const string& key) {
        pair<int, bool> r = findSlot(key);
        if (!r.second) return false;
        state[(size_t)r.first] = DELETED;     // 墓碑：設回 EMPTY 會截斷探測鏈
        keys[(size_t)r.first] = "";
        vals[(size_t)r.first] = 0;
        size--;
        return true;
    }

    vector<pair<string, long long>> items() const {
        vector<pair<string, long long>> out;
        for (int i = 0; i < cap; ++i)
            if (state[(size_t)i] == OCCUPIED) out.push_back({keys[(size_t)i], vals[(size_t)i]});
        sort(out.begin(), out.end());
        return out;
    }
};

// ---------------------------------------------------------------- LRU Cache
struct LruNode {
    string key;
    long long val = 0;
    int prev = -1;
    int next = -1;
    bool sentinel = false;
};

class LRUCache {
public:
    int capacity;
    vector<LruNode> nodes;                   // 0 = head 哨兵, 1 = tail 哨兵
    unordered_map<string, int> table;        // key -> 節點下標

    explicit LRUCache(int cap = 2) : capacity(cap > 0 ? cap : 0) {
        nodes.push_back(LruNode());          // 0: head（其後是最久未使用）
        nodes.push_back(LruNode());          // 1: tail（其前是最近使用）
        nodes[0].sentinel = true;
        nodes[1].sentinel = true;
        nodes[0].next = 1;
        nodes[1].prev = 0;
    }

    void detach(int id) {
        int p = nodes[(size_t)id].prev, n = nodes[(size_t)id].next;
        nodes[(size_t)p].next = n;
        nodes[(size_t)n].prev = p;
    }

    void attachTail(int id) {
        int last = nodes[1].prev;
        nodes[(size_t)last].next = id;
        nodes[(size_t)id].prev = last;
        nodes[(size_t)id].next = 1;
        nodes[1].prev = id;
    }

    long long get(const string& key) {
        auto it = table.find(key);
        if (it == table.end()) return -1;
        int id = it->second;
        detach(id);
        attachTail(id);
        return nodes[(size_t)id].val;
    }

    void put(const string& key, long long value) {
        if (capacity == 0) return;
        auto it = table.find(key);
        if (it != table.end()) {
            int id = it->second;
            nodes[(size_t)id].val = value;
            detach(id);
            attachTail(id);
            return;
        }
        if ((int)table.size() >= capacity) {
            int victim = nodes[0].next;
            detach(victim);
            table.erase(nodes[(size_t)victim].key);
        }
        int id = (int)nodes.size();
        nodes.push_back(LruNode());
        nodes[(size_t)id].key = key;
        nodes[(size_t)id].val = value;
        table[key] = id;
        attachTail(id);
    }

    int size() const { return (int)table.size(); }

    // 由最久未使用到最近使用
    vector<pair<string, long long>> snapshot() const {
        vector<pair<string, long long>> res;
        int cur = nodes[0].next;
        while (cur != 1) {
            res.push_back({nodes[(size_t)cur].key, nodes[(size_t)cur].val});
            cur = nodes[(size_t)cur].next;
        }
        return res;
    }
};

// ---------------------------------------------------------------- 小工具
// 只接受 [+-]?digits，與 Python 版的 parse_int 完全一致；其它 token 一律當成「輸入結束」。
static bool tryLL(const string& s, long long& out) {
    if (s.empty()) return false;
    size_t i = 0;
    bool neg = false;
    if (s[0] == '-' || s[0] == '+') {
        neg = (s[0] == '-');
        i = 1;
    }
    if (i >= s.size()) return false;
    long long val = 0;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') return false;
        val = val * 10 + (s[i] - '0');
        if (val > 4000000000000000000LL) return false;      // 溢出保護
    }
    out = neg ? -val : val;
    return true;
}

static string joinItems(const vector<pair<string, long long>>& v) {
    ostringstream os;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ' ';
        os << v[i].first << ':' << v[i].second;
    }
    return os.str();
}

static mt19937 rngEngine;
static int rndInt(int lo, int hi) { return lo + (int)(rngEngine() % (unsigned)(hi - lo + 1)); }

// ---------------------------------------------------------------- IO 模式
static void runIo(const vector<string>& toks) {
    size_t pos = 0;
    auto nxt = [&]() -> const string* {
        if (pos < toks.size()) return &toks[pos++];
        return nullptr;
    };
    auto nextInt = [&]() -> optional<long long> {
        const string* v = nxt();
        if (!v) return nullopt;
        long long out = 0;
        if (!tryLL(*v, out)) return nullopt;
        return out;
    };
    auto nextIntOr = [&](int defVal) -> int {
        optional<long long> v = nextInt();
        return v ? (int)*v : defVal;
    };

    const string* mt = nxt();
    string mode = mt ? *mt : string("linear");
    if (mode != "linear" && mode != "quadratic" && mode != "double") mode = "linear";
    int cap = nextIntOr(8);
    int m = nextIntOr(0);

    OpenAddressingHashTable ht(cap, mode);
    vector<string> getLines;
    bool stopped = false;
    for (int i = 0; i < m; ++i) {
        const string* op = nxt();
        if (!op) { stopped = true; break; }
        if (*op == "put") {
            const string* k = nxt();
            optional<long long> v = nextInt();
            if (!k || !v) { stopped = true; break; }
            ht.put(*k, *v);
        } else if (*op == "get") {
            const string* k = nxt();
            if (!k) { stopped = true; break; }
            pair<bool, long long> r = ht.get(*k);
            getLines.push_back(r.first ? ("1 " + to_string(r.second)) : string("0 -1"));
        } else if (*op == "del") {
            const string* k = nxt();
            if (!k) { stopped = true; break; }
            ht.erase(*k);
        }
    }

    int lruCap = stopped ? 0 : nextIntOr(0);
    int q = stopped ? 0 : nextIntOr(0);

    LRUCache cache(lruCap);
    vector<string> lruLines;
    for (int i = 0; i < q; ++i) {
        const string* op = nxt();
        if (!op) break;
        if (*op == "put") {
            const string* k = nxt();
            optional<long long> v = nextInt();
            if (!k || !v) break;
            cache.put(*k, *v);
        } else if (*op == "get") {
            const string* k = nxt();
            if (!k) break;
            lruLines.push_back(to_string(cache.get(*k)));
        }
    }

    cout << ht.cap << '\n';
    cout << ht.size << '\n';
    cout << ht.probes << '\n';
    cout << joinItems(ht.items()) << '\n';
    for (const string& line : getLines) cout << line << '\n';
    cout << cache.size() << '\n';
    cout << joinItems(cache.snapshot()) << '\n';
    for (const string& line : lruLines) cout << line << '\n';
}

// ---------------------------------------------------------------- 測試
static void runTests() {
    // 質數工具
    {
        int in[] = {0, 1, 2, 3, 4, 10, 11, 12};
        int want[] = {2, 2, 2, 3, 5, 11, 11, 13};
        for (int i = 0; i < 8; ++i) assert(nextPrime(in[i]) == want[i]);
    }

    // 墓碑不能截斷探測鏈
    for (const string& mode : {string("linear"), string("quadratic"), string("double")}) {
        OpenAddressingHashTable ht(5, mode);
        const char* ks[] = {"aa", "bb", "cc", "dd"};
        for (int i = 0; i < 4; ++i) ht.put(ks[i], i + 1);
        assert(ht.size == 4);
        assert(ht.get("cc") == make_pair(true, 3LL));
        assert(ht.erase("bb"));
        assert(ht.size == 3);
        assert(ht.get("bb") == make_pair(false, -1LL));
        assert(ht.get("cc") == make_pair(true, 3LL));   // 關鍵：刪掉 bb 後 cc 仍查得到
        assert(ht.get("dd") == make_pair(true, 4LL));
        assert(!ht.erase("bb"));                        // 重複刪除
        ht.put("bb", 22);                               // 重新插入複用墓碑
        assert(ht.get("bb") == make_pair(true, 22LL));
        assert(ht.size == 4);
    }

    // 更新既有 key 不增加 size
    {
        OpenAddressingHashTable ht(8, "linear");
        ht.put("x", 1);
        ht.put("x", 2);
        assert(ht.size == 1);
        assert(ht.get("x") == make_pair(true, 2LL));
    }

    // 三種策略的大量插入 / 刪除 / 再插入
    {
        vector<string> keys;
        for (int i = 0; i < 60; ++i) {
            ostringstream os;
            os << "k" << (i / 100) << (i / 10) % 10 << i % 10;
            keys.push_back(os.str());
        }
        for (const string& mode : {string("linear"), string("quadratic"), string("double")}) {
            OpenAddressingHashTable ht(4, mode);
            for (int i = 0; i < (int)keys.size(); ++i) ht.put(keys[(size_t)i], i);
            assert(ht.size == 60);
            assert((long long)ht.used * 2 <= ht.cap);    // 負載因子 <= 0.5
            for (int i = 0; i < (int)keys.size(); ++i)
                assert(ht.get(keys[(size_t)i]) == make_pair(true, (long long)i));
            for (int i = 0; i < 30; ++i) assert(ht.erase(keys[(size_t)i]));
            assert(ht.size == 30);
            for (int i = 0; i < (int)keys.size(); ++i) {
                if (i < 30) assert(ht.get(keys[(size_t)i]) == make_pair(false, -1LL));
                else assert(ht.get(keys[(size_t)i]) == make_pair(true, (long long)i));
            }
        }
    }

    // 與 map 隨機對拍
    rngEngine.seed(20261005);
    for (const string& mode : {string("linear"), string("quadratic"), string("double")}) {
        for (int t = 0; t < 30; ++t) {
            int cap = rndInt(1, 12);
            OpenAddressingHashTable ht(cap, mode);
            map<string, long long> ref;
            for (int step = 0; step < 300; ++step) {
                string key = "s" + to_string(rndInt(0, 40));
                int op = rndInt(0, 2);
                if (op == 0) {
                    long long val = rndInt(-50, 50);
                    ht.put(key, val);
                    ref[key] = val;
                } else if (op == 1) {
                    pair<bool, long long> r = ht.get(key);
                    auto it = ref.find(key);
                    if (it != ref.end()) assert(r.first && r.second == it->second);
                    else assert(r == make_pair(false, -1LL));
                } else {
                    bool gone = ht.erase(key);
                    auto it = ref.find(key);
                    assert(gone == (it != ref.end()));
                    if (it != ref.end()) ref.erase(it);
                }
                assert(ht.size == (int)ref.size());
                assert((long long)ht.used * 2 <= ht.cap);
                vector<pair<string, long long>> want(ref.begin(), ref.end());
                assert(ht.items() == want);
            }
        }
    }

    // LRU：基本行為
    {
        LRUCache cache(2);
        assert(cache.get("a") == -1);
        cache.put("a", 1);
        cache.put("b", 2);
        assert(cache.get("a") == 1);
        cache.put("c", 3);                              // 淘汰最久未使用的 b
        assert(cache.get("b") == -1);
        assert(cache.get("a") == 1);
        assert(cache.get("c") == 3);
        assert(cache.size() == 2);
        assert((cache.snapshot() == vector<pair<string, long long>>{{"a", 1}, {"c", 3}}));
    }

    // LRU：get 刷新順序
    {
        LRUCache cache(3);
        cache.put("a", 1);
        cache.put("b", 2);
        cache.put("c", 3);
        assert(cache.get("a") == 1);
        cache.put("d", 4);                              // 淘汰 b
        assert(cache.get("b") == -1);
        assert((cache.snapshot() == vector<pair<string, long long>>{{"c", 3}, {"a", 1}, {"d", 4}}));
    }

    // LRU：put 既有 key 只更新值、不淘汰
    {
        LRUCache cache(2);
        cache.put("a", 1);
        cache.put("b", 2);
        cache.put("a", 10);
        assert(cache.size() == 2);
        assert(cache.get("a") == 10);
        assert((cache.snapshot() == vector<pair<string, long long>>{{"b", 2}, {"a", 10}}));
    }

    // LRU：容量 0
    {
        LRUCache cache(0);
        cache.put("a", 1);
        assert(cache.size() == 0);
        assert(cache.get("a") == -1);
        assert(cache.snapshot().empty());
    }

    // LRU：容量 1
    {
        LRUCache cache(1);
        cache.put("a", 1);
        cache.put("b", 2);
        assert(cache.size() == 1);
        assert(cache.get("a") == -1);
        assert(cache.get("b") == 2);
    }

    // LRU 隨機對拍：list 當樸素基準
    for (int t = 0; t < 30; ++t) {
        int cap = rndInt(1, 8);
        LRUCache cache(cap);
        vector<string> order;                           // 最久未使用 -> 最近使用
        map<string, long long> vals;
        for (int step = 0; step < 200; ++step) {
            string key = "x" + to_string(rndInt(0, 12));
            if (rndInt(0, 1) == 0) {
                long long val = rndInt(0, 99);
                cache.put(key, val);
                if (cap == 0) continue;
                auto it = vals.find(key);
                if (it != vals.end()) {
                    order.erase(find(order.begin(), order.end(), key));
                    it->second = val;
                } else {
                    if ((int)order.size() >= cap) {
                        vals.erase(order.front());
                        order.erase(order.begin());
                    }
                    vals[key] = val;
                }
                order.push_back(key);
            } else {
                long long got = cache.get(key);
                auto it = vals.find(key);
                if (it != vals.end()) {
                    assert(got == it->second);
                    order.erase(find(order.begin(), order.end(), key));
                    order.push_back(key);
                } else {
                    assert(got == -1);
                }
            }
            assert(cache.size() == (int)vals.size());
            vector<pair<string, long long>> want;
            for (const string& k : order) want.push_back({k, vals[k]});
            assert(cache.snapshot() == want);
        }
    }

    cout << "all tests passed" << '\n';
}

int main() {
    vector<string> toks;
    string t;
    while (cin >> t) toks.push_back(t);
    if (toks.empty()) runTests();
    else runIo(toks);
    return 0;
}
