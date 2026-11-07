// 同餘最短路（Congruence Shortest Path）
//
// 給定正整數 b_1..b_n 與多個查詢 x，求最小的 y >= x 使得
// y = a_1*b_1 + ... + a_n*b_n（a_i >= 0）。不存在則輸出 -1。
//
// 記 m = min(b_i)，節點為剩餘類 0..m-1，邊 r -> (r+b_i) mod m，權 +b_i，
// 從 0 跑 Dijkstra 得 dist[r] = 可湊出且 ≡ r (mod m) 的最小值。
// m 本身是某個 b_i，故類內可湊值恰為 {dist[r] + t*m}；查詢 y >= x 不要求
// 與 x 同餘，答案對全部剩餘類取 min（見 solution.py 的推導）。
// 複雜度 O(n*m*log m + q*m)。
//
// IO 格式與 solution.py 完全一致（見同目錄 README.md）。
#include <cstdio>
#include <cstring>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

typedef long long ll;
static const ll INF = 1LL << 62;

static std::vector<ll> dijkstraResidues(const std::vector<int>& bs, int m) {
    std::vector<ll> dist(m, INF);
    dist[0] = 0;
    // 小根堆：（距離, 剩餘類）
    std::priority_queue<std::pair<ll, int>,
                        std::vector<std::pair<ll, int>>,
                        std::greater<std::pair<ll, int>>> pq;
    pq.push({0, 0});
    while (!pq.empty()) {
        std::pair<ll, int> top = pq.top();
        pq.pop();
        ll d = top.first;
        int r = top.second;
        if (d > dist[r]) continue;
        for (size_t i = 0; i < bs.size(); ++i) {
            int b = bs[i];
            int nr = (int)(((ll)r + b) % m);
            ll nd = d + b;
            if (nd < dist[nr]) {
                dist[nr] = nd;
                pq.push({nd, nr});
            }
        }
    }
    return dist;
}

static ll solveQuery(ll x, const std::vector<ll>& dist, int m) {
    ll best = INF;
    for (int r = 0; r < m; ++r) {
        ll f = dist[r];
        if (f >= INF) continue;
        ll c;
        if (f >= x) {
            c = f;
        } else {
            ll t = (x - f + m - 1) / m;  // ceil((x - f) / m) >= 1
            c = f + t * m;
        }
        if (c < best) best = c;
    }
    return best >= INF ? -1 : best;
}

// ---------- 內建測試 ----------
struct RNG {
    unsigned int s;
    explicit RNG(unsigned int seed) : s(seed) {}
    int next(int m) {
        s = (s * 1103515245u + 12345u) & 0x7FFFFFFFu;
        return (int)(s % (unsigned int)m);
    }
};

static std::vector<ll> brute(const std::vector<int>& bs, int m, int lim) {
    std::vector<char> reach(lim + 1, 0);
    reach[0] = 1;
    for (int y = 1; y <= lim; ++y)
        for (size_t i = 0; i < bs.size(); ++i)
            if (y >= bs[i] && reach[y - bs[i]]) { reach[y] = 1; break; }
    std::vector<ll> bf(m, INF);
    for (int y = 0; y <= lim; ++y)
        if (reach[y]) {
            int r = y % m;
            if (y < bf[r]) bf[r] = y;
        }
    return bf;
}

static void runTests() {
    RNG rng(20261106u);
    std::vector<std::vector<int>> fixed = {
        {3, 5}, {4, 6}, {7}, {6, 10, 15}, {5, 7, 11},
    };
    for (size_t c = 0; c < fixed.size(); ++c) {
        std::vector<int> bs = fixed[c];
        int m = bs[0];
        std::vector<ll> dist = dijkstraResidues(bs, m);
        int lim = m * bs[bs.size() - 1];
        std::vector<ll> bf = brute(bs, m, lim);
        for (int r = 0; r < m; ++r)
            if (dist[r] != bf[r]) { printf("FAIL dist\n"); exit(1); }
        for (int x = 0; x <= lim; ++x) {
            ll want = -1;  // 暴力查詢
            for (int y = x; y <= lim; ++y) {
                ll fr = bf[y % m];
                if (fr < INF && y >= fr) { want = y; break; }
            }
            ll got = solveQuery(x, dist, m);
            if (got != want) { printf("FAIL query bs=%d x=%d want=%lld got=%lld\n",
                                      (int)c, x, want, got); exit(1); }
        }
    }
    for (int t = 0; t < 200; ++t) {
        int n = rng.next(5) + 1;
        std::vector<int> bs;
        for (int i = 0; i < n; ++i) bs.push_back(rng.next(60) + 1);
        // 與 Python 版相同：直接排序，m = bs[0]
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (bs[j] < bs[i]) std::swap(bs[i], bs[j]);
        int m = bs[0];
        std::vector<ll> dist = dijkstraResidues(bs, m);
        int lim = m * bs[bs.size() - 1];
        std::vector<ll> bf = brute(bs, m, lim);
        for (int r = 0; r < m; ++r)
            if (dist[r] != bf[r]) { printf("FAIL dist rand\n"); exit(1); }
        for (int k = 0; k < 10; ++k) {
            int x = rng.next(lim + 1);
            ll want = -1;
            for (int y = x; y <= lim; ++y) {
                ll fr = bf[y % m];
                if (fr < INF && y >= fr) { want = y; break; }
            }
            ll got = solveQuery(x, dist, m);
            if (got != want) { printf("FAIL query rand x=%d\n", x); exit(1); }
        }
    }
    printf("all tests passed\n");
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    long long first;
    if (scanf("%lld", &first) != 1) {
        runTests();
        return 0;
    }
    int n = (int)first;
    std::vector<int> bs(n);
    for (int i = 0; i < n; ++i)
        if (scanf("%d", &bs[i]) != 1) { runTests(); return 0; }
    long long q;
    if (scanf("%lld", &q) != 1) { runTests(); return 0; }
    int m = bs[0];
    for (int i = 1; i < n; ++i)
        if (bs[i] < m) m = bs[i];
    std::vector<ll> dist = dijkstraResidues(bs, m);
    std::string out;
    out.reserve((size_t)q * 8);
    char buf[32];
    for (long long t = 0; t < q; ++t) {
        long long x;
        if (scanf("%lld", &x) != 1) break;
        sprintf(buf, "%lld\n", solveQuery(x, dist, m));
        out += buf;
    }
    fputs(out.c_str(), stdout);
    return 0;
}
