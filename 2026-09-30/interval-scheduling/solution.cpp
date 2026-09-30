// 區間問題與貪心：最多不重疊區間 / 最少移除 / 最少會議室 / 合併區間 / 引爆氣球
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution && ./solution
//
// 思路：
//   區間題 90% 的第一步都是排序，排完之後關係變成單調，一次線性掃描就能解決。
//   1. 最多不重疊區間：按**右端點升序**貪心，每次選結束最早的那個（交換論證：
//      用結束最早的替換最優解的第一個，不會與後續衝突，規模不變）。
//   2. 最少移除 = n − 最多保留。
//   3. 最少會議室：按左端點掃描，小根堆存「正在使用的房間的結束時間」，
//      先把所有 end <= l 的房間釋放，再壓入新的，堆的最大容量即答案。
//      也可以用掃描線（左端 +1、右端 −1，同座標先處理 −1）得到同樣結果。
//   4. 合併區間：按左端點掃描，l <= 當前右端就合併（端點相接也要合併）。
//   5. 引爆氣球：按右端點升序，箭射在右端點 x，只有 l > x 時才補新箭（閉區間）。
//
//   「端點相接是否算重疊」是最容易踩的坑，本題統一約定：
//     - 不重疊 / 會議室：半開語義 [l, r)，端點相接**不算**衝突，判定用 >= / <=；
//     - 合併區間 / 引爆氣球：閉區間 [l, r]，端點相接**算**重疊，判定用 <= / >。
//   行爲與 Python 版逐字節一致。
//
// 輸入（空白分隔）：n / l1 r1 / l2 r2 / …
// 輸出：最多不重疊數 / 最少移除數 / 最少會議室數 / 合併後區間數 / 合併結果 / 最少箭數

#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using PII = pair<int, int>;

// ---------- 各問解法 ----------

// 1. 最多互不重疊區間數：按右端點升序貪心（端點相接不算衝突）
static int maxNonOverlap(vector<PII> v) {
    sort(v.begin(), v.end(), [](const PII& a, const PII& b) {
        if (a.second != b.second) return a.second < b.second;
        return a.first < b.first;
    });
    int cnt = 0;
    bool has = false;
    int lastEnd = 0;
    for (auto& e : v) {
        if (!has || e.first >= lastEnd) {
            cnt++;
            lastEnd = e.second;
            has = true;
        }
    }
    return cnt;
}

// 同上，但用 O(n^2) 動態規劃暴力求解，測試裡當貪心的對拍基準
static int maxNonOverlapDp(vector<PII> v) {
    sort(v.begin(), v.end(), [](const PII& a, const PII& b) {
        if (a.second != b.second) return a.second < b.second;
        return a.first < b.first;
    });
    int n = (int)v.size();
    vector<int> dp(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        int l = v[i - 1].first;
        int p = 0;
        for (int j = 0; j < i - 1; ++j)
            if (v[j].second <= l) p++;   // 與第 i 個不衝突的前綴個數
        dp[i] = max(dp[i - 1], 1 + dp[p]);
    }
    return dp[n];
}

// 2. 最少移除 = n − 最多保留
static int minRemoved(const vector<PII>& v) {
    return (int)v.size() - maxNonOverlap(v);
}

// 3a. 最少會議室數：小根堆
static int minMeetingRoomsHeap(vector<PII> v) {
    sort(v.begin(), v.end());
    priority_queue<int, vector<int>, greater<int>> pq;   // 最早空出來的房間
    int busy = 0;
    for (auto& e : v) {
        int l = e.first, r = e.second;
        if (l == r) continue;                 // [l, l) 是空區間，不佔用房間
        while (!pq.empty() && pq.top() <= l) pq.pop();
        pq.push(r);
        busy = max(busy, (int)pq.size());
    }
    return busy;
}

// 3b. 最少會議室數：掃描線（同座標先處理結束）
static int minMeetingRoomsSweep(const vector<PII>& v) {
    vector<PII> ev;
    for (auto& e : v) {
        if (e.first == e.second) continue;
        ev.push_back({e.first, 1});
        ev.push_back({e.second, -1});
    }
    sort(ev.begin(), ev.end());               // -1 < 1，同座標先退房再入住
    int cur = 0, best = 0;
    for (auto& e : ev) {
        cur += e.second;
        best = max(best, cur);
    }
    return best;
}

// 4. 合併區間
static vector<PII> mergeIntervals(vector<PII> v) {
    sort(v.begin(), v.end());
    vector<PII> res;
    for (auto& e : v) {
        if (!res.empty() && e.first <= res.back().second)
            res.back().second = max(res.back().second, e.second);
        else
            res.push_back(e);
    }
    return res;
}

// 5. 引爆氣球的最少箭數（閉區間）
static int minArrows(vector<PII> v) {
    sort(v.begin(), v.end(), [](const PII& a, const PII& b) {
        if (a.second != b.second) return a.second < b.second;
        return a.first < b.first;
    });
    int arrows = 0;
    bool has = false;
    int lastX = 0;
    for (auto& e : v) {
        if (!has || e.first > lastX) {
            arrows++;
            lastX = e.second;
            has = true;
        }
    }
    return arrows;
}

// 同上，枚舉候選位置的子集求最優解（座標範圍必須很小），只用於測試對拍
static int minArrowsBrute(const vector<PII>& v) {
    if (v.empty()) return 0;
    vector<int> cand;
    for (auto& e : v) {
        cand.push_back(e.first);
        cand.push_back(e.second);
    }
    sort(cand.begin(), cand.end());
    cand.erase(unique(cand.begin(), cand.end()), cand.end());
    int c = (int)cand.size(), n = (int)v.size();
    int full = (1 << n) - 1;
    // hit[i] = 箭射在 cand[i] 處能打掉的氣球集合（位元遮罩）
    vector<int> hit(c, 0);
    for (int i = 0; i < c; ++i)
        for (int j = 0; j < n; ++j)
            if (v[j].first <= cand[i] && cand[i] <= v[j].second) hit[i] |= (1 << j);
    int best = c;
    for (int mask = 1; mask < (1 << c); ++mask) {
        int cover = 0, bits = 0;
        for (int i = 0; i < c; ++i)
            if (mask & (1 << i)) {
                cover |= hit[i];
                bits++;
            }
        if (cover == full) best = min(best, bits);
    }
    return best;
}

// ---------- 輸入輸出 ----------

static void runIo(const string& data) {
    istringstream iss(data);
    vector<int> tk;
    int x;
    while (iss >> x) tk.push_back(x);
    int n = tk.empty() ? 0 : tk[0];
    vector<PII> iv;
    for (int i = 0; i < n && 1 + 2 * i + 1 < (int)tk.size(); ++i) {
        int l = tk[1 + 2 * i], r = tk[1 + 2 * i + 1];
        if (l > r) swap(l, r);
        iv.push_back({l, r});
    }
    vector<PII> mg = mergeIntervals(iv);
    ostringstream out;
    out << maxNonOverlap(iv) << '\n';
    out << minRemoved(iv) << '\n';
    out << minMeetingRoomsHeap(iv) << '\n';
    out << mg.size() << '\n';
    for (size_t i = 0; i < mg.size(); ++i) {
        if (i) out << ' ';
        out << mg[i].first << ' ' << mg[i].second;
    }
    out << '\n';
    out << minArrows(iv) << '\n';
    cout << out.str();
}

// ---------- 內置測試 ----------

static bool covered(const vector<PII>& iv, int x) {
    for (auto& e : iv)
        if (e.first <= x && x <= e.second) return true;
    return false;
}

static void runTests() {
    vector<PII> a = {{1, 2}, {2, 3}, {3, 4}, {1, 3}};
    assert(maxNonOverlap(a) == 3);
    assert(minRemoved(a) == 1);
    assert(minMeetingRoomsHeap(a) == 2);
    assert(minMeetingRoomsSweep(a) == 2);
    assert((mergeIntervals(a) == vector<PII>{{1, 4}}));
    assert(minArrows(a) == 2);

    vector<PII> b = {{1, 4}, {2, 3}, {3, 5}, {7, 9}};
    assert(maxNonOverlap(b) == 3);
    assert(minRemoved(b) == 1);
    assert(minMeetingRoomsHeap(b) == 2);
    assert((mergeIntervals(b) == vector<PII>{{1, 5}, {7, 9}}));
    assert(minArrows(b) == 2);

    vector<PII> empty;
    assert(maxNonOverlap(empty) == 0);
    assert(minRemoved(empty) == 0);
    assert(minMeetingRoomsHeap(empty) == 0);
    assert(mergeIntervals(empty).empty());
    assert(minArrows(empty) == 0);
    assert(minArrowsBrute(empty) == 0);

    vector<PII> c = {{5, 5}, {5, 5}, {6, 6}};
    assert(maxNonOverlap(c) == 3);
    assert(minMeetingRoomsHeap(c) == 0);      // 半開語義下空區間不佔用會議室
    assert((mergeIntervals(c) == vector<PII>{{5, 5}, {6, 6}}));
    assert(minArrows(c) == 2);

    mt19937 rng(20260930);
    for (int t = 0; t < 600; ++t) {
        int n = (int)(rng() % 9);
        vector<PII> iv;
        for (int i = 0; i < n; ++i) {
            int l = (int)(rng() % 9);
            int r = l + (int)(rng() % 5);
            iv.push_back({l, r});
        }

        assert(maxNonOverlap(iv) == maxNonOverlapDp(iv));       // 貪心 vs DP
        assert(minRemoved(iv) == n - maxNonOverlap(iv));
        assert(minMeetingRoomsHeap(iv) == minMeetingRoomsSweep(iv));

        vector<PII> mg = mergeIntervals(iv);
        assert(is_sorted(mg.begin(), mg.end()));
        for (size_t i = 1; i < mg.size(); ++i) assert(mg[i - 1].second < mg[i].first);
        for (int x = -1; x <= 14; ++x)
            assert(covered(mg, x) == covered(iv, x));           // 覆蓋點集一致

        int ans = minArrows(iv);
        assert(ans <= maxNonOverlap(iv));
        // 貪心構造出的箭位置確實能打掉全部氣球
        vector<PII> sv = iv;
        sort(sv.begin(), sv.end(), [](const PII& p, const PII& q) {
            if (p.second != q.second) return p.second < q.second;
            return p.first < q.first;
        });
        vector<int> xs;
        bool has = false;
        int last = 0;
        for (auto& e : sv) {
            if (!has || e.first > last) {
                xs.push_back(e.second);
                last = e.second;
                has = true;
            }
        }
        assert((int)xs.size() == ans);
        for (auto& e : iv) {
            bool ok = false;
            for (int x : xs)
                if (e.first <= x && x <= e.second) ok = true;
            assert(ok);
        }
    }

    // 與子集枚舉的最優解對拍（座標範圍取小一些，控制枚舉量）
    for (int t = 0; t < 200; ++t) {
        int n = (int)(rng() % 8);
        vector<PII> iv;
        for (int i = 0; i < n; ++i) {
            int l = (int)(rng() % 7);
            int r = l + (int)(rng() % 4);
            iv.push_back({l, r});
        }
        assert(minArrowsBrute(iv) == minArrows(iv));
    }

    cout << "all tests passed" << '\n';
}

int main() {
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
