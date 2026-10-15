// Z 函數（Z Algorithm）與擴展 KMP
//
// 題意：
//     給定字串 s（文本）與 p（模式串），求：
//       1. s 的 Z 陣列：z[i] = s 與 s[i:] 的最長公共前綴長度（定義 z[0] = |s|）
//       2. p 的 Z 陣列
//       3. p 在 s 中所有出現的起始位置（0-indexed）
//       4. 擴展 KMP 的 ext 陣列：ext[i] = p 與 s[i:] 的最長公共前綴長度
//       5. s 的最短週期（最小的 k 使 s 由某個長度 k 的塊重複整數次構成）
//     以上全部同時用「O(n) 演算法」與「O(n²) 暴力」計算並互相驗證。
//
// 思路：
//     ### Z 函數：用已經算過的區間當跳板
//     維護一個「Z-box」[l, r]，它是當前已知的、最靠右的匹配區間，
//     也就是 s[l..r] == s[0..r-l]。算 z[i] 時：
//       - 若 i > r：Z-box 蓋不到 i，只能從 0 開始老實比（while 迴圈）。
//       - 若 i ≤ r：因為 s[i..r] 已經等於 s[0..r-i]，
//         可以直接複製 z[i-l]（同一個相對位置先前算過的值），
//         但最多只能複製到 r 為止 → 取 min(r - i + 1, z[i - l])，再繼續往後比。
//     關鍵是每個字元最多被 while 迴圈比對一次：一旦匹配超過 r，
//     就把 Z-box 右端推進到新位置，而 r 只增不減。故總複雜度 O(n)。
//     這跟 KMP 的 prefix function 是對偶的兩套方法：KMP 預處理模式串、
//     掃描時不回頭；Z 函數則把「自己和自己的每個後綴」一次算完。
//
//     ### 擴展 KMP：一次 Z 搞定所有後綴
//     要求 p 與 s 的每個後綴 s[i:] 的 LCP，只要在
//     p + '#' + s（'#' 是不出現在字串裡的哨兵）上跑一次 Z 函數：
//     位置 |p|+1+i 對應 s[i:]，取其 z 值再跟 |p| 取 min 即可（哨兵會擋住越界）。
//     這比「對每個 i 都跑一次 KMP」快得多，是 O(|p| + |s|)。
//
//     ### 最短週期
//     若 s 有週期 k（且 k | n），則 s[k:] == s[:n-k]，也就是 z[k] == n - k，
//     反過來也成立。所以枚舉 n 的因數 k，取第一個滿足 z[k] == n-k 的就是答案。
//     這是 Z 函數在「字串週期性」問題上的標準應用（冪次判定、字串壓縮）。
//
// 複雜度：
//     Z 陣列                    O(n) 時間 / O(n) 空間
//     擴展 KMP（所有後綴 LCP）  O(|p| + |s|) 時間 / O(|p| + |s|) 空間
//     模式匹配（所有出現位置）  O(|p| + |s| + 出現次數)
//     最短週期                  O(n)（枚舉因數 + Z 陣列）

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

const int MAXL = 2000;      // 單串最大長度（IO 模式下的截斷上限）
const char SENTINEL = '#';  // 哨兵字元：不在輸入中出現

// ------------------------------------------------------- 核心演算法

std::vector<int> zAlgorithm(const std::string& s) {
    int n = static_cast<int>(s.size());
    std::vector<int> z(n, 0);
    if (n == 0) return z;
    z[0] = n;
    int l = 0, r = 0;  // 當前最右的 Z-box [l, r]（閉區間）
    for (int i = 1; i < n; ++i) {
        if (i <= r) {
            // s[i..r] 已知等於 s[0..r-i]，直接複製舊答案，但不要超出 r
            z[i] = std::min(r - i + 1, z[i - l]);
        }
        while (i + z[i] < n && s[z[i]] == s[i + z[i]]) {
            ++z[i];
        }
        if (i + z[i] - 1 > r) {  // 匹配突破了舊的右端，推進 Z-box
            l = i;
            r = i + z[i] - 1;
        }
    }
    return z;
}

std::vector<int> extendedKmp(const std::string& text, const std::string& pat) {
    int n = static_cast<int>(text.size());
    int m = static_cast<int>(pat.size());
    std::vector<int> ext(n, 0);
    if (n == 0 || m == 0) return ext;
    std::vector<int> z = zAlgorithm(pat + SENTINEL + text);
    for (int i = 0; i < n; ++i) {
        ext[i] = std::min(z[m + 1 + i], m);
    }
    return ext;
}

std::vector<int> zMatch(const std::string& text, const std::string& pat) {
    int n = static_cast<int>(text.size());
    int m = static_cast<int>(pat.size());
    std::vector<int> occ;
    if (m == 0) {
        for (int i = 0; i <= n; ++i) occ.push_back(i);
        return occ;
    }
    std::vector<int> ext = extendedKmp(text, pat);
    for (int i = 0; i <= n - m; ++i) {
        if (ext[i] == m) occ.push_back(i);
    }
    return occ;
}

int smallestPeriod(const std::string& s) {
    int n = static_cast<int>(s.size());
    if (n == 0) return 0;
    std::vector<int> z = zAlgorithm(s);
    for (int k = 1; k < n; ++k) {  // k = n 必然是週期，留作兜底（z[n] 不存在）
        if (n % k == 0 && z[k] == n - k) return k;
    }
    return n;
}

// ------------------------------------------------------- 暴力對拍

std::vector<int> zBrute(const std::string& s) {
    int n = static_cast<int>(s.size());
    std::vector<int> res(n, 0);
    for (int i = 0; i < n; ++i) {
        int k = 0;
        while (i + k < n && s[k] == s[i + k]) ++k;
        res[i] = k;
    }
    if (n) res[0] = n;
    return res;
}

std::vector<int> extBrute(const std::string& text, const std::string& pat) {
    int n = static_cast<int>(text.size());
    int m = static_cast<int>(pat.size());
    std::vector<int> res(n, 0);
    for (int i = 0; i < n; ++i) {
        int k = 0;
        while (k < m && i + k < n && pat[k] == text[i + k]) ++k;
        res[i] = k;
    }
    return res;
}

std::vector<int> matchBrute(const std::string& text, const std::string& pat) {
    int n = static_cast<int>(text.size());
    int m = static_cast<int>(pat.size());
    std::vector<int> occ;
    if (m == 0) {
        for (int i = 0; i <= n; ++i) occ.push_back(i);
        return occ;
    }
    for (int i = 0; i <= n - m; ++i) {
        if (text.compare(i, m, pat) == 0) occ.push_back(i);
    }
    return occ;
}

int periodBrute(const std::string& s) {
    int n = static_cast<int>(s.size());
    for (int k = 1; k <= n; ++k) {
        if (n % k) continue;
        bool ok = true;
        for (int i = 0; i < n && ok; ++i) {
            if (s[i] != s[i % k]) ok = false;
        }
        if (ok) return k;
    }
    return n;
}

// ------------------------------------------------------- 輸出

std::string join(const std::vector<int>& v) {
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}

std::string run(const std::string& text, const std::string& pat) {
    std::vector<int> zs = zAlgorithm(text);
    std::vector<int> zp = zAlgorithm(pat);
    std::vector<int> occ = zMatch(text, pat);
    std::vector<int> ext = extendedKmp(text, pat);
    int per = smallestPeriod(text);

    int ok = 1;
    if (zs != zBrute(text)) ok = 0;
    if (zp != zBrute(pat)) ok = 0;
    if (occ != matchBrute(text, pat)) ok = 0;
    if (ext != extBrute(text, pat)) ok = 0;
    if (per != periodBrute(text)) ok = 0;
    // 出現位置必須與 ext 陣列自洽：ext[i] == |p| ⟺ i 是合法起點且完全匹配
    int m = static_cast<int>(pat.size());
    if (m > 0) {
        std::vector<int> derived;
        for (int i = 0; i <= static_cast<int>(text.size()) - m; ++i) {
            if (ext[i] == m) derived.push_back(i);
        }
        if (derived != occ) ok = 0;
    }

    std::ostringstream oss;
    oss << join(zs) << '\n'
        << join(zp) << '\n'
        << join(occ) << '\n'
        << join(ext) << '\n'
        << per << '\n'
        << ok;
    return oss.str();
}

// ------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例 ----
    assert((zAlgorithm("aaaaa") == std::vector<int>{5, 4, 3, 2, 1}));
    assert((zAlgorithm("ababab") == std::vector<int>{6, 0, 4, 0, 2, 0}));
    assert(zAlgorithm("").empty());
    assert((zAlgorithm("a") == std::vector<int>{1}));
    assert((zAlgorithm("abcababc") == std::vector<int>{8, 0, 0, 2, 0, 3, 0, 0}));

    assert((zMatch("abababab", "abab") == std::vector<int>{0, 2, 4}));
    assert((zMatch("aaaaa", "aa") == std::vector<int>{0, 1, 2, 3}));
    assert(zMatch("abc", "xyz").empty());
    assert((zMatch("abc", "") == std::vector<int>{0, 1, 2, 3}));
    assert(zMatch("", "a").empty());
    assert((zMatch("", "") == std::vector<int>{0}));

    assert((extendedKmp("aaaaa", "aa") == std::vector<int>{2, 2, 2, 2, 1}));
    assert((extendedKmp("abcabc", "abc") == std::vector<int>{3, 0, 0, 3, 0, 0}));
    assert((extendedKmp("abc", "") == std::vector<int>{0, 0, 0}));

    assert(smallestPeriod("ababab") == 2);
    assert(smallestPeriod("aaaa") == 1);
    assert(smallestPeriod("abcab") == 5);
    assert(smallestPeriod("abcabcabc") == 3);
    assert(smallestPeriod("") == 0);

    // ---- 隨機對拍 ----
    std::mt19937 rng(20261015);
    const char* alphabets[4] = {"ab", "abc", "a", "abcdefgh"};
    for (int iter = 0; iter < 600; ++iter) {
        std::string alpha = alphabets[rng() % 4];
        int tl = static_cast<int>(rng() % 61);
        int pl = static_cast<int>(rng() % 13);
        std::string text, pat;
        for (int i = 0; i < tl; ++i) text.push_back(alpha[rng() % alpha.size()]);
        for (int i = 0; i < pl; ++i) pat.push_back(alpha[rng() % alpha.size()]);
        assert(zAlgorithm(text) == zBrute(text));
        assert(zAlgorithm(pat) == zBrute(pat));
        assert(extendedKmp(text, pat) == extBrute(text, pat));
        assert(zMatch(text, pat) == matchBrute(text, pat));
        assert(smallestPeriod(text) == periodBrute(text));
        std::string out = run(text, pat);
        assert(out.size() >= 2 && out[out.size() - 1] == '1');
    }

    // ---- 長串只跑線性演算法（不跑 O(n²) 暴力）----
    std::string longText;
    for (int i = 0; i < 800; ++i) longText += "ab";
    assert(smallestPeriod(longText) == 2);
    std::vector<int> expect;
    for (int i = 0; i <= static_cast<int>(longText.size()) - 4; i += 2) {
        expect.push_back(i);
    }
    assert(zMatch(longText, "abab") == expect);
    assert(zAlgorithm(std::string(1500, 'a'))[1] == 1499);
}

// ------------------------------------------------------- IO

bool isWord(const std::string& t) {
    if (t.empty()) return false;
    for (char c : t) {
        if (c < 'a' || c > 'z') return false;
    }
    return true;
}

std::string clean(const std::string& t) {
    if (!isWord(t)) return "";
    return t.substr(0, MAXL);
}

void runIo(const std::string& raw) {
    std::istringstream iss(raw);
    std::string t1, t2;
    iss >> t1 >> t2;
    std::string text = clean(t1), pat = clean(t2);
    std::cout << run(text, pat) << "\n";
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (char c : raw) {
        unsigned char u = static_cast<unsigned char>(c);
        if (u > ' ') { blank = false; break; }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
