/**
 * Title: String hashing
 * Description: Polynomial hash of `s[l, r)` in $O(1)$:
 *   $sum_(i = l)^(r - 1) s_i B^(r - 1 - i) mod M$ with the prime
 *   $M = 2^61 - 1$ and a random base $B$ shared by all instances,
 *   so substrings of different strings compare too. Different
 *   strings of length $n$ collide with probability at most
 *   $n / M approx n dot 4.3 dot 10^(-19)$.
 * Precondition: Values in $[0, M)$; cast chars above 127 to
 *   `unsigned char`. Compare hashes of equal lengths only.
 * Time: $O(n)$ to build, $O(1)$ per `get`.
 * Usage:
 *   StrHash hs(s); // string or vector<int>
 *   bool eq = hs.get(a, a + len) == hs.get(b, b + len);
 * Tweak: Hash of `x` followed by `y`, where `y` has length
 *   `k <= n`: `hs.norm(hs.mul(x, hs.pw[k]) + y)`. Keys of mixed
 *   lengths: `pair(hash, len)`.
 * Example:
 *   // Longest common prefix of s[a, n) and s[b, n):
 *   int lo = 0, hi = n - max(a, b);
 *   while (lo < hi) {
 *     int mid = (lo + hi + 1) / 2;
 *     if (hs.get(a, a + mid) == hs.get(b, b + mid)) lo = mid;
 *     else hi = mid - 1;
 *   }
 *   // s[a, n) < s[b, n) iff lo < n - b && (lo == n - a
 *   //   || s[a + lo] < s[b + lo]).
 * Notes: Why the bound holds: for two different strings the
 *   difference of the polynomials has at most $n - 1$ roots among
 *   the $M$ bases, so $q$ comparisons fail with probability at most
 *   $q n / M$. Never hash modulo $2^64$ by overflow: a Thue-Morse
 *   string of length 2048 and its complement collide for every odd
 *   base, and even bases only see the last 64 characters. Fixed
 *   bases are hackable, and a single modulus near $10^9$ fails by
 *   the birthday bound once about $10^5$ hashes are compared all
 *   against all. Palindrome test of `s[l, r)`: compare with
 *   `rev.get(n - r, n - l)` on the reversed string. Under point
 *   updates keep the hashes in a Fenwick or segment tree. LCP in
 *   $O(log n)$ (Example) sorts suffixes or rotations in
 *   $O(n log^2 n)$, a cheap substitute for a suffix array. 2D
 *   patterns: hash each row, then hash the column of row hashes
 *   with a second base. Cues: "count distinct substrings of length
 *   `k`" (sort the hashes), many "is `s[l, r)` equal to `t[x, y)`"
 *   queries, huge implicitly defined strings (a block repeated $q$
 *   times hashes to a geometric sum).
 * Source: Karp, Rabin (1987); $2^61 - 1$ multiply-reduce as in
 *   common competitive programming libraries.
 * Status: stress-tested against exact strings and an independent
 *   evaluation
 */
#pragma once

struct StrHash {
  typedef uint64_t ull;
  static constexpr ull M = (1ULL << 61) - 1;
  static ull norm(ull x) { return x >= M ? x - M : x; }
  static ull mul(ull a, ull b) {
    __uint128_t c = (__uint128_t)a * b;
    return norm((ull)(c >> 61) + (ull)(c & M));
  }
  static inline const ull B = [] {
    auto seed = chrono::steady_clock::now().time_since_epoch();
    return mt19937_64(seed.count())() % M;
  }();
  vector<ull> h, pw;
  template<class S> StrHash(const S& s)
      : h(sz(s) + 1), pw(sz(s) + 1, 1) {
    rep(i, 0, sz(s)) {
      h[i + 1] = norm(mul(h[i], B) + s[i]);
      pw[i + 1] = mul(pw[i], B);
    }
  }
  ull get(int l, int r) { // s[l, r)
    return norm(h[r] + M - mul(h[l], pw[r - l]));
  }
};
