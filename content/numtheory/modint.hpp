/**
 * Title: Modular integer
 * Description: Integer modulo `MOD` with `+ - * /`, unary minus,
 *   `^` for powers (a negative exponent inverts first), `inv()`, `==`
 *   and stream output. Everything except output is `constexpr`, so
 *   small tables (up to about $10^4$ entries) can be built at compile
 *   time.
 * Precondition: `MOD` prime (for `/`, `inv` and negative powers) and
 *   below $2^31$, so products fit in `ll`. The constructor takes any
 *   `ll`, negative too.
 * Time: $O(1)$; `/` and `inv` take $O(log "MOD")$, `^` takes
 *   $O(log |e|)$.
 * Usage:
 *   Mint a = 5, b = (a ^ 3) / 2 - 1;
 *   cout << b * b.inv() << ' ' << b.x << '\n';
 * Tweak: Runtime modulus: declare `ll md = 1;` (nonzero, so a stray
 *   global `Mint` cannot divide by zero), use it instead of `MOD` and
 *   drop `constexpr`. Set `md` before making any `Mint`: a value made
 *   earlier, such as a global, stays reduced modulo the old `md`. For
 *   a composite `md`, let `inv()` return `inv_mod(x, md)` from
 *   Extended Euclid; it needs $gcd(x, "md") = 1$.
 * Example:
 *   constexpr auto fact = [] {
 *     array<Mint, 21> f{1};
 *     rep(i, 1, 21) f[i] = f[i - 1] * i;
 *     return f;
 *   }(); // fact[20] = 20! mod MOD, at compile time
 *   Mint s = 0;
 *   for (int k : {1, 2, 3}) s += Mint(k).inv(); // 1 + 1/2 + 1/3
 *   assert(s * 6 == 11);
 * Notes: `^` binds more weakly than `+`, `*` and even `==`: always
 *   write `(a ^ e)` in parentheses. `x / 0` silently gives 0. Read
 *   input into an `ll` and convert; compare with `.x` for `<`. Build
 *   big tables at run time: a compile-time factorial table of
 *   $2 dot 10^5$ entries takes g++ 15 s and 1.2 GB to compile.
 * Source: kactl ModularArithmetic (CC0), AtCoder Library modint
 *   (CC0).
 * Status: stress-tested against `ll` arithmetic, constexpr checked
 */
#pragma once
#include "modpow.hpp"

struct Mint {
  using M = Mint;
  ll x;
  constexpr Mint(ll v = 0) : x((v %= MOD) < 0 ? v + MOD : v) {}
  constexpr M& operator+=(M o) {
    if ((x += o.x) >= MOD) x -= MOD;
    return *this;
  }
  constexpr M& operator-=(M o) {
    if ((x -= o.x) < 0) x += MOD;
    return *this;
  }
  constexpr M& operator*=(M o) { return *this = x * o.x; }
  constexpr M& operator/=(M o) { return *this *= o.inv(); }
  friend constexpr M operator+(M a, M b) { return a += b; }
  friend constexpr M operator-(M a, M b) { return a -= b; }
  friend constexpr M operator*(M a, M b) { return a *= b; }
  friend constexpr M operator/(M a, M b) { return a /= b; }
  constexpr M operator-() const { return M() - *this; }
  constexpr M operator^(ll e) const {
    if (e < 0) return inv() ^ -e;
    M r = 1, b = *this;
    for (; e; e >>= 1, b *= b)
      if (e & 1) r *= b;
    return r;
  }
  constexpr M inv() const { return *this ^ (MOD - 2); }
  constexpr bool operator==(const M&) const = default;
};
ostream& operator<<(ostream& o, Mint a) { return o << a.x; }
