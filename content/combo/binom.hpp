/**
 * Title: Binomial coefficients
 * Description: Factorials modulo `MOD` up to `n` give $binom(n, k)$
 *   (0 outside $0 <= k <= n$), $n!\/(n - k)!$, multinomials, Catalan
 *   numbers and $1\/i$ in $O(1)$. `Lucas` gives $binom(n, k)$ modulo
 *   a small prime `p` for huge `n` and `k`.
 * Precondition: `MOD` prime and $> n$. Arguments at most `n` (in
 *   `multinom` their sum): `catalan(k)` needs $2k <= n$, `inv(i)`
 *   needs $i >= 1$. `Lucas`: `p` prime up to about $10^7$;
 *   $n, k <= 10^18$.
 * Time: $O(n)$ to build, $O(1)$ per query. `Lucas`: $O(p)$ to
 *   build, $O(log_p n)$ per query.
 * Usage:
 *   Binom bc(2 * n);
 *   bc.C(5, 2), bc.perm(5, 2), bc.catalan(3); // 10 20 5
 *   bc.multinom({2, 1, 1}), bc.inv(4);         // 12, 1/4
 *   Lucas lu(10007);
 *   lu.C(ll(1e18), 12345);
 * Tweak: Change `MOD` for another prime. `Binom<int> bc(n)` halves
 *   the memory: 80 MB instead of 160 MB for $n = 10^7$ (products
 *   stay `ll`). A prime power or composite modulus has no inverse
 *   factorials: use Pascal's triangle in $O(n^2)$, or CRT over prime
 *   moduli.
 * Example:
 *   // Monotone paths (0, 0) -> (h, w) avoiding blocked cells,
 *   // sorted by (row, col); the target is the last cell.
 *   auto paths = [&](pii a, pii b) {
 *     int dr = b.first - a.first, dc = b.second - a.second;
 *     return bc.C(dr + dc, dr); // 0 if b is not reachable
 *   };
 *   cells.push_back({h, w});
 *   vector<ll> dp(sz(cells));
 *   rep(i, 0, sz(cells)) {
 *     dp[i] = paths({0, 0}, cells[i]);
 *     rep(j, 0, i) dp[i] = (dp[i] + MOD - dp[j] *
 *                           paths(cells[j], cells[i]) % MOD) % MOD;
 *   }
 *   // dp.back(): good paths. Stars and bars: ways to put n
 *   // equal balls into k boxes = bc.C(n + k - 1, k - 1).
 * Notes: Size the table for the largest argument: Catalan and
 *   stars and bars read up to $2n$ or $n + k$. For huge `n` and small
 *   `k` multiply $n (n - 1) dots (n - k + 1)$ by `ifact[k]` in $O(k)$.
 *   `C` returns 0 out of range, which handles unreachable cells and
 *   negative counts in inclusion-exclusion sums, but stars and bars
 *   with $n = k = 0$ needs the answer 1 by hand. A modulus $<= n$
 *   (like 10007) makes factorials vanish: use `Lucas`. Mod 2,
 *   $binom(n, k)$ is odd iff `(k & ~n) == 0`. A member or parameter
 *   needs the full type, `Binom<>`.
 * Source: Standard; Wilson's theorem seeds the inverse table of
 *   `Lucas`.
 * Status: stress-tested against Pascal's triangle and exact
 *   binomials
 */
#pragma once
#include "../numtheory/modpow.hpp"

template<class T = ll> struct Binom {
  vector<T> fact, ifact;
  Binom(int n) : fact(n + 1, 1), ifact(n + 1) {
    rep(i, 1, n + 1) fact[i] = 1LL * fact[i - 1] * i % MOD;
    ifact[n] = modinv(fact[n]);
    for (int i = n; i; i--)
      ifact[i - 1] = 1LL * ifact[i] * i % MOD;
  }
  ll C(int n, int k) {
    if (k < 0 || k > n) return 0;
    return 1LL * fact[n] * ifact[k] % MOD * ifact[n - k] % MOD;
  }
  ll perm(int n, int k) {
    if (k < 0 || k > n) return 0;
    return 1LL * fact[n] * ifact[n - k] % MOD;
  }
  ll multinom(const vi& ks) { // (sum k)! / prod k!
    ll r = 1;
    int s = 0;
    for (int k : ks) r = r * ifact[k] % MOD, s += k;
    return r * fact[s] % MOD;
  }
  ll catalan(int k) {
    return (C(2 * k, k) - C(2 * k, k + 1) + MOD) % MOD;
  }
  ll inv(int i) { return 1LL * ifact[i] * fact[i - 1] % MOD; }
};

struct Lucas {
  ll p;
  vector<ll> fact, ifact;
  // Wilson: (p - 1)! = -1 mod p seeds ifact[p - 1].
  Lucas(ll prime) : p(prime), fact(p, 1), ifact(p, p - 1) {
    rep(i, 1, p) fact[i] = fact[i - 1] * i % p;
    for (ll i = p - 1; i; i--) ifact[i - 1] = ifact[i] * i % p;
  }
  ll C(ll n, ll k) {
    if (k < 0 || k > n) return 0;
    ll r = 1;
    for (; k; n /= p, k /= p) {
      ll a = n % p, b = k % p;
      if (b > a) return 0;
      r = r * fact[a] % p * ifact[b] % p * ifact[a - b] % p;
    }
    return r;
  }
};
