/**
 * Title: Discrete logarithm
 * Description: The smallest $x >= 0$ with
 *   $a^x equiv b thick (mod m)$, or -1 if there is none, for any
 *   modulus, with $0^0 = 1$. Divides out the common factors of `a`
 *   and `m`, then runs baby-step giant-step with a hash map.
 * Precondition: $1 <= m <= 3 dot 10^9$, so products fit in `ll`;
 *   `a` and `b` may be any integers.
 * Time: $O(sqrt(m))$ hash operations: about 5 ms for $m = 10^9$.
 * Usage:
 *   ll x = discrete_log(3, 13, 17); // 4: 3^4 = 81 = 13 mod 17
 * Tweak: For $m$ up to about $10^12$ multiply with `mulmod` from
 *   Primality and factoring. Many queries with one `a` coprime to
 *   `m`: map $a^j -> j$ for $j < T$ once, keeping the smallest `j`
 *   for repeated values; a query sets $c = b$ and multiplies `c` by
 *   $a^(-T)$ until it is in the map, at most $m\/T$ times. After `i`
 *   steps $c = a^j$ and the answer is $x = i T + j$.
 * Example:
 *   discrete_log(2, 0, 8);  // 3
 *   discrete_log(2, 1, 8);  // 0
 *   discrete_log(2, 3, 8);  // -1
 *   discrete_log(2, 8, 12); // 3: 2^x mod 12 = 1, 2, 4, 8, 4, 8..
 *   // order of a, gcd(a, m) = 1 (inv_mod from Extended Euclid)
 *   ll ord = 1 + discrete_log(a, inv_mod(a, m), m);
 * Notes: When $gcd(a, m) = 1$ all solutions are
 *   $x + t op("ord")(a)$. Cues: "after how many multiplications by
 *   `a`", or a generator $x_(i + 1) = a x_i + c$ modulo a prime: with
 *   $d = c\/(a - 1)$ the values $x_i + d = a^i (x_0 + d)$ form a
 *   geometric sequence (handle $a = 1$ and $x_0 + d equiv 0$ apart).
 *   With a primitive root `g` of a prime `p`, logs turn products into
 *   sums modulo $p - 1$: $x^k equiv b$ becomes the linear congruence
 *   $k y equiv log_g b$ for $x = g^y$. Modulo $2^k$ with large `k`
 *   skip BSGS. Even `a`: $a^x equiv 0$ once $x >= k$, so try
 *   $x <= k$. Odd `a`: the order is a power of 2 and
 *   $v_2 (a^(2^j) - 1)$ strictly grows with `j`, so fix the bits of `x` from
 *   low to high (UCup S4-9 A "Cipher").
 * Source: Shanks (1971); gcd reduction as on cp-algorithms; kactl
 *   ModLog (CC0).
 * Status: stress-tested against iterating the powers
 */
#pragma once

ll discrete_log(ll a, ll b, ll m) {
  a = (a % m + m) % m, b = (b % m + m) % m;
  ll k = 1 % m, add = 0, g;
  while ((g = gcd(a, m)) > 1) { // now k * a^(x - add) = b
    if (b == k) return add;
    if (b % g) return -1;
    b /= g, m /= g, add++, k = k * (a / g) % m;
  }
  ll n = sqrt(m) + 1, an = 1, cur = b;
  unordered_map<ll, ll> baby; // b a^q -> largest q <= n
  rep(q, 0, n + 1) baby[cur] = q, cur = cur * a % m;
  rep(i, 0, n) an = an * a % m;
  cur = k;
  rep(p, 1, n + 1) {
    cur = cur * an % m;
    if (auto it = baby.find(cur); it != baby.end())
      return n * p - it->second + add;
  }
  return -1;
}
