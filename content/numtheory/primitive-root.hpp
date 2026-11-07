/**
 * Title: Primitive root
 * Description: The smallest primitive root `g` of a prime `p`.
 * Precondition: `p` prime, $p < 2^63$.
 * Time: One `factor(p - 1)` plus $O(g log p)$ per prime factor;
 *   the smallest root is at most 113 for $p < 10^8$.
 * Usage:
 *   ll g = primitive_root(998244353); // 3
 * Uses: numtheory/factor
 * Example:
 *   primitive_root(1'000'000'007); // 5
 *   primitive_root((1LL << 61) - 1); // 37
 *   // order of a mod p: the smallest e >= 1 with a^e = 1
 *   vector<ll> f = factor(p - 1);
 *   ll e = p - 1;
 *   for (ll q : f)
 *     if (powmod(a, e / q, p) == 1) e /= q;
 *   // k-th root of unity (NTT): needs k | p - 1
 *   ll w = powmod(g, (p - 1) / k, p);
 * Notes: Primitive roots modulo `n` exist iff `n` is 1, 2, 4,
 *   $p^k$ or $2 p^k$ for an odd prime `p`. If `g` is one mod `p`,
 *   then `g` or $g + p$ is one mod every $p^k$ ($g + p$ exactly
 *   when $g^(p - 1) equiv 1 thick (mod p^2)$); mod $2 p^k$ take the
 *   odd one of `g` and $g + p^k$. There are $phi(p - 1)$ primitive
 *   roots, and $g^k$ is one iff $gcd(k, p - 1) = 1$. The order of
 *   any `a` divides $p - 1$, and $x^k equiv 1$ has $gcd(k, p - 1)$
 *   solutions. Cues: NTT roots of unity, discrete logs that turn
 *   products into sums (a product over a multiset becomes a
 *   convolution of exponents, UCup S4-21 J), counting `k`-th power
 *   residues.
 * Source: Classic (Gauss); cp-algorithms "Primitive root".
 * Status: stress-tested against brute-force orders
 */
#pragma once
#include "factor.hpp"

ll primitive_root(ll p) {
  if (p == 2) return 1;
  vector<ll> f = factor(p - 1);
  f.erase(unique(all(f)), f.end());
  for (ll g = 2;; g++) {
    bool ok = true;
    for (ll q : f) ok &= powmod(g, (p - 1) / q, p) != 1;
    if (ok) return g;
  }
}
