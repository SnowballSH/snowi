/**
 * Title: Modular square root
 * Description: A square root of `a` modulo a prime `p` by
 *   Tonelli-Shanks: returns the smaller root $r <= p - r$ with
 *   $r^2 equiv a thick (mod p)$, or -1 if `a` is not a square.
 * Precondition: `p` prime, $p < 2^63$; `a` any integer.
 * Time: $O(log^2 p)$ worst case, when $2^k | p - 1$ for large $k$;
 *   about 3 µs for $p = 998244353$.
 * Usage:
 *   ll r = mod_sqrt(a, p); // -1, or roots r and p - r
 * Tweak: For $p equiv 3 thick (mod 4)$ the root is simply
 *   $a^((p + 1)\/4)$ (check its square). `powmod` and `mulmod` come
 *   from Primality and factoring; only those two are needed.
 * Uses: numtheory/factor
 * Example:
 *   mod_sqrt(2, 7); // 3: 3^2 = 9 = 2, and 4^2 = 16 = 2
 *   mod_sqrt(3, 7); // -1: the squares mod 7 are 0, 1, 2, 4
 *   mod_sqrt(0, 5); // 0
 *   mod_sqrt(2, 998244353); // 116195171
 * Notes: The roots are `r` and $p - r$, distinct unless $r = 0$ or
 *   $p = 2$.
 *   Legendre symbol: $a^((p - 1)\/2) mod p$ is 1 for a nonzero
 *   square, $p - 1$ for a non-square, 0 for $a equiv 0$; half of the
 *   nonzero residues are squares. Solve $x^2 + b x + c equiv 0$ for
 *   odd `p` as $x = (-b plus.minus r)\/2$ with $r^2 equiv b^2 - 4c$.
 *   Modulo $p^k$ ($p$ odd, $p divides.not a$) lift a root with
 *   $r <- r - (r^2 - a)\/(2r)$, which doubles the correct digits;
 *   for a composite modulus combine prime powers with CRT. Cues: the
 *   constant term of a power series square root, "is $x^2 equiv a$
 *   solvable", $sqrt(5)$ for Binet's Fibonacci formula (for odd `p`
 *   it exists iff $p equiv plus.minus 1 thick (mod 5)$ or $p = 5$).
 * Source: Tonelli (1891), Shanks (1973).
 * Status: stress-tested against all squares modulo small primes
 */
#pragma once
#include "factor.hpp"

ll mod_sqrt(ll a, ll p) {
  if ((a %= p) < 0) a += p;
  if (a < 2) return a; // also p = 2
  if (powmod(a, (p - 1) / 2, p) != 1) return -1;
  ll s = __builtin_ctzll(p - 1), q = (p - 1) >> s, z = 2;
  while (powmod(z, (p - 1) / 2, p) == 1) z++;
  ll c = powmod(z, q, p), x = powmod(a, (q + 1) / 2, p),
     t = powmod(a, q, p);
  while (t != 1) { // x^2 = a t, t has order 2^i < 2^s
    ll i = 0, u = t;
    for (; u != 1; i++) u = mulmod(u, u, p);
    ll b = powmod(c, 1LL << (s - i - 1), p);
    x = mulmod(x, b, p), c = mulmod(b, b, p);
    t = mulmod(t, c, p), s = i;
  }
  return min(x, p - x);
}
