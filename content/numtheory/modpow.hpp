/**
 * Title: Modular arithmetic
 * Description: Fast power and inverse modulo the prime `MOD`. Other
 *   snippets include this file for `MOD` and `modpow`.
 * Precondition: `b >= 0`, `e >= 0`; `MOD` prime for `modinv`.
 * Time: $O(log e)$.
 * Usage:
 *   ll x = modpow(3, 1'000'000'000'000LL), y = modinv(x);
 * Tweak: For a modulus above $2^31$ multiply in `__int128`.
 * Status: stress-tested against repeated multiplication
 */
#pragma once

const ll MOD = 998244353;
ll modpow(ll b, ll e) {
  ll r = 1;
  for (b %= MOD; e; e >>= 1, b = b * b % MOD)
    if (e & 1) r = r * b % MOD;
  return r;
}
ll modinv(ll a) { return modpow(a, MOD - 2); }
