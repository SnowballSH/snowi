/**
 * Title: Hash map
 * Description: Seeded splitmix64 for `unordered_map` and
 *   `gp_hash_table`: the default integer hash is the identity and
 *   `gp_hash_table` masks it, so structured keys (multiples of
 *   $2^k$) collide.
 * Precondition: Integer keys (up to 64 bits) or pairs of them.
 *   With `gp_hash_table`, do not define `_GLIBCXX_DEBUG`: it then
 *   checks the whole table on every operation (400 keys: a minute).
 * Time: $O(1)$ expected per operation. $10^6$ random keys:
 *   `reserve` halves the time of `unordered_map`, `gp_hash_table`
 *   takes about a third of it (presizing saves a third more), and
 *   `map` takes 4 times as long.
 * Usage:
 *   unordered_map<ll, int, SafeHash> cnt;
 *   cnt.reserve(1 << 20); // before inserting ~1e6 keys
 *   gp_hash_table<ll, int, SafeHash> g; // faster
 *   unordered_map<pii, int, SafeHash> grid; // any pair
 * Tweak: Tuples: nest pairs (`pair<pii, ll>`) or chain through
 *   the seeded `(*this)` as the pair overload does. Preset the size
 *   of a `gp_hash_table` with the constructor arguments
 *   `({}, {}, {}, {}, {1 << 20})`.
 * Notes: Use it on judges with hacks (Codeforces) or when keys are
 *   structured (multiples of a power of two, grid coordinates):
 *   `gp_hash_table` hashes with a mask, so even ordinary keys like
 *   `x << 20` collide without a mixing hash. `gp_hash_table`
 *   has no `count` (use `find`) and no `reserve`. Iteration order
 *   is random and changes between runs with this seed. The pair
 *   overload is a template, so `pair<ll, ll>` is never truncated
 *   to `pii`. Never combine parts with the unseeded `mix`: with
 *   `mix(a) ^ b`, keys `(a, mix(a) ^ C)` all collide. A sorted
 *   vector with `lower_bound`, or coordinate compression into an
 *   array, is often faster than any hash map.
 * Source: splitmix64 (Vigna); the seeded hash from neal's
 *   Codeforces blog on hacking `unordered_map`; kactl (CC0).
 * Status: stress-tested against map; spread checked on bad keys
 */
#pragma once
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

struct SafeHash {
  static uint64_t mix(uint64_t x) { // splitmix64
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }
  size_t operator()(uint64_t x) const {
    static const uint64_t seed =
        chrono::steady_clock::now().time_since_epoch().count();
    return mix(x + seed);
  }
  template<class A, class B>
  size_t operator()(const pair<A, B>& p) const {
    return (*this)((*this)(p.first) ^ p.second);
  }
};
