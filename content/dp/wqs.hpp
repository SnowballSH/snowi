/**
 * Title: Aliens trick (WQS binary search)
 * Description: Minimum cost $F(k)$ using exactly `k` items when
 *   $F$ is convex in `k`, from a solver of the relaxed problem
 *   with no count limit where each item costs an extra `lam`.
 *   Finds the smallest `lam` whose solver count is at most `k`;
 *   then $F(k) = v - "lam" dot k$.
 * Precondition: `solve(lam)` returns `{v, c}`: the minimum `v` of
 *   cost + `lam` * count, and the smallest count `c` among
 *   solutions with value `v`. `k` is a feasible count, and every
 *   slope $F(j + 1) - F(j)$ lies in $[-"hi", -"lo"]$. $F$ fits in
 *   `ll`. Slopes can be as large as $F$ itself, so if
 *   `lam * count` can reach $2^63$, `v` must be `__int128`.
 * Time: $O(log("hi" - "lo"))$ calls to `solve`.
 * Usage:
 *   ll best = wqs(k, -C, C, solve); // |slopes| <= C
 * Tweak: Ties toward more items: search the largest `lam` with
 *   count `>= k`. At most `k` items: `lo = 0`. Smallest count in
 *   the solver: minimise $v dot 2^20 + c$ as one `ll` ($|v| < 2^42$).
 * Example:
 *   // Largest sum of exactly k pairwise non-adjacent a[i].
 *   // Pairs compare (value, count): min breaks ties toward
 *   // fewer items in every transition. For a = 0, M, 0, ...,
 *   // M, 0 the last slope is the whole sum: with n = 2e5 and
 *   // a[i] <= 1e9, lam * count reaches 1e19, so use __int128.
 *   typedef pair<__int128, ll> P;
 *   auto solve = [&](ll lam) {
 *     P skip{0, 0}, take{LLONG_MAX, 0};
 *     for (ll x : a) {
 *       P t{skip.first - x + lam, skip.second + 1};
 *       skip = min(skip, take), take = t;
 *     }
 *     return min(skip, take);
 *   };
 *   ll best = -wqs(k, -C, C, solve); // C = sum of |a[i]|
 * Notes: Why the tie rule: with ties toward fewer items, the
 *   reported count is the smallest optimal count, which never
 *   increases with `lam`, and it first drops to `k` or below at
 *   $"lam" = F(k) - F(k + 1)$, where `k` and `k + 1` are both
 *   optimal. So `k` is optimal there and $F(k) = v - "lam" k$ even
 *   if the solver reports fewer items; always use `k`, not `c`.
 *   Where $F$ is straight, many counts tie at one `lam`; a solver
 *   that breaks ties inconsistently (one comparison on values
 *   alone) makes the search overshoot to a `lam` where `k` is not
 *   optimal (see Tweak). Concavity sources: min-cost flow (`k`
 *   units), matroids (spanning tree with `k` white edges),
 *   partitions with a quadrangle-inequality cost (`dp/dc-opt`),
 *   exchange arguments.
 *   If unsure, print $F$ from a brute force and check
 *   $2 F(k) <= F(k - 1) + F(k + 1)$. Real-valued slopes: 60 rounds
 *   of bisection on a `double` `lam`. All $F(k)$ at once: query
 *   `lam` at chord slopes of known points; each call finds a new
 *   vertex or proves a segment straight.
 * Source: IOI 2016 "Aliens"; Packet 07 (UCup).
 * Status: stress-tested against exhaustive search
 */
#pragma once

template<class S> ll wqs(ll k, ll lo, ll hi, S solve) {
  while (lo < hi) {
    ll lam = lo + (hi - lo) / 2;
    if (solve(lam).second <= k) hi = lam;
    else lo = lam + 1;
  }
  return solve(lo).first - (__int128)lo * k;
}
