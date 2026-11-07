/**
 * Title: Ordered set
 * Description: A `set` with order statistics, from the GNU
 *   policy-based data structures: `find_by_order(k)` is an
 *   iterator to the `k`-th smallest key (0-indexed, `end()` if
 *   `k >= size()`) and `order_of_key(x)` counts the keys `< x`.
 * Precondition: Keys are distinct under `less<T>`: inserting an
 *   existing key does nothing. Do not define `_GLIBCXX_DEBUG`:
 *   pb_ds then checks the whole tree on every operation (400
 *   inserts: 30 s).
 * Time: $O(log n)$ per operation.
 * Usage:
 *   Tree<int> t; t.insert(x); t.erase(x);
 *   int k = t.order_of_key(x); // keys < x
 *   int y = *t.find_by_order(k); // k-th smallest
 *   Tree<pii> ms; int id = 0; // multiset of ints
 *   ms.insert({x, id++});
 *   ms.erase(ms.lower_bound({x, INT_MIN})); // needs x in ms
 *   ms.order_of_key({x, INT_MIN}); // values < x
 * Tweak: Map with order statistics: `null_type` becomes the mapped
 *   type, then `t[key] = v` works. Never `less_equal` for a
 *   multiset: `find`, `erase` and `lower_bound` then misbehave.
 * Example:
 *   // Inversions online: elements before a[i] that exceed it.
 *   Tree<pii> seen;
 *   ll inv = 0;
 *   rep(i, 0, n) {
 *     inv += i - seen.order_of_key({a[i], INT_MAX});
 *     seen.insert({a[i], i});
 *   }
 *   // Key ranges of t and t2 must not overlap; O(log n):
 *   t.join(t2); // t2 becomes empty
 *   t.split(x, t2); // t2 is cleared, gets keys > x; O(|t2|)
 * Notes: Think of it for "the k-th smallest of a changing set",
 *   "how many stored values are below x", the median of a sliding
 *   window, or Josephus-style "remove the k-th remaining". Speed
 *   is close to `set`. If all values are known offline, a Fenwick
 *   tree over compressed values (`lower_bound` gives the k-th)
 *   does the same with less memory. `join` throws if the key
 *   ranges overlap; `split` counts the part it moves, so it is
 *   linear in that part, not logarithmic.
 * Source: kactl OrderStatisticTree (CC0); GNU libstdc++ pb_ds.
 * Status: stress-tested against a sorted vector
 */
#pragma once
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

template<class T> using Tree =
    tree<T, null_type, less<T>, rb_tree_tag,
         tree_order_statistics_node_update>;
