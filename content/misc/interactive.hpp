/**
 * Title: Interactive I/O
 * Description: Talks to an interactive judge. `ask(x, y)` prints
 *   the query line `? x y`, flushes, and reads one integer reply;
 *   `reply()` reads the next one; `answer(x)` prints `! x`. If the
 *   judge replies `-1` (an invalid query, or too many) or closes
 *   the pipe, the program exits at once with status 0, so the
 *   verdict is Wrong Answer instead of a time limit or a runtime
 *   error. `reply` turns off the template's `cin.exceptions`.
 * Precondition: Replies are integers separated by any whitespace.
 *   Every line must be flushed (`endl` does it); debug output goes
 *   to `cerr` only.
 * Time: Each query is a round trip through the judge, measured at
 *   about 8 µs of CPU and 15 to 25 µs of wall time: $10^5$
 *   queries are fine.
 * Usage:
 *   ll r = ask(u, k); // prints "? u k", reads one reply
 *   ll a = ask(i), b = reply(); // a reply with two numbers
 *   answer(x); // multi-test with a verdict: reply() after it
 * Tweak: Word replies (`<`, `YES`): let `reply` read and return a
 *   `string`; test `r == "-1"`. If `-1` is a legal reply, drop the
 *   check. To send a list, define `ostream& operator<<(ostream&,
 *   const vi&)` printing the elements with spaces before this
 *   snippet: `say` does not see one defined after it.
 * Example:
 *   // Hidden x in [1, n]; "? y" replies 1 iff x <= y.
 *   int n = reply();
 *   ll lo = 1, hi = n; // invariant: lo <= x <= hi
 *   while (lo < hi) {
 *     ll mid = lo + (hi - lo) / 2; // never asks about hi
 *     if (ask(mid)) hi = mid;
 *     else lo = mid + 1;
 *   }
 *   answer(lo);
 * Notes: Count first: $N$ hidden possibilities and $k$ possible
 *   replies need $ceil(log_k N)$ queries in the worst case. A limit
 *   equal to that bound means every query must split the
 *   candidates evenly; about $1.7$ times it allows 1/3-2/3 splits
 *   (centroids, subtrees of a binary tree); a limit linear in $n$
 *   means $O(1)$ queries per element. Binary search needs an
 *   invariant, not a monotone array: keep `lo` and `hi` with
 *   opposite properties (`a[lo] = 0`, `a[hi] = 1`) and let `mid`
 *   replace the one it matches; this finds a sign change or a
 *   local minimum (compare `a[mid]` with `a[mid + 1]`). Never
 *   spend a query on a value you know, such as an endpoint.
 *   Adaptive judges may change the hidden object as long as the
 *   replies stay consistent: answer only when one candidate is
 *   left, and do not rely on average-case (random) query counts.
 *   Sorting with a comparator that asks the judge: `stable_sort`
 *   (merge sort, about $n log_2 n$ comparisons); `sort` uses more
 *   and can run off the array if the replies are inconsistent.
 *   Local testing: an interactor reads queries on stdin, writes
 *   replies on stdout and the verdict on stderr; connect the two
 *   with a named pipe (command in Contest setup). A hang means a
 *   missing flush or a read the other side never answers. The
 *   best local judge is adversarial: it keeps all candidates and
 *   replies so that the larger part survives. Count your queries
 *   and assert the limit.
 *   Group testing ("does X contain a marked item?"): test one
 *   half, keep it if positive and the other half if not, until
 *   one item is left; remove it and test the rest again: $k$
 *   marked items cost at most $k (ceil(log_2 n) + 1) + 1$
 *   tests. An independence oracle ("is X + v independent?")
 *   answers "does v have a neighbour in X": grow a maximal
 *   independent set, find each rejected vertex's neighbours in
 *   it by group testing (skip its first test, on the whole set:
 *   it is known positive), recurse on the rejected ones; at
 *   most $n + m (2 + ceil(log_2 n))$ queries. When every reply
 *   is a known function of the hidden object, keep all
 *   consistent candidates and ask the query that splits them
 *   best; telling $k$ behaviours apart needs at most $k - 1$
 *   tests (add one only when it splits a group).
 * Source: UCup interactive training packet.
 * Status: run against a judge process through pipes
 */
#pragma once

ll reply() {
  ll r;
  cin.exceptions(cin.goodbit);         // EOF must not throw
  if (!(cin >> r) || r == -1) exit(0); // rejected or closed
  return r;
}
template<class... T> void say(char c, T... xs) {
  cout << c;
  ((cout << ' ' << xs), ...);
  cout << endl; // endl flushes
}
template<class... T> auto ask(T... xs) {
  say('?', xs...);
  return reply();
}
template<class... T> void answer(T... xs) { say('!', xs...); }
