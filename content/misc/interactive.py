"""
Title: Interactive I/O
Description: Talks to an interactive judge. `ask(x, y)` prints
  the query line `? x y`, flushes, and reads one integer reply;
  `reply()` reads the next one; `answer(x)` prints `! x`. If the
  judge replies `-1` (an invalid query, or too many) or closes
  the pipe, the program exits at once with status 0, so the
  verdict is Wrong Answer instead of a time limit or a runtime
  error.
Precondition: Replies are integers separated by any whitespace.
  Every line must be flushed (`flush=True` does it); debug output
  goes to `sys.stderr` only. Read the judge only with `reply()` or
  `input()`, never with `sys.stdin.buffer` (as the template does):
  that waits forever.
Time: Each query is a round trip through the judge, measured at
  about 20 µs of CPU and 25 to 40 µs of wall time: $10^5$
  queries take 3 to 4 s of wall time.
Usage:
  n = reply()  # the judge's input, too
  r = ask(u, k)  # prints "? u k", reads one reply
  a, b = ask(i), reply()  # a reply with two numbers
  answer(x)  # multi-test with a verdict: reply() after it
Tweak: Word replies (`<`, `YES`): return `r`, not `int(r)`; the
  `-1` check still works. If `-1` is a legal reply, drop the
  check. Lists: `ask(*v)`, or `ask(len(v), *v)`.
Example:
  >>> def solve(n):  # hidden x in [1, n]; "? y": is x <= y?
  ...   lo, hi = 1, n  # invariant: lo <= x <= hi
  ...   while lo < hi:
  ...     mid = (lo + hi) // 2  # never asks about hi
  ...     if ask(mid):
  ...       hi = mid
  ...     else:
  ...       lo = mid + 1
  ...   answer(lo)
  >>> import io, sys  # test offline with canned replies:
  >>> sys.stdin = io.StringIO("1 0")  # x = 2
  >>> solve(4)
  ? 2
  ? 1
  ! 2
Notes: Count first: $N$ hidden possibilities and $k$ possible
  replies need $ceil(log_k N)$ queries in the worst case. A limit
  equal to that bound means every query must split the
  candidates evenly; about $1.7$ times it allows 1/3-2/3 splits
  (centroids, subtrees of a binary tree); a limit linear in $n$
  means $O(1)$ queries per element. Binary search needs an
  invariant, not a monotone array: keep `lo` and `hi` with
  opposite properties (`a[lo] = 0`, `a[hi] = 1`) and let `mid`
  replace the one it matches; this finds a sign change or a
  local minimum (compare `a[mid]` with `a[mid + 1]`). Never
  spend a query on a value you know, such as an endpoint.
  Adaptive judges may change the hidden object as long as the
  replies stay consistent: answer only when one candidate is
  left, and do not rely on average-case (random) query counts.
  Sorting with a comparator that asks the judge:
  `sorted(a, key=cmp_to_key(cmp))` (from `functools`) uses about
  $n log_2 n$ comparisons. Local testing: an interactor reads
  queries on stdin, writes replies on stdout and the verdict on
  stderr; run both with `mkfifo io; ./judge < io | pypy3 a.py >
  io`. A hang means a missing flush or a read the other side
  never answers. The best local judge is adversarial: it keeps
  all candidates and replies so that the larger part survives.
  Count your queries and assert the limit. `input()` and `reply()`
  share one buffer, but `sys.stdin.buffer.readline()` does not:
  after `input()` it waits for lines that `input()` already took.
  Group testing ("does X contain a marked item?"): test one
  half, keep it if positive and the other half if not, until one
  item is left; remove it and test the rest again: $k$ marked
  items cost at most $k (ceil(log_2 n) + 1) + 1$ tests. An
  independence oracle ("is X + v independent?") answers "does v
  have a neighbour in X": grow a maximal independent set, find
  each rejected vertex's neighbours in it by group testing (skip
  its first test, on the whole set: it is known positive),
  recurse on the rejected ones; at most
  $n + m (2 + ceil(log_2 n))$ queries. When every reply is a
  known function of the hidden object, keep all consistent
  candidates and ask the query that splits them best; telling
  $k$ behaviours apart needs at most $k - 1$ tests (add one only
  when it splits a group).
Source: UCup interactive training packet.
Status: run against a judge process through pipes
"""

import sys

judge_tokens = []


def reply():
  while not judge_tokens:
    line = sys.stdin.readline()  # the buffer input() uses
    if not line:
      sys.exit(0)  # the judge closed the pipe
    judge_tokens.extend(line.split()[::-1])
  r = judge_tokens.pop()
  if r == "-1":
    sys.exit(0)  # the judge rejected a query
  return int(r)


def ask(*xs):
  print("?", *xs, flush=True)
  return reply()


def answer(*xs):
  print("!", *xs, flush=True)
