#import "/typst/notes.typ": *

=== Before submitting
- Run the samples and your own tiny and extreme cases: $n = 1$,
  all values equal, sorted input, the largest values.
- Time a maximal case#py[ (with #cmd-inline("{pypy}"), not #cmd-inline("{python}"))] and
  estimate the memory.
- Can anything overflow? Is the output format exact (case,
  spaces, precision, one answer per test)? Is it the right file?

=== Wrong answer
- Read the statement again: limits, indexing, several test cases,
  what to print when there is no answer.
- Reset all global state between test cases.
- Overflow, uninitialised values, swapped `n` and `m` or `i` and
  `j`, an off-by-one at a boundary, a wrong loop direction.
- Stress test against a brute force (Setup); run a small case by
  hand. Add assertions.
- Is the idea right? Explain it to a teammate; recheck the special
  cases (empty, $n = 1$, ties, everything equal, no solution).
- Floating point: compare with a tolerance, print enough digits.
  Floors, ceilings and comparisons of fractions: stay in integers
  (cross-multiply, divide once at the end); a double misjudges an
  exact boundary.
- Open or closed ends, half-integer points: double every
  coordinate so that all relevant points are integers.
- Still stuck: let a teammate read the code, take a short walk,
  or rewrite the solution from scratch.

=== Runtime error
- An index out of range (the debug build catches it), an empty
  container (#cpp[`v.back()`, `*s.begin()`]#py[`a[-1]`, `max([])`]),
  division or modulo by 0, an assertion.
- Deep or infinite recursion; too much memory can also end as a
  runtime error.
#cpp[
  - Iterators or references invalidated by `push_back`, `erase`
    or `insert`.
]

=== Time limit
- An infinite loop: a pointer that does not move, a binary search
  that does not shrink.
- Complexity against the limits, then constants: copies of large
  objects#cpp[ (pass by reference)], hash maps where sorting would
  do, I/O.
- Many test cases: never reset an array of the maximum size per
  test; clear only what the test used.

=== Memory limit
- #cpp[`int` is 4 bytes, `ll` 8; a `vector<vector<int>>` row has
  24 bytes of overhead plus the allocation.]#py[A PyPy list of
  ints, or of floats, takes 8 bytes per element; with `None` mixed
  in (or floats with ints of $2^31$ or more) it takes about 32. A
  list of $10^6$ lists of 10 ints took 160 MB.]
- Store only what later steps need; reuse buffers across tests.

#cpp[
=== C++ pitfalls
- `int` products overflow before they are stored: with `int a, b`,
  `ll c = a * b` is wrong; write `1LL * a * b`.
  `accumulate(all(v), 0)` sums in `int`: pass `0LL`.
- `ll` ends at $9.2 dot 10^18$: $2 dot 10^5 dot 10^9 dot 10^9$
  needs `__int128` (which `cout` cannot print) or a test by
  division: $a b > c$ iff `a > c / b` for $b > 0$, $c >= 0$.
- `1 << k` is an `int`: write `1LL << k` for $k >= 31$.
  `__builtin_popcount` takes 32 bits: use `__builtin_popcountll`.
- `v.size() - 1` is $2^64 - 1$ for an empty `v`: write
  `sz(v) - 1`, which is a signed `int`.
- `endl` flushes: $10^6$ lines took 10 times longer than with
  `'\n'`. Do not mix `cin`/`cout` with `scanf`/`printf` after
  `sync_with_stdio(0)`.
- Local variables and arrays start uninitialised; globals are 0.
- A sort comparator must be a strict weak order: with `<=`, 100
  equal elements made `sort` read out of bounds.
- `-7 % 3 == -1` and `-7 / 2 == -3`: normalise with
  `(a % m + m) % m` and round negative divisions down by hand.
- `(ll)sqrt(x)` and `pow` are inexact above $2^53$:
  `(ll)sqrt(999999999999999999)` is $10^9$. Fix roots up and down
  with integer loops; compute powers with integers.
- `m[k]` inserts `k` when you only read it: use `find` or
  `count`. `min(a, b)` needs one type: `min<ll>(a, b)`.
]

#py[
=== PyPy pitfalls
- Recursion: PyPy crashes (exit code 139) after about $10^4$
  nested calls, even with `sys.setrecursionlimit`. The thread below
  reached $3 dot 10^5$ levels, but each level costs about 1.2 KB of
  memory and deep calls run several times slower. Prefer an
  explicit stack: record a DFS order, then process it in reverse
  for subtree values.
#code("py", "import sys, threading  # main() may recurse deeply
sys.setrecursionlimit(1 << 30)
threading.stack_size(1 << 28)
threading.Thread(target=main).start()")
- `s += c` in a loop is quadratic: $10^6$ characters took over a
  minute. Append to a list and `"".join` it at the end.
- Dicts and sets with tuple keys are about 3 times slower than
  with ints: encode `(u, v)` as `u * N + v`.
- Ints of $2^63$ and more are big ints, about 10 times slower:
  reduce mod $p$ after every product (two factors below $2^31$
  stay small).
- `input()` per line is about twice as slow as reading everything
  with `sys.stdin.buffer.read().split()`; write the output with
  one `"\n".join`.
- `list.pop(0)` and `x in list` take linear time: use
  `collections.deque` and sets.
- A list mixing ints with `None` is about 4 times larger and
  several times slower than a list of ints only: use `-1` for
  "none". Ints below $2^31$ mix with floats at little cost.
]
