#import "/typst/notes.typ": *

#full[
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
- An index out of range#cpp[ (the debug build catches it)], an empty
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
]

#full[
=== Memory limit
#cpp[
  - A `vector<vector<int>>` row has 24 bytes of overhead plus the
    allocation.
]
- Store only what later steps need; reuse buffers across tests.
]

#cpp[
=== C++ pitfalls
- Overflow-free test: $a b > c$ iff `a > c / b` ($b > 0$,
  $c >= 0$).
- `pow` and `(ll)sqrt(x)` are inexact above $2^53$: $x = 10^18 - 1$
  gives $10^9$. Fix roots up and down with integer loops; compute
  powers with integers.
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
  stay small); keep bitmasks within 62 bits.
- Memory: a list of ints, or of floats, takes 8 bytes per
  element; with `None` mixed in (use `-1` instead), or floats with
  ints of $2^31$ or more, about 32. A list of $10^6$ lists of 10
  ints took 160 MB.
]
