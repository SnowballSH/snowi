# Writing snowi content

Everything in the PDFs is generated from files under `content/`. This file is
the contract those files follow; `make lint` enforces most of it.

## Layout

```
content/<section>/section.toml   order of items and notes, ICPC selection
content/<section>/<item>.hpp     C++ snippet (a .cpp only for full programs)
content/<section>/<item>.py      PyPy snippet
content/<section>/<note>.typ     prose: theorems, formulas, guides
tests/<section>/<item>.test.cpp  C++ stress test
tests/<section>/test_<item>.py   PyPy stress test ('-' becomes '_')
```

Item ids are `<section>/<name>`, lowercase with hyphens (`ds/lazy-segtree`).

## section.toml

```toml
title = "Data structures"

[[entry]]
item = "segtree"          # code item: segtree.hpp + segtree.py
icpc = true               # also printed in the 25-page ICPC packet

[[entry]]
item = "ordered-set"
langs = ["cpp"]           # only when a language genuinely has no
                          # version (a note explains the alternative)

[[entry]]
note = "monoids"          # monoids.typ
title = "Designing a monoid"
icpc = false
```

Entries print in this order. `icpc` defaults to false; `icpc = ["cpp"]`
puts an item in only the C++ ICPC packet (PyPy code is longer, so the two
25-page packets fill at different rates).

## Header block

Every snippet starts with a header. C++ uses a `/** ... */` block, Python a
module docstring. The fields are the same in both languages:

```cpp
/**
 * Title: Segment tree
 * Description: Point assign and range fold over a monoid. Also finds
 *   the largest `r` with `pred(query(l, r))` true.
 * Precondition: `op` is associative and `id` is its identity.
 * Time: $O(log n)$ per operation.
 * Usage:
 *   SegTree st(n);
 *   st.set(i, 5); st.query(l, r); // [l, r)
 * Tweak: For range sum use `T = ll`, `id = 0`, `op = a + b`.
 * Example:
 *   ...longer example, printed only in the full packet...
 * Notes: Pitfalls and problem cues, printed only in the full packet.
 * Source: Idea from the AtCoder Library (CC0).
 * Status: stress-tested against brute force
 */
#pragma once
```

```python
"""
Title: Segment tree
Description: Point assign and range fold over a monoid.
Time: $O(log n)$ per operation.
Usage:
  st = SegTree([0] * n, max, -INF)
  st.set(i, 5); st.query(l, r)  # [l, r)
Example:
  >>> st = SegTree([3, 1, 4], max, -1)
  >>> st.query(0, 2)
  3
Status: stress-tested against brute force
"""
```

| Field | Required | Printed in | Content |
|---|---|---|---|
| Title | yes | both | Short name. Same in both languages. |
| Description | yes | both | What it computes, in one to three plain sentences. |
| Precondition | no | both | Input requirements: indexing, ranges, graph shape. |
| Time | no | both | Complexity, e.g. `$O(n log n)$`. Add memory if notable. |
| Usage | yes | both | Code (verbatim, monospace). The minimal calls. |
| Tweak | no | both | How to change types, monoids, the modulus. |
| Uses | no | both | Other items this depends on, e.g. `ds/dsu`. C++ local includes are added automatically. |
| Example | no | full | Code. In Python, `>>>` lines are run as doctests. |
| Notes | no | full | Pitfalls, recognition cues, variants. |
| Source | no | full | Where the idea comes from. |
| Status | yes | never | How it is verified. |

A field starts at a line `Name: text`. Following lines belong to it until the
next field. In text fields lines are joined with spaces. In code fields
(Usage, Example) the lines are kept and dedented.

Text fields are plain text with two kinds of markup: `` `code` `` spans and
`$math$` spans in Typst math syntax (`$O(n log n)$`, `$sum_i a_i$`,
`$n^2$`, `$floor(n/2)$`). Nothing else is interpreted.

## Code style

- Lines are at most 64 characters, including in Usage and Example. The
  packet layout fits exactly 64 monospace characters per column.
- Run `make format`. C++ follows `.clang-format`: 2-space indent, short
  blocks on one line. Python follows `ruff.toml`: 2-space indent.
- C++ targets g++ 13 with `-std=gnu++20`. Snippets assume the contest
  template (`content/contest/template.cpp`): `bits/stdc++.h`,
  `using namespace std`, `rep`, `all`, `sz`, `ll`, `pii`, `vi`.
  Use other snippets with `#include "../ds/dsu.hpp"`; the include is not
  printed, the dependency is listed instead.
- Python targets PyPy 7.3.15 (Python 3.9): no `match`, no `X | Y` types,
  no new 3.10+ library functions. Snippets must be pasteable: no imports
  between snippets (declare them with `Uses:`), no top-level work, no stdin
  reads. Import stdlib modules inside the snippet when needed.
- Readable names over single letters, except loop indices and the usual
  `n`, `m`, `u`, `v`. Concrete types are fine; say in Tweak how to change
  them.
- Comment only non-obvious lines. Explanations belong in Description or
  Notes.
- PyPy speed: flat lists of ints, no recursion deeper than about 1000
  (use explicit stacks), avoid per-element tuples and closures in hot
  loops, read input with `sys.stdin.buffer`.

## Tests

Every code item has a stress test per language that compares the snippet
with an independent brute force on many random small inputs, plus edge
cases (n = 0 or 1, equal values, extremes) and one large input as a speed
sanity check.

- C++: `#include "../testlib.hpp"` then the snippet. `testlib.hpp`
  provides the template prelude, `rnd(lo, hi)`, `rnd_vec`, `CHECK(cond)`
  and `CHECK_EQ(a, b)`. Print one summary line, return 0 on success.
  Tests compile with `-std=gnu++20 -O1 -g -Wall -Wextra -Wshadow
  -fsanitize=address,undefined -D_GLIBCXX_ASSERTIONS` (warnings fail the
  build) and must finish in about 10 seconds.
- Python: `from testlib import load, rnd, check, check_eq`;
  `lib = load("ds/segtree")` executes the snippet (and its `Uses:`) into
  one namespace, exactly like pasting. Runs under `pypy3`, about 20 seconds.
- A test must fail when the snippet is wrong. After writing one, break the
  snippet on purpose (an off-by-one, a wrong comparison) and confirm the
  test catches it.

Run one item with `make test ITEM=ds/segtree`, everything with `make test`.

## Notes (.typ)

Notes are Typst markup included into the packet. Start with
`#import "/typst/notes.typ": *`. The note's title comes from section.toml,
so do not repeat it. Use `===` for subheadings and these helpers:

- `#full[...]` shows only in the full packet, `#icpc[...]` only in the
  ICPC one. `#cpp[...]` and `#py[...]` show only in that language's PDFs.
- `#fact[Name][Statement]` for a compact named theorem or formula.
- `#code("cpp", "...")` for a short snippet inside a note.

Keep it terse and in plain language: state the fact, the condition, and
when it helps. Prefer formulas and short examples over paragraphs.
