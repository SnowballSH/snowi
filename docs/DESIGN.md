# Design

## One source, four packets

Each item is one directory entry with a C++ file and a PyPy file. The
documentation lives in a header block at the top of each file, so a file
copied out of the repository still explains itself, and so the PDF cannot
drift from the code: `tools/build.py` parses the headers, strips them,
and hands Typst a manifest (`build/manifest.json`) with the fields, the
printable code and a typing hash.

`section.toml` fixes the order and whether an entry is in the 25-page
ICPC packet. The full packet prints everything, plus the `Example`,
`Notes` and `Source` fields.

## Typst

- `typst/bundle.typ` uses the experimental bundle export of Typst 0.15:
  one compile emits all four PDFs, sharing parsing and layout work.
  `make pdf-separate` compiles them one by one with `sys.inputs` instead,
  should bundle export change.
- Notes are Typst files included directly; `#full[...]`, `#icpc[...]`,
  `#cpp[...]` and `#py[...]` read the packet's state, so one note serves
  all four packets.
- The code font size is computed in a `context` block from the column
  width so that exactly `line-limit` (64) characters fit; the linter
  enforces 64 columns, so code never wraps.
- Inconsolata is a variable font; the packets use weight 500, which
  prints darker than the regular cut at small sizes.
- The ICPC page limit is checked inside Typst (`counter(page).final()`),
  so an over-long packet fails the build instead of printing.
- The syntax theme (`typst/theme.tmTheme`) uses dark colors plus bold
  keywords and italic comments, so it still reads on a grayscale printer.

## Typing hash

Each item prints the first six hex digits of an MD5 of its code with
whitespace (C++, kactl's `cpp -dD -P -fpreprocessed | tr -d '[:space:]'`)
or comments and blank lines (PyPy, a token stream with indentation as
structure) removed. Retyped code can be checked against it; the command
is in the Contest section.

## Tests

Every snippet has a stress test per language against an independent
brute force. C++ tests build with ASan, UBSan and libstdc++ assertions
(full `_GLIBCXX_DEBUG` combined with ASan was about 300 times slower);
`<bits/stdc++.h>` is precompiled once per flag set. PyPy tests load
snippets with `testlib.load`, which executes the files and their declared
dependencies into one namespace, exactly as pasting them would. Python
`Example` fields are doctests.

## Style

Code is formatted by clang-format and ruff with 64 columns and two-space
indentation, so it fits the three-column landscape layout. C++ follows
kactl's terse conventions with clearer names; PyPy code avoids recursion
and uses flat lists for speed.
