# snowi

A competitive programming library and notebook for ICPC and UCup, in C++
and PyPy, in the spirit of kactl. Every snippet is stress-tested, and the
printed packets are generated from the sources with Typst.

## Packets

`make pdf` writes four PDFs to `build/pdf/`:

| File | Contents |
|---|---|
| `snowi-icpc-cpp.pdf` | ICPC team reference document, C++ (at most 25 pages) |
| `snowi-icpc-py.pdf` | ICPC team reference document, PyPy (at most 25 pages) |
| `snowi-full-cpp.pdf` | Everything: all snippets, examples, notes, guides (C++) |
| `snowi-full-py.pdf` | Everything (PyPy) |

The ICPC packets follow the North America / World Finals rules: letter
paper, single-sided, at most 25 pages, university and team name upper
left, page number upper right. Set the names in `config/packets.toml`. The
build fails if an ICPC packet goes over the page limit.

## Requirements

- Typst 0.15 or newer (`typst`), g++ 13, PyPy 7.3.15 (`pypy3`, Python 3.9;
  the ICPC World Finals version), Python 3.11+ for the tools.
- For `make lint`: `clang-format` and `ruff` (`uv tool install clang-format ruff`).
- Tools are looked up on `PATH`, then in `~/.local/bin`; override with
  `TYPST=...` and `PYPY=...`.

## Commands

```
make pdf                      build the four PDFs (one Typst bundle compile)
make test                     all stress tests: C++ (ASan/UBSan) and PyPy
make test ITEM=ds/segtree     one item; ITEM=ds for a section; ONLY=cpp|py
make lint                     headers, 64-column code, clang-format, ruff
make format                   apply clang-format and ruff format
make check                    lint, test and pdf
```

## Layout

```
config/packets.toml   section order, team names, layout, page limit
content/<section>/    snippets (.hpp, .py), notes (.typ), section.toml
tests/<section>/      stress tests (<item>.test.cpp, test_<item>.py)
tools/                build.py, test.py, lint.py and the snowi package
typst/                packet layout, bundle entry point, print theme
fonts/                Inconsolata (OFL), the code font
```

`CONTRIBUTING.md` specifies the snippet header format, the code style
and the test conventions. `docs/DESIGN.md` explains the design choices.

## Using a snippet

C++ snippets assume `content/contest/template.cpp` (kactl's template:
`bits/stdc++.h`, `rep`, `all`, `sz`, `ll`, `pii`, `vi`). Paste the code
below the template. PyPy snippets are self-contained; paste them above
`main()` in `content/contest/template.py`. Each item lists its
dependencies under "Uses".

## Licence

The code is original, written for this repository, with ideas credited in
each snippet's `Source` field (mostly kactl, the AtCoder Library and maspy's
library, all CC0). Inconsolata is under the SIL Open Font License.
