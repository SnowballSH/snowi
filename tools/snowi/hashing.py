"""Typing-check hashes, identical to the commands printed in the packet.

C++:  cpp -dD -P -fpreprocessed | tr -d '[:space:]' | md5sum | cut -c-6
PyPy: tokens without comments; indentation kept as structure markers.

The C++ hash is computed here without a preprocessor, so the build runs
where GCC's cpp is missing (macOS ships Clang's, which has no
-fpreprocessed); `tools/lint.py` checks it against GCC's cpp when present.
"""

from __future__ import annotations

import hashlib
import io
import re
import tokenize

STRUCTURE = {tokenize.INDENT: "{", tokenize.DEDENT: "}", tokenize.NEWLINE: ";"}
SKIPPED = {tokenize.COMMENT, tokenize.NL, tokenize.ENCODING, tokenize.ENDMARKER}
ASCII_SPACE = re.compile(r"[ \t\n\v\f\r]")  # what tr -d '[:space:]' removes
# Lines cpp drops: the null directive anywhere, and #pragma once (with
# anything after it) when its # is in column 0.
DROPPED_LINES = re.compile(r"^(?:[ \t]*#[ \t]*|#[ \t]*pragma[ \t]+once(?![A-Za-z0-9_]).*)$", re.M)


def strip_cpp_comments(code: str) -> str:
    """The text `cpp -dD -P -fpreprocessed` keeps: comments removed, nothing
    expanded, directives passed through. Raises ValueError where cpp fails."""
    out = []
    i, n = 0, len(code)
    while i < n:
        if code.startswith("//", i):
            end = code.find("\n", i)
            i = n if end < 0 else end
        elif code.startswith("/*", i):
            end = code.find("*/", i + 2)
            if end < 0:
                raise ValueError("unterminated comment")
            out.append(" ")
            i = end + 2
        elif code[i] in "\"'":
            # A literal ends at its closing quote, or unterminated at a newline;
            # a backslash escapes the next character, but never a newline
            # (-fpreprocessed does no line splicing).
            j = i + 1
            while j < n and code[j] not in (code[i], "\n"):
                j += 2 if code.startswith("\\", j) and code[j + 1 : j + 2] != "\n" else 1
            end = j + 1 if j < n and code[j] == code[i] else min(j, n)
            out.append(code[i:end])
            i = end
        else:
            out.append(code[i])
            i += 1
    return DROPPED_LINES.sub("", "".join(out))


def cpp_hash(code: str) -> str:
    compact = ASCII_SPACE.sub("", strip_cpp_comments(code))
    return hashlib.md5(compact.encode()).hexdigest()[:6]


def py_hash(code: str) -> str:
    tokens = tokenize.tokenize(io.BytesIO(code.encode()).readline)
    parts = [STRUCTURE.get(tok.type, tok.string) for tok in tokens if tok.type not in SKIPPED]
    return hashlib.md5("".join(parts).encode()).hexdigest()[:6]


def code_hash(lang: str, code: str) -> str:
    return cpp_hash(code) if lang == "cpp" else py_hash(code)
