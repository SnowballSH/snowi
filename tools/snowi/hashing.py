"""Typing-check hashes, identical to the commands printed in the packet.

C++:  cpp -dD -P -fpreprocessed | tr -d '[:space:]' | md5sum | cut -c-6
PyPy: tokens without comments; indentation kept as structure markers.
"""

from __future__ import annotations

import hashlib
import io
import subprocess
import tokenize

STRUCTURE = {tokenize.INDENT: "{", tokenize.DEDENT: "}", tokenize.NEWLINE: ";"}
SKIPPED = {tokenize.COMMENT, tokenize.NL, tokenize.ENCODING, tokenize.ENDMARKER}


def cpp_hash(code: str) -> str:
    stripped = subprocess.run(
        ["cpp", "-dD", "-P", "-fpreprocessed"],
        input=code,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    compact = "".join(stripped.split())
    return hashlib.md5(compact.encode()).hexdigest()[:6]


def py_hash(code: str) -> str:
    tokens = tokenize.tokenize(io.BytesIO(code.encode()).readline)
    parts = [STRUCTURE.get(tok.type, tok.string) for tok in tokens if tok.type not in SKIPPED]
    return hashlib.md5("".join(parts).encode()).hexdigest()[:6]


def code_hash(lang: str, code: str) -> str:
    return cpp_hash(code) if lang == "cpp" else py_hash(code)
