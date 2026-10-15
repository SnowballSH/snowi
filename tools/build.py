#!/usr/bin/env python3
"""Generate build/manifest.json from content/ and compile the packets.

One Typst bundle compile emits all four PDFs (ICPC/full x C++/PyPy) into
build/pdf/. Use --no-bundle to compile them one by one instead.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from snowi.hashing import code_hash  # noqa: E402
from snowi.meta import ROOT, ContentError, Entry, all_sections  # noqa: E402

BUILD = ROOT / "build"
PDF_DIR = BUILD / "pdf"
VARIANTS = ("icpc", "full")
LANG_NAMES = {"cpp": "C++", "py": "PyPy"}


def find_typst() -> str:
    for candidate in (os.environ.get("TYPST"), shutil.which("typst")):
        if candidate:
            return candidate
    local = Path.home() / ".local" / "bin" / "typst"
    if local.exists():
        return str(local)
    sys.exit("typst not found: install it or set TYPST=/path/to/typst")


def entry_json(entry: Entry, allow_missing: bool) -> dict | None:
    base = {"id": entry.id, "kind": entry.kind, "title": entry.title, "icpc": entry.icpc}
    if entry.kind == "note":
        if entry.note_path is None:
            if allow_missing:
                return None
            raise ContentError(f"{entry.id}: missing {entry.id}.typ")
        return base | {"path": "/" + entry.note_path.relative_to(ROOT).as_posix()}
    missing = [lang for lang in entry.langs if lang not in entry.snippets]
    if missing and not allow_missing:
        raise ContentError(f"{entry.id}: missing {', '.join(missing)} source")
    code = {}
    for lang, snippet in entry.snippets.items():
        code[lang] = {
            "file": snippet.path.name,
            "fields": snippet.fields,
            "code": snippet.code,
            "hash": code_hash(lang, snippet.code),
            "uses": snippet.uses,
            "lines": snippet.code.count("\n") + 1,
        }
    if not code:
        return None
    return base | {"code": code}


def write_manifest(allow_missing: bool) -> dict:
    config = tomllib.loads((ROOT / "config" / "packets.toml").read_text())
    sections = []
    for section in all_sections():
        entries = [entry_json(e, allow_missing) for e in section.entries]
        sections.append(
            {
                "id": section.name,
                "title": section.title,
                "entries": [e for e in entries if e is not None],
            }
        )
    manifest = {"config": config, "sections": sections, "langs": LANG_NAMES}
    BUILD.mkdir(exist_ok=True)
    (BUILD / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def typst_base(typst: str) -> list[str]:
    return [typst, "compile", "--root", str(ROOT), "--font-path", str(ROOT / "fonts")]


def compile_bundle(typst: str) -> None:
    out = BUILD / "bundle"
    cmd = typst_base(typst) + [
        "--features",
        "bundle",
        "--format",
        "bundle",
        str(ROOT / "typst" / "bundle.typ"),
        str(out),
    ]
    run(cmd)
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    for pdf in out.glob("*.pdf"):
        shutil.copy2(pdf, PDF_DIR / pdf.name)


def compile_single(typst: str, variant: str, lang: str) -> None:
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    name = f"snowi-{variant}-{lang}.pdf"
    cmd = typst_base(typst) + [
        "--input",
        f"variant={variant}",
        "--input",
        f"lang={lang}",
        str(ROOT / "typst" / "single.typ"),
        str(PDF_DIR / name),
    ]
    run(cmd)


def run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    noise = ("bundle export is experimental", "hint: its behaviour", "hint: do not rely")
    errors = "\n".join(
        line for line in result.stderr.splitlines() if not any(n in line for n in noise)
    ).strip()
    if result.returncode != 0:
        sys.exit(f"typst failed:\n{errors}")
    if errors:
        print(errors, file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-missing", action="store_true", help="skip unwritten entries")
    parser.add_argument("--no-bundle", action="store_true", help="compile PDFs separately")
    parser.add_argument("--manifest-only", action="store_true")
    args = parser.parse_args()
    try:
        write_manifest(args.allow_missing)
    except ContentError as err:
        sys.exit(f"content error: {err}")
    if args.manifest_only:
        return
    typst = find_typst()
    if args.no_bundle:
        jobs = [(v, lang) for v in VARIANTS for lang in LANG_NAMES]
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(lambda job: compile_single(typst, *job), jobs))
    else:
        compile_bundle(typst)
    for pdf in sorted(PDF_DIR.glob("*.pdf")):
        print(f"{pdf.relative_to(ROOT)}  {pdf.stat().st_size // 1024} KiB")


if __name__ == "__main__":
    main()
