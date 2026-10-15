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
    base = {"id": entry.id, "kind": entry.kind, "title": entry.title, "icpc": list(entry.icpc)}
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


def write_manifest(out: Path, allow_missing: bool, only: set[str] | None) -> dict:
    config = tomllib.loads((ROOT / "config" / "packets.toml").read_text())
    sections = []
    for section in all_sections():
        if only and section.name not in only:
            continue
        entries = [entry_json(e, allow_missing) for e in section.entries]
        sections.append(
            {
                "id": section.name,
                "title": section.title,
                "entries": [e for e in entries if e is not None],
            }
        )
    manifest = {"config": config, "sections": sections, "langs": LANG_NAMES}
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def manifest_input(out: Path) -> list[str]:
    return ["--input", "manifest=/" + (out / "manifest.json").relative_to(ROOT).as_posix()]


def typst_base(typst: str) -> list[str]:
    # Tags (accessibility structure) are ~75% of the file size of these
    # code-dense packets and add nothing to a printed notebook.
    root = ["--root", str(ROOT), "--font-path", str(ROOT / "fonts")]
    return [typst, "compile", *root, "--no-pdf-tags"]


def compile_bundle(typst: str, out: Path) -> None:
    bundle = out / "bundle"
    cmd = (
        typst_base(typst)
        + manifest_input(out)
        + [
            "--features",
            "bundle",
            "--format",
            "bundle",
            str(ROOT / "typst" / "bundle.typ"),
            str(bundle),
        ]
    )
    run(cmd)
    (out / "pdf").mkdir(parents=True, exist_ok=True)
    for pdf in bundle.glob("*.pdf"):
        shutil.copy2(pdf, out / "pdf" / pdf.name)


def compile_single(typst: str, out: Path, variant: str, lang: str, png: bool = False) -> None:
    (out / "pdf").mkdir(parents=True, exist_ok=True)
    name = f"snowi-{variant}-{lang}" + ("-{p}.png" if png else ".pdf")
    target = out / ("png" if png else "pdf") / name
    target.parent.mkdir(parents=True, exist_ok=True)
    cmd = (
        typst_base(typst)
        + manifest_input(out)
        + (["--ppi", "110"] if png else [])
        + [
            "--input",
            f"variant={variant}",
            "--input",
            f"lang={lang}",
            str(ROOT / "typst" / "single.typ"),
            str(target),
        ]
    )
    run(cmd)


def page_report(typst: str, out: Path, variant: str, lang: str, manifest: dict) -> None:
    """Print how many columns each entry takes in one packet."""
    expr = (
        "(entries: query(<snowi-pos>).map(it => it.value), geo: query(<snowi-pages>).first().value)"
    )
    cmd = [typst, "eval", "--root", str(ROOT), "--font-path", str(ROOT / "fonts")]
    cmd += manifest_input(out) + ["--input", f"variant={variant}", "--input", f"lang={lang}"]
    cmd += [expr, "--in", str(ROOT / "typst" / "single.typ")]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"typst eval failed:\n{result.stderr}")
    data = json.loads(result.stdout)
    geo = data["geo"]
    cols = geo["columns"]
    col_w = (geo["width"] - 2 * geo["margin-x"] - (cols - 1) * geo["gutter"]) / cols
    col_h = geo["height"] - geo["top"] - geo["bottom"]

    def position(entry: dict) -> float:
        col = min(cols - 1, int((entry["x"] - geo["margin-x"]) / (col_w + geo["gutter"]) + 0.5))
        return (entry["page"] - 1) * cols + col + (entry["y"] - geo["top"]) / col_h

    starts = [(e["id"], position(e)) for e in data["entries"]]
    ends = [p for _, p in starts[1:]] + [geo["pages"] * cols]
    titles = {e["id"]: e["title"] for s in manifest["sections"] for e in s["entries"]}
    per_section: dict[str, float] = {}
    for (entry_id, start), end in zip(starts, ends, strict=True):
        size = end - start
        per_section[entry_id.split("/")[0]] = per_section.get(entry_id.split("/")[0], 0) + size
        print(
            f"{size:5.2f} col  p{int(start // cols) + 1:<3} {entry_id:28} {titles.get(entry_id, '')}"
        )
    print("-- per section (columns):")
    for section, size in per_section.items():
        print(f"{size:6.2f}  {section}")
    print(f"total: {geo['pages']} pages, {geo['pages'] * cols} columns")


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
    parser.add_argument("--out", type=Path, default=BUILD, help="output directory under the repo")
    parser.add_argument("--only", help="comma-separated sections (preview builds)")
    parser.add_argument(
        "--png", metavar="VARIANT-LANG", help="render PNG pages of one packet, e.g. full-cpp"
    )
    parser.add_argument(
        "--report", metavar="VARIANT-LANG", help="columns used per entry, e.g. icpc-py"
    )
    args = parser.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(ROOT):
        sys.exit("--out must be inside the repository (Typst reads files under --root)")
    only = set(args.only.split(",")) if args.only else None
    try:
        manifest = write_manifest(out, args.allow_missing, only)
    except ContentError as err:
        sys.exit(f"content error: {err}")
    except (OSError, subprocess.CalledProcessError) as err:
        cpp = os.environ.get("CPP", "cpp")
        sys.exit(f"typing hash failed: the C preprocessor CPP={cpp!r} did not run ({err})")
    if args.manifest_only:
        return
    typst = find_typst()
    if args.report:
        variant, lang = args.report.split("-")
        page_report(typst, out, variant, lang, manifest)
        return
    if args.png:
        variant, lang = args.png.split("-")
        compile_single(typst, out, variant, lang, png=True)
        for page in sorted((out / "png").glob(f"snowi-{args.png}-*.png")):
            print(page.relative_to(ROOT))
        return
    if args.no_bundle:
        jobs = [(v, lang) for v in VARIANTS for lang in LANG_NAMES]
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(lambda job: compile_single(typst, out, *job), jobs))
    else:
        compile_bundle(typst, out)
    for pdf in sorted((out / "pdf").glob("*.pdf")):
        print(f"{pdf.relative_to(ROOT)}  {pdf.stat().st_size // 1024} KiB")


if __name__ == "__main__":
    main()
