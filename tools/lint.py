#!/usr/bin/env python3
"""Check content: headers, 64-column code, formatting, tests present.

tools/lint.py            everything
tools/lint.py ds         one section or item
tools/lint.py --fix      apply clang-format and ruff format first
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from snowi.meta import (  # noqa: E402
    CODE_FIELDS,
    CONTENT,
    LINE_LIMIT,
    ROOT,
    TESTS,
    ContentError,
    all_sections,
    test_paths,
)


def tool(name: str) -> str:
    found = shutil.which(name) or str(Path.home() / ".local" / "bin" / name)
    if not Path(found).exists():
        sys.exit(f"{name} not found (uv tool install {name})")
    return found


def matches(item_id: str, selector: str | None) -> bool:
    return not selector or item_id == selector or item_id.startswith(selector + "/")


def check_lines(label: str, text: str, problems: list[str]) -> None:
    for number, line in enumerate(text.split("\n"), 1):
        if len(line) > LINE_LIMIT:
            problems.append(f"{label}:{number}: {len(line)} > {LINE_LIMIT} columns")
        if "\t" in line:
            problems.append(f"{label}:{number}: tab character")
        if line != line.rstrip():
            problems.append(f"{label}:{number}: trailing whitespace")


def check_content(selector: str | None) -> list[str]:
    problems: list[str] = []
    try:
        sections = all_sections()
    except ContentError as err:
        return [str(err)]
    listed: set[Path] = set()
    for section in sections:
        for entry in section.entries:
            if entry.kind == "note":
                if entry.note_path:
                    listed.add(entry.note_path)
                elif matches(entry.id, selector):
                    problems.append(f"{entry.id}: missing note file")
                continue
            for lang in entry.langs:
                snippet = entry.snippets.get(lang)
                if not matches(entry.id, selector):
                    if snippet:
                        listed.add(snippet.path)
                    continue
                if snippet is None:
                    problems.append(f"{entry.id}: missing {lang} source")
                    continue
                listed.add(snippet.path)
                check_lines(snippet.rel, snippet.code, problems)
                for name in CODE_FIELDS:
                    if name in snippet.fields:
                        check_lines(f"{snippet.rel} {name}", snippet.fields[name], problems)
                if not test_paths(entry)[lang].exists():
                    problems.append(
                        f"{entry.id}: missing {lang} test {test_paths(entry)[lang].relative_to(ROOT)}"
                    )
    if not selector:
        for path in sorted(CONTENT.rglob("*")):
            if path.suffix in (".hpp", ".cpp", ".py", ".typ") and path not in listed:
                problems.append(f"{path.relative_to(ROOT)}: not listed in section.toml")
    return problems


def scoped_files(selector: str | None, suffixes: tuple[str, ...]) -> list[str]:
    files = []
    for base in (CONTENT, TESTS):
        for path in base.rglob("*"):
            if path.suffix not in suffixes or not path.is_file():
                continue
            rel = path.relative_to(base).as_posix()
            if selector and not (
                rel.startswith(selector + "/") or selector in rel.replace("_", "-")
            ):
                continue
            files.append(str(path))
    return sorted(files)


def run(cmd: list[str], problems: list[str], label: str) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if proc.returncode != 0:
        problems.append(f"{label}:\n{(proc.stdout + proc.stderr).strip()}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("selector", nargs="?")
    parser.add_argument("--fix", action="store_true")
    args = parser.parse_args()

    cpp_files = scoped_files(args.selector, (".hpp", ".cpp"))
    py_files = scoped_files(args.selector, (".py",))
    if not args.selector:
        py_files += [str(p) for p in (ROOT / "tools").rglob("*.py")]
    clang, ruff = tool("clang-format"), tool("ruff")
    if args.fix:
        if cpp_files:
            subprocess.run([clang, "-i", *cpp_files], check=True)
        if py_files:
            subprocess.run([ruff, "format", "-q", *py_files], check=True)
            subprocess.run([ruff, "check", "-q", "--fix", *py_files])

    problems = check_content(args.selector)
    if cpp_files:
        run([clang, "--dry-run", "-Werror", *cpp_files], problems, "clang-format")
    if py_files:
        run([ruff, "format", "--check", "-q", *py_files], problems, "ruff format")
        run([ruff, "check", "-q", *py_files], problems, "ruff check")
    for problem in problems:
        print(problem)
    print(f"lint: {len(problems)} problem(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
