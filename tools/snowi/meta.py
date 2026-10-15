"""Parse content/: sections, snippet headers and notes."""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from textwrap import dedent

ROOT = Path(__file__).resolve().parents[2]
CONTENT = ROOT / "content"
TESTS = ROOT / "tests"

LANGS = ("cpp", "py")
EXTENSIONS = {"cpp": (".hpp", ".cpp"), "py": (".py",)}
LINE_LIMIT = 64

TEXT_FIELDS = (
    "Title",
    "Description",
    "Precondition",
    "Time",
    "Tweak",
    "Uses",
    "Notes",
    "Source",
    "Status",
)
CODE_FIELDS = ("Usage", "Example")
FIELDS = TEXT_FIELDS + CODE_FIELDS
REQUIRED = ("Title", "Description", "Usage", "Status")

FIELD_RE = re.compile(rf"^({'|'.join(FIELDS)}):(?: (.*))?$")
INCLUDE_RE = re.compile(r'^#include "([^"]+)"\s*$')
NOQA_RE = re.compile(r"\s*(# noqa\b.*|// NOLINT\b.*)$", re.M)


class ContentError(Exception):
    pass


@dataclass
class Snippet:
    lang: str
    path: Path
    fields: dict[str, str]
    code: str
    uses: list[str]

    @property
    def rel(self) -> str:
        return self.path.relative_to(ROOT).as_posix()


@dataclass
class Entry:
    section: str
    name: str
    kind: str  # "item" or "note"
    icpc: bool
    title: str
    langs: tuple[str, ...] = LANGS
    snippets: dict[str, Snippet] = field(default_factory=dict)
    note_path: Path | None = None

    @property
    def id(self) -> str:
        return f"{self.section}/{self.name}"


@dataclass
class Section:
    name: str
    title: str
    entries: list[Entry]


def _split_header(lines: list[str], lang: str, path: Path) -> tuple[list[str], int]:
    """Return the header's inner lines and the index of the first body line."""
    if lang == "cpp":
        if not lines or lines[0].strip() != "/**":
            raise ContentError(f"{path}: must start with a '/**' header")
        for i in range(1, len(lines)):
            if lines[i].strip() == "*/":
                inner = [re.sub(r"^ ?\* ?", "", ln, count=1) for ln in lines[1:i]]
                return inner, i + 1
        raise ContentError(f"{path}: header is not closed with ' */'")
    if not lines or lines[0].strip() != '"""':
        raise ContentError(f'{path}: must start with a \'"""\' docstring header')
    for i in range(1, len(lines)):
        if lines[i].strip() == '"""':
            return lines[1:i], i + 1
    raise ContentError(f'{path}: docstring header is not closed with \'"""\'')


def parse_fields(inner: list[str], path: Path) -> dict[str, str]:
    raw: dict[str, list[str]] = {}
    current: str | None = None
    for line in inner:
        match = FIELD_RE.match(line)
        if match:
            current = match.group(1)
            if current in raw:
                raise ContentError(f"{path}: duplicate field {current}")
            raw[current] = [match.group(2) or ""]
        elif current is None:
            if line.strip():
                raise ContentError(f"{path}: text before first field: {line!r}")
        else:
            raw[current].append(line)
    fields: dict[str, str] = {}
    for name, chunk in raw.items():
        if name in CODE_FIELDS:
            first, rest = chunk[0], dedent("\n".join(chunk[1:])).split("\n")
            body = ([first] if first.strip() else []) + rest
            fields[name] = "\n".join(body).strip("\n")
        else:
            fields[name] = " ".join(s.strip() for s in chunk if s.strip())
    missing = [f for f in REQUIRED if not fields.get(f)]
    if missing:
        raise ContentError(f"{path}: missing fields {', '.join(missing)}")
    return fields


def _resolve_include(path: Path, target: str) -> str:
    resolved = (path.parent / target).resolve()
    try:
        rel = resolved.relative_to(CONTENT)
    except ValueError as exc:
        raise ContentError(f"{path}: include outside content/: {target}") from exc
    return rel.with_suffix("").as_posix()


def parse_snippet(path: Path, lang: str) -> Snippet:
    lines = path.read_text().split("\n")
    inner, start = _split_header(lines, lang, path)
    fields = parse_fields(inner, path)
    body = lines[start:]
    uses = [u.strip() for u in fields.get("Uses", "").split(",") if u.strip()]
    if lang == "cpp":
        kept = []
        for line in body:
            match = INCLUDE_RE.match(line)
            if match:
                uses.append(_resolve_include(path, match.group(1)))
            elif line.strip() != "#pragma once":
                kept.append(line)
        body = kept
    code = printable("\n".join(body))
    if not code:
        raise ContentError(f"{path}: no code after the header")
    return Snippet(lang, path, fields, code, list(dict.fromkeys(uses)))


def printable(code: str) -> str:
    """Code as printed: lint pragmas dropped, blank-line runs collapsed."""
    code = NOQA_RE.sub("", code)
    code = re.sub(r"\n\s*\n(\s*\n)+", "\n\n", code)
    return code.strip("\n")


def _snippet_path(section: str, name: str, lang: str) -> Path | None:
    for ext in EXTENSIONS[lang]:
        candidate = CONTENT / section / f"{name}{ext}"
        if candidate.exists():
            return candidate
    return None


def load_sections(order: list[str]) -> list[Section]:
    sections = []
    for name in order:
        toml_path = CONTENT / name / "section.toml"
        if not toml_path.exists():
            raise ContentError(f"missing {toml_path.relative_to(ROOT)}")
        data = tomllib.loads(toml_path.read_text())
        entries = [_load_entry(name, raw, toml_path) for raw in data.get("entry", [])]
        sections.append(Section(name, data["title"], entries))
    return sections


def _load_entry(section: str, raw: dict, toml_path: Path) -> Entry:
    icpc = bool(raw.get("icpc", False))
    if "note" in raw:
        path = CONTENT / section / f"{raw['note']}.typ"
        if "title" not in raw:
            raise ContentError(f"{toml_path}: note {raw['note']} needs a title")
        entry = Entry(section, raw["note"], "note", icpc, raw["title"])
        entry.note_path = path if path.exists() else None
        return entry
    name = raw["item"]
    langs = tuple(raw.get("langs", LANGS))
    entry = Entry(section, name, "item", icpc, raw.get("title", name), langs)
    for lang in langs:
        path = _snippet_path(section, name, lang)
        if path is not None:
            entry.snippets[lang] = parse_snippet(path, lang)
    titles = {s.fields["Title"] for s in entry.snippets.values()}
    if len(titles) > 1:
        raise ContentError(f"{entry.id}: titles differ between languages: {titles}")
    if titles:
        entry.title = titles.pop()
    return entry


def test_paths(entry: Entry) -> dict[str, Path]:
    base = TESTS / entry.section
    return {
        "cpp": base / f"{entry.name}.test.cpp",
        "py": base / f"test_{entry.name.replace('-', '_')}.py",
    }


def section_order() -> list[str]:
    config = tomllib.loads((ROOT / "config" / "packets.toml").read_text())
    return list(config["sections"])


def all_sections() -> list[Section]:
    return load_sections(section_order())
