#!/usr/bin/env python3
"""Run stress tests: C++ (sanitized g++) and PyPy, plus Python doctests.

tools/test.py                 everything
tools/test.py ds/segtree      one item (both languages)
tools/test.py ds --lang py    one section, one language
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import hashlib
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from snowi.meta import CONTENT, ROOT, TESTS, ContentError, all_sections, load_sections, test_paths  # noqa: E402

CXX_FLAGS = [
    "-std=gnu++20",
    "-O1",
    "-g",
    "-Wall",
    "-Wextra",
    "-Wshadow",
    "-Wno-unused-function",
    "-fsanitize=address,undefined",
    "-fno-sanitize-recover=all",
    # Bounds checks on operator[]; full _GLIBCXX_DEBUG with ASan is ~300x slower.
    "-D_GLIBCXX_ASSERTIONS",
]
# Binaries live under a hash of the flags, so changing flags rebuilds.
BUILD = ROOT / "build" / "tests" / hashlib.md5(" ".join(CXX_FLAGS).encode()).hexdigest()[:8]
TIMEOUT = {"cpp": 120, "py": 240, "doc": 120}


@dataclass
class Job:
    kind: str  # cpp, py, doc
    item: str
    path: Path


@dataclass
class Result:
    job: Job
    ok: bool
    seconds: float
    output: str


def find_pypy() -> str:
    for candidate in (os.environ.get("PYPY"), shutil.which("pypy3")):
        if candidate:
            return candidate
    local = Path.home() / ".local" / "bin" / "pypy3"
    if local.exists():
        return str(local)
    sys.exit("pypy3 not found: install it or set PYPY=/path/to/pypy3")


def item_of(path: Path) -> str:
    rel = path.relative_to(TESTS)
    name = rel.name
    if name.endswith(".test.cpp"):
        name = name[: -len(".test.cpp")]
    else:
        name = name[len("test_") : -len(".py")].replace("_", "-")
    return f"{rel.parent.as_posix()}/{name}"


def has_doctest(source: Path) -> bool:
    parts = source.read_text().split('"""', 2)
    return len(parts) == 3 and ">>>" in parts[1]


def discover(selector: str | None, lang: str | None) -> list[Job]:
    jobs = []
    if lang in (None, "cpp"):
        jobs += [Job("cpp", item_of(p), p) for p in sorted(TESTS.rglob("*.test.cpp"))]
    if lang in (None, "py"):
        jobs += [Job("py", item_of(p), p) for p in sorted(TESTS.rglob("test_*.py"))]
        for source in sorted(CONTENT.rglob("*.py")):
            if has_doctest(source):
                item = source.relative_to(CONTENT).with_suffix("").as_posix()
                jobs.append(Job("doc", item, source))
    if selector:
        jobs = [j for j in jobs if j.item == selector or j.item.startswith(selector + "/")]
    return jobs


def deps_digest(depfile: Path) -> str | None:
    """MD5 over the contents of every file the binary was built from."""
    if not depfile.exists():
        return None
    deps = depfile.read_text().replace("\\\n", " ").split(":", 1)[1].split()
    digest = hashlib.md5(" ".join(CXX_FLAGS).encode())
    for dep in deps:
        path = Path(dep)
        if not path.exists():
            return None
        digest.update(path.read_bytes())
    return digest.hexdigest()


def up_to_date(binary: Path, depfile: Path, stamp: Path) -> bool:
    """Content-based, so restoring a file with an old mtime still rebuilds."""
    if not binary.exists() or not stamp.exists():
        return False
    return deps_digest(depfile) == stamp.read_text()


def ensure_pch() -> Path:
    """Precompile <bits/stdc++.h> with the test flags; returns the header to -include."""
    header = BUILD / "pch" / "stdc.hpp"
    gch = header.with_suffix(".hpp.gch")
    if not gch.exists():
        header.parent.mkdir(parents=True, exist_ok=True)
        header.write_text("#include <bits/stdc++.h>\n")
        partial = gch.with_suffix(f".tmp{os.getpid()}")
        cmd = ["g++", *CXX_FLAGS, "-x", "c++-header", str(header), "-o", str(partial)]
        with compile_slot():
            subprocess.run(cmd, check=True)
        os.replace(partial, gch)  # atomic: parallel runs never see half a file
    return header


COMPILE_SLOTS = 1  # machine-wide: each sanitized g++ needs up to ~450 MB


@contextlib.contextmanager
def compile_slot():
    """Hold one of COMPILE_SLOTS file locks shared by every test.py process."""
    lock_dir = ROOT / "build" / "locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    handles = [open(lock_dir / f"compile-{i}.lock", "w") for i in range(COMPILE_SLOTS)]
    try:
        while True:
            for handle in handles:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    continue
                yield
                return
            time.sleep(0.2)
    finally:
        for handle in handles:
            handle.close()


def run_cpp(job: Job, pch: Path) -> Result:
    start = time.perf_counter()
    binary = BUILD / (job.item + ".bin")
    depfile, stamp = binary.with_suffix(".d"), binary.with_suffix(".sum")
    binary.parent.mkdir(parents=True, exist_ok=True)
    if not up_to_date(binary, depfile, stamp):
        cmd = ["g++", *CXX_FLAGS, "-include", str(pch), "-MMD", "-MF", str(depfile)]
        cmd += [str(job.path), "-o", str(binary)]
        with compile_slot():
            build = subprocess.run(cmd, capture_output=True, text=True)
        if build.returncode != 0 or "warning:" in build.stderr:
            binary.unlink(missing_ok=True)
            return Result(job, False, time.perf_counter() - start, build.stderr)
        stamp.write_text(deps_digest(depfile) or "")
    # Fake stacks for use-after-return make small calls ~1000x slower here.
    asan = "detect_leaks=0:detect_stack_use_after_return=0"
    env = os.environ | {"ASAN_OPTIONS": asan, "UBSAN_OPTIONS": "print_stacktrace=1"}
    return execute(job, [str(binary)], env, start)


def run_py(job: Job, pypy: str) -> Result:
    start = time.perf_counter()
    env = os.environ | {"PYTHONPATH": str(TESTS)}
    return execute(job, [pypy, str(job.path)], env, start)


def run_doc(job: Job, pypy: str) -> Result:
    start = time.perf_counter()
    env = os.environ | {"PYTHONPATH": str(TESTS)}
    return execute(job, [pypy, str(TESTS / "run_doctests.py"), job.item], env, start)


MEMORY_LIMIT = "1536M"  # per test; a runaway (e.g. a mutation test) dies alone


def memory_capped(cmd: list[str]) -> list[str]:
    """Run cmd in a transient systemd scope with a hard memory cap, if available."""
    if shutil.which("systemd-run") is None:
        return cmd
    scope = ["systemd-run", "--user", "--scope", "--quiet"]
    return scope + ["-p", f"MemoryMax={MEMORY_LIMIT}", "-p", "MemorySwapMax=0", *cmd]


def execute(job: Job, cmd: list[str], env: dict, start: float) -> Result:
    cmd = memory_capped(cmd)
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, env=env, cwd=ROOT, timeout=TIMEOUT[job.kind]
        )
        output = (proc.stdout + proc.stderr).strip()
        ok = proc.returncode == 0
        if proc.returncode in (-9, 137):
            output += f"\nkilled: over the {MEMORY_LIMIT} memory limit (or SIGKILL)"
    except subprocess.TimeoutExpired:
        output, ok = f"timeout after {TIMEOUT[job.kind]}s", False
    return Result(job, ok, time.perf_counter() - start, output)


def missing_tests(selector: str | None) -> list[str]:
    try:
        sections = load_sections([selector.split("/")[0]]) if selector else all_sections()
    except ContentError as err:
        return [f"(cannot check: {err})"]
    missing = []
    for section in sections:
        for entry in section.entries:
            if entry.kind != "item":
                continue
            if selector and not (entry.id == selector or entry.id.startswith(selector + "/")):
                continue
            paths = test_paths(entry)
            for lang in entry.snippets:
                if not paths[lang].exists():
                    missing.append(f"{entry.id} ({lang})")
    return missing


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("selector", nargs="?", help="item id or section")
    parser.add_argument("--lang", choices=["cpp", "py"])
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 2)
    parser.add_argument("-v", "--verbose", action="store_true", help="show output of passing tests")
    args = parser.parse_args()

    pypy = find_pypy() if args.lang != "cpp" else ""
    jobs = discover(args.selector, args.lang)
    if not jobs:
        sys.exit(f"no tests match {args.selector!r}")

    pch = ensure_pch() if any(j.kind == "cpp" for j in jobs) else Path()

    def run(job: Job) -> Result:
        if job.kind == "cpp":
            return run_cpp(job, pch)
        return run_py(job, pypy) if job.kind == "py" else run_doc(job, pypy)

    failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for result in pool.map(run, jobs):
            status = "PASS" if result.ok else "FAIL"
            print(
                f"{status} {result.job.kind:3} {result.job.item:32} {result.seconds:6.1f}s",
                flush=True,
            )
            if not result.ok or args.verbose:
                tail = "\n".join(result.output.splitlines()[-40:])
                print("    " + tail.replace("\n", "\n    "))
            failed += not result.ok
    missing = missing_tests(args.selector)
    for name in missing:
        print(f"MISSING test for {name}")
    print(f"{len(jobs) - failed}/{len(jobs)} passed, {len(missing)} missing")
    sys.exit(1 if failed or missing else 0)


if __name__ == "__main__":
    main()
