"""Shared helpers for PyPy stress tests.

lib = load("ds/segtree") executes the snippet, after its `Uses:`
dependencies, into one namespace, the same as pasting them into a file.
"""

import os
import random
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
CONTENT = ROOT / "content"

rng = random.Random(int(os.environ.get("SNOWI_SEED", "20261015")))


def rnd(lo, hi):
  """Uniform integer in [lo, hi]."""
  return rng.randint(lo, hi)


def rnd_list(n, lo, hi):
  return [rng.randint(lo, hi) for _ in range(n)]


def iters(n):
  """Scale iteration counts: SNOWI_SCALE=0.1 for a quick run."""
  return max(1, int(n * float(os.environ.get("SNOWI_SCALE", "1"))))


def header_uses(source):
  match = re.search(r"^Uses: (.*)$", source, re.M)
  if not match:
    return []
  return [u.strip() for u in match.group(1).split(",") if u.strip()]


def load(*item_ids):
  """Exec snippets (dependencies first) into one shared namespace."""
  namespace = {"__name__": "snowi_snippet"}
  done = set()

  def visit(item_id):
    if item_id in done:
      return
    done.add(item_id)
    path = CONTENT / (item_id + ".py")
    source = path.read_text()
    for dep in header_uses(source):
      visit(dep)
    exec(compile(source, str(path), "exec"), namespace)

  for item_id in item_ids:
    visit(item_id)
  return SimpleNamespace(**namespace)


def check(cond, msg=""):
  if not cond:
    raise AssertionError(msg or "check failed")


def check_eq(got, want, msg=""):
  if got != want:
    raise AssertionError(f"{msg} got {got!r}, want {want!r}")


class Timer:
  """with Timer("label"): ...  prints elapsed seconds to stderr."""

  def __init__(self, label):
    self.label = label

  def __enter__(self):
    self.start = time.perf_counter()
    return self

  def __exit__(self, *exc):
    self.seconds = time.perf_counter() - self.start
    print(f"  {self.label}: {self.seconds:.2f}s", file=sys.stderr)
