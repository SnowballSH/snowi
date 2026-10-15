"""Runs a solution built on the snippet in a child process and
plays the judge through two pipes, as a real interactor does.

Protocol: T, then per test n. Hidden x in [1, n]; "? y" replies 1
iff x <= y; "? y z" replies "y + z" and "y * z" on one line;
"! x d" is right if d is the first number minus the second (so the
order of the tokens matters), and the judge confirms with 1.
"""

import io
import os
import re
import select
import subprocess
import sys

from testlib import CONTENT, Timer, check, check_eq, iters, rnd

# READ_T is how the solution reads T: with reply(), or with input()
# as contestants often do. The judge sends T and the first n before
# it reads anything, so a reply() that bypasses input()'s buffer
# would wait forever.
SOLUTION = """
from testlib import load

lib = load("misc/interactive")
for _ in range(READ_T):
  n = lib.reply()
  lo, hi = 1, n
  while lo < hi:
    mid = (lo + hi) // 2
    if lib.ask(mid):
      hi = mid
    else:
      lo = mid + 1
  s, p = lib.ask(lo, 7), lib.reply()
  lib.answer(lo, s - p)
  lib.reply()
"""
WAIT = 10.0  # seconds; a missing flush shows up as a timeout


class Judge:
  def __init__(self, read_t):
    self.proc = subprocess.Popen(
      [sys.executable, "-c", SOLUTION.replace("READ_T", read_t)],
      stdin=subprocess.PIPE,
      stdout=subprocess.PIPE,
    )
    self.out = self.proc.stdout.fileno()
    self.buf = b""

  def send(self, line):
    self.proc.stdin.write(line.encode() + b"\n")
    self.proc.stdin.flush()

  def line(self):
    """The next line the child wrote, or None on timeout or EOF."""
    while b"\n" not in self.buf:
      ready, _, _ = select.select([self.out], [], [], WAIT)
      chunk = os.read(self.out, 1 << 16) if ready else b""
      if not chunk:
        return None
      self.buf += chunk
    line, self.buf = self.buf.split(b"\n", 1)
    return line.decode()

  def finish(self):
    """Exit code of the child."""
    if not self.proc.stdin.closed:
      self.proc.stdin.close()
    try:
      return self.proc.wait(WAIT)
    except subprocess.TimeoutExpired:
      return self.fail("child did not exit")

  def fail(self, msg):
    self.proc.kill()
    self.proc.wait()
    return msg


LINE = re.compile(r"([?!])((?: -?(?:0|[1-9][0-9]*)){1,2})")


def parse(line):
  """("?" or "!", [numbers]) for a well-formed line, else (None, [])."""
  m = LINE.fullmatch(line)
  return (m[1], [int(v) for v in m[2].split()]) if m else (None, [])


def session(ns, reject, hang_up, stats, read_t="lib.reply()"):
  """Plays one session; returns an error message, "" if all went
  well. With reject = k >= 0, the judge answers the k-th line the
  child sends with -1 (or, with hang_up, closes the pipe), and the
  child must exit with status 0 without writing more."""
  j = Judge(read_t)
  lines = 0
  j.send(str(len(ns)))
  for n in ns:
    x, used = rnd(1, n), 0
    j.send(str(n))
    while True:
      got = j.line()
      if got is None:
        return j.fail("no flushed line from the child")
      if lines == reject:
        if not hang_up:
          j.send("-1")
        j.proc.stdin.close()
        more = j.line()
        if more is not None:
          return j.fail("child wrote after -1: " + more)
        return "" if j.finish() == 0 else "bad exit after -1"
      lines += 1
      kind, nums = parse(got)
      if kind == "?" and len(nums) == 1:
        y = nums[0]
        if not 1 <= y <= n:
          return j.fail("query out of range: " + got)
        used += 1
        if used > (n - 1).bit_length():
          return j.fail("too many queries")
        j.send("1" if x <= y else "0")
      elif kind == "?" and len(nums) == 2:
        y, z = nums
        j.send(f"{y + z} {y * z}")
      elif kind == "!" and len(nums) == 2:
        if nums != [x, (x + 7) - 7 * x]:
          return j.fail("wrong answer: " + got)
        j.send("1")
        break
      else:
        return j.fail("malformed line: " + repr(got))
    stats[0] += used + 2
  more = j.line()
  if more is not None:
    return j.fail("extra output: " + more)
  return "" if j.finish() == 0 else "bad exit status"


def word_replies():
  """The Tweak for word replies (return r, not int(r)) keeps the -1
  guard: reply() returns the tokens, then exits on -1 and on EOF."""
  src = (CONTENT / "misc/interactive.py").read_text()
  check_eq(src.count("return int(r)"), 1, "Tweak target")
  for tail in ("-1\nYES\n", ""):
    namespace = {"__name__": "snowi_snippet"}
    exec(src.replace("return int(r)", "return r"), namespace)
    stdin, sys.stdin = sys.stdin, io.StringIO("YES <\n" + tail)
    try:
      check_eq([namespace["reply"]() for _ in range(2)], ["YES", "<"])
      namespace["reply"]()
      check(False, "reply() did not exit on " + repr(tail))
    except SystemExit as e:
      check_eq(e.code, 0)
    finally:
      sys.stdin = stdin


word_replies()
stats = [0]
for _ in range(iters(40)):
  ns = [rnd(1, 20) if rnd(0, 1) else rnd(1, 10**18) for _ in range(rnd(1, 4))]
  reject = -1 if rnd(0, 2) else rnd(0, 6)
  read_t = "lib.reply()" if rnd(0, 1) else "int(input())"
  check_eq(session(ns, reject, rnd(0, 1), stats, read_t), "")

with Timer("3000 tests, about 1e5 round trips"):
  stats[0] = 0
  check_eq(session([10**9] * 3000, -1, False, stats), "")
print(f"  {stats[0]} round trips", file=sys.stderr)
print("OK misc/interactive")
