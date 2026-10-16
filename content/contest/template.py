"""
Title: Contest template
Description: `main()` makes locals fast in PyPy.
Usage:
  x = nxt()  # next int token
Status: run on sample input by tests/contest/test_template.py
"""

import sys
from bisect import bisect_left, bisect_right  # noqa: F401
from collections import Counter, defaultdict, deque  # noqa: F401
from heapq import heappop, heappush  # noqa: F401
from math import gcd, inf, isqrt  # noqa: F401


def main():
  data = sys.stdin.buffer.read().split()
  pos = 0

  def nxt():
    nonlocal pos
    pos += 1
    return int(data[pos - 1])

  n = nxt()
  a = [nxt() for _ in range(n)]
  out = [sum(a)]
  sys.stdout.write("\n".join(map(str, out)) + "\n")


main()
