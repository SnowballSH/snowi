from itertools import product

from testlib import Timer, check, check_eq, iters, load, rnd, rnd_list

lib = load("dp/slope-trick")
INF = float("inf")


class Grid:
  """f tabulated on [-G, G]; breakpoints stay well inside, so f is
  linear beyond the grid and value() extrapolates it exactly."""

  def __init__(self, G):
    self.G, self.v = G, [0] * (2 * G + 1)

  def at(self, x):
    return self.v[x + self.G]

  def value(self, x):
    G, v = self.G, self.v
    if x < -G:
      return v[0] + (v[1] - v[0]) * (x + G)
    if x > G:
      return v[-1] + (v[-1] - v[-2]) * (x - G)
    return v[x + G]

  def add(self, f):
    self.v = [y + f(x - self.G) for x, y in enumerate(self.v)]

  def prefix_min(self):
    for i in range(1, len(self.v)):
      self.v[i] = min(self.v[i], self.v[i - 1])

  def suffix_min(self):
    for i in range(len(self.v) - 2, -1, -1):
      self.v[i] = min(self.v[i], self.v[i + 1])

  def slide(self, a, b):
    G = self.G
    self.v = [
      min(self.value(y) for y in range(x - b, x - a + 1))
      for x in range(-G, G + 1)
    ]


def compare(st, g, W):
  check(st.top_l() <= st.top_r())
  xs = range(-W, W + 1) if W <= 30 else [rnd(-W, W) for _ in range(20)]
  for x in xs:
    check_eq(st.eval(x), g.at(x), x)
  mn = min(g.v)
  check_eq(st.min_f, mn)
  lo = g.v.index(mn)
  hi = len(g.v) - 1 - g.v[::-1].index(mn)
  # An argmin reaching the grid edge with slope 0 is unbounded.
  check_eq(st.top_l(), -INF if lo == 0 and g.v[1] == mn else lo - g.G)
  check_eq(st.top_r(), INF if hi == 2 * g.G and g.v[-2] == mn else hi - g.G)


def stress(rounds, ops, c, d):
  for _ in range(iters(rounds)):
    st, g = lib.SlopeTrick(), Grid(c + ops * 2 * d + 40)
    for op in range(ops):
      a, kind = rnd(-c, c), rnd(0, 7)
      if kind == 0:
        k = rnd(-5, 5)
        st.add_const(k)
        g.add(lambda x, k=k: k)
      elif kind == 1:
        st.add_x_minus_a(a)
        g.add(lambda x, a=a: max(0, x - a))
      elif kind == 2:
        st.add_a_minus_x(a)
        g.add(lambda x, a=a: max(0, a - x))
      elif kind <= 4:
        st.add_abs(a)
        g.add(lambda x, a=a: abs(x - a))
      elif kind == 5:
        st.prefix_min()
        g.prefix_min()
      elif kind == 6:
        st.suffix_min()
        g.suffix_min()
      else:
        lo = rnd(-d, d)
        hi = rnd(lo, d)
        st.slide(lo, hi)
        g.slide(lo, hi)
      compare(st, g, c + (op + 1) * d + 3)


stress(800, 12, 3, 1)
stress(400, 25, 10, 3)

for _ in range(iters(300)):  # extremes: hinges only, explicit sum
  st, hs, e = lib.SlopeTrick(), [], 10**17
  for _ in range(20):
    t, a = rnd(0, 2), rnd(-e, e)
    hs.append((t, a))
    (st.add_x_minus_a, st.add_a_minus_x, st.add_abs)[t](a)

  def f(x, hs=hs):
    return sum(
      (max(0, x - a) if t != 1 else 0) + (max(0, a - x) if t != 0 else 0)
      for t, a in hs
    )

  for _, a in hs:
    check_eq(st.eval(a), f(a))
  check_eq(st.min_f, min(f(a) for _, a in hs))
  for x in (st.top_l(), st.top_r()):
    if abs(x) != INF:
      check_eq(f(x), st.min_f)


def bounded_steps(rounds):
  """Regression for the documented variants: fewest +-1 steps so that
  |b[i] - b[i-1]| <= d (slide(-d, d) in place of prefix_min), or
  0 <= b[i] - b[i-1] <= d (slide(0, d)), against every b in
  [min a, max a]^n (clamping b there keeps it feasible)."""
  for _ in range(iters(rounds)):
    a, d = rnd_list(rnd(1, 4), 0, 6), rnd(0, 2)
    for lo in (-d, 0):
      st = lib.SlopeTrick()
      for x in a:
        st.slide(lo, d)
        st.add_abs(x)
      want = min(
        sum(abs(x - y) for x, y in zip(a, b))
        for b in product(range(7), repeat=len(a))
        if all(lo <= b[i] - b[i - 1] <= d for i in range(1, len(a)))
      )
      check_eq(st.min_f, want, (a, d, lo))


bounded_steps(400)

st = lib.SlopeTrick()  # empty function
check_eq((st.min_f, st.top_l(), st.top_r()), (0, -INF, INF))
st.prefix_min()
st.suffix_min()
st.slide(-3, 5)
check_eq(st.eval(12345), 0)


def brute_isotonic(a, V):
  """Fewest +-1 steps making a non-decreasing, DP over values."""
  dp = [0] * (V + 1)
  for x in a:
    for v in range(1, V + 1):
      dp[v] = min(dp[v], dp[v - 1])
    dp = [y + abs(x - v) for v, y in enumerate(dp)]
  return min(dp)


n = 200_000
a = rnd_list(n, 0, 30)
with Timer("2e5 isotonic steps"):
  st = lib.SlopeTrick()
  for x in a:
    st.prefix_min()
    st.add_abs(x)
check_eq(st.min_f, brute_isotonic(a, 30))
with Timer("2e5 mixed operations"):
  mix = lib.SlopeTrick()
  for _ in range(n):
    v, kind = rnd(-(10**9), 10**9), rnd(0, 9)
    if kind == 0:
      mix.prefix_min()
    elif kind == 1:
      mix.suffix_min()
    elif kind == 2:
      mix.slide(-abs(v) // 1000, abs(v) // 1000)
    else:
      mix.add_abs(v)
if mix.top_l() != -INF:
  check_eq(mix.eval(mix.top_l()), mix.min_f)
for _ in range(iters(300)):  # Tweak: start pinned at x0
  a, x0, d = rnd_list(rnd(1, 3), -4, 6), rnd(-3, 3), rnd(0, 2)
  st = lib.SlopeTrick()
  for _ in range(50):
    st.add_abs(x0)
  for x in a:
    st.slide(-d, d)
    st.add_abs(x)
  want = min(
    sum(abs(x - y) for x, y in zip(a, b))
    for b in product(range(-12, 13), repeat=len(a))
    if abs(b[0] - x0) <= d
    and all(abs(b[i] - b[i - 1]) <= d for i in range(1, len(a)))
  )
  check_eq(st.min_f, want, (a, x0, d))
print("OK dp/slope-trick")
