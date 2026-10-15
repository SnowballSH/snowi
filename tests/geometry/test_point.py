import math
from functools import cmp_to_key

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("geometry/point")


def angle_key(v):
  """Angle in [0, 2 pi) of the reduced direction: ties are exact."""
  g = math.gcd(v[0], v[1])
  a = math.atan2(v[1] // g, v[0] // g)
  return a + 2 * math.pi if a < 0 else a


def brute_turn(o, a, b):
  """Orientation of o -> a -> b from angles, no cross products."""
  if a == o or b == o:
    return 0
  s = math.sin(
    math.atan2(b[1] - o[1], b[0] - o[0]) - math.atan2(a[1] - o[1], a[0] - o[0])
  )
  return 0 if abs(s) < 1e-9 else (1 if s > 0 else -1)


def quadrant(v):
  x, y = v[0], v[1]
  if x > 0 and y >= 0:
    return 0
  if x <= 0 and y > 0:
    return 1
  if x < 0 and y <= 0:
    return 2
  return 3


def exact_less(a, b):
  qa, qb = quadrant(a), quadrant(b)
  if qa != qb:
    return qa < qb
  return a[0] * b[1] - a[1] * b[0] > 0


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def rnd_nonzero(c):
  while True:
    v = rnd_pt(c)
    if v != (0, 0):
      return v


for _ in range(iters(20000)):
  o, a, b = rnd_pt(20), rnd_pt(20), rnd_pt(20)
  check_eq(lib.sub(a, b), (a[0] - b[0], a[1] - b[1]))
  check_eq(lib.dot(a, b), a[0] * b[0] + a[1] * b[1])
  check_eq(lib.dist2(a, b), (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)
  t = lib.cross3(o, a, b)
  check_eq((t > 0) - (t < 0), brute_turn(o, a, b))
  c = lib.cross(a, b)
  check_eq((c > 0) - (c < 0), brute_turn((0, 0), a, b))
  check_eq(t, lib.cross(lib.sub(a, o), lib.sub(b, o)))
  if a == (0, 0) or b == (0, 0):
    continue
  ka, kb = angle_key(a), angle_key(b)
  check_eq(lib.angle_cmp(a, b), (ka > kb) - (ka < kb))
  check_eq(lib.half(a), ka >= math.pi)

for _ in range(iters(300)):
  n = rnd(0, 40)
  v = [rnd_nonzero(rnd(1, 6)) + (i,) for i in range(n)]
  w = v[:]
  lib.angle_sort(w)
  check_eq(sorted(w, key=lambda p: p[2]), v)  # a permutation
  check_eq([angle_key(p) for p in w], sorted(angle_key(p) for p in v))

# Extremes: Python ints are exact; atan2 alone cannot separate these.
C = 10**9
u, w = (C, C - 1), (C - 1, C - 2)
check_eq(lib.cross(u, w), -1)
check_eq(lib.angle_cmp(w, u), -1)
v = [u, w] * 3 + [(C, -1), (1, 0), (-C, 1), (-1, 0), (C, 1)]
lib.angle_sort(v)
check_eq(v, sorted(v, key=cmp_to_key(lib.angle_cmp)))
check_eq(v[:2], [(1, 0), (C, 1)])
check_eq(v[-1], (C, -1))
for _ in range(iters(2000)):
  a, b = rnd_nonzero(C), rnd_nonzero(C)
  check_eq(lib.angle_cmp(a, b) < 0, exact_less(a, b))
  check_eq(lib.angle_cmp(a, a), 0)

# Near-parallel directions at 1e9 share one atan2 float (cross = j - k)
# and come with copies and with directions that atan2 does separate.
pool = [(C - k, C - k - 1) for k in range(4)]
pool += [(-x, -y) for x, y in pool] + [(C, -1), (C, -2), (1, 0), (-1, 0)]
for _ in range(iters(1000)):
  v = [pool[rnd(0, len(pool) - 1)] + (i,) for i in range(rnd(0, 30))]
  w = v[:]
  lib.angle_sort(w)
  check_eq(sorted(w, key=lambda p: p[2]), v)
  check(all(not exact_less(w[i + 1], w[i]) for i in range(len(w) - 1)), w)

# Copies of two directions with equal atan2 floats must not make the
# repair quadratic (it took 5 s at n = 16000, 33 s at n = 40000).
m, a, b = 20_000, (C, C - 1), (C - 1, C - 2)
v = [a] * m + [b] * m + [(C, 1), (1, C)] * (m // 2)
with Timer("angle_sort 6e4: two float-tied directions, copies") as t:
  lib.angle_sort(v)
check(t.seconds < 5, "angle_sort repair is quadratic")
check_eq(v, [(C, 1)] * (m // 2) + [b] * m + [a] * m + [(1, C)] * (m // 2))

n = 200_000
v = [rnd_nonzero(C) for _ in range(n)]
v += [(2 * x, 2 * y) for x, y in v[:1000]]  # repeated directions
with Timer("angle_sort 2e5 vectors"):
  lib.angle_sort(v)
for i in range(len(v) - 1):
  check(not exact_less(v[i + 1], v[i]), "not sorted")
print("OK geometry/point")
