from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("linalg/xor-basis")
XorBasis, PrefixXorBasis = lib.XorBasis, lib.PrefixXorBasis


def span_of(a):
  """All subset XORs of a, sorted."""
  s = {0}
  for x in a:
    s |= {v ^ x for v in s}
  return sorted(s)


def random_value(width, prev):
  if prev and rnd(0, 3) == 0:  # dependent on earlier ones
    x = 0
    for v in prev:
      if rnd(0, 1):
        x ^= v
    return x
  return rnd(0, (1 << width) - 1)


for _ in range(iters(1500)):
  width = rnd(1, 6) if rnd(0, 4) else 60
  bits = 60 if width == 60 or rnd(0, 1) else width  # Tweak: bits
  n = rnd(0, 10)
  xb = XorBasis(bits)
  a = []
  for _ in range(n):
    x = random_value(width, a)
    before = len(span_of(a))
    a.append(x)
    s = span_of(a)
    check_eq(xb.insert(x), len(s) > before)
    check_eq(1 << xb.rank, len(s))
    check_eq([xb.kth(k) for k in range(len(s))], s)
    for q in range(5):
      y = 0 if q == 0 else random_value(width, a)
      check_eq(xb.max_xor(y), max(y ^ v for v in s))
      check_eq(xb.min_xor(y), min(y ^ v for v in s))
      check_eq(xb.contains(y), y in s)

ones = (1 << 60) - 1
full = XorBasis()  # all-ones minus one bit: independent (60 even)
check(all(full.insert(ones ^ 1 << b) for b in range(60)))
check(full.rank == 60 and not full.insert(1 << 59))
check_eq((full.max_xor(), full.min_xor(ones)), (ones, 0))
check_eq([full.kth(k) for k in (0, 1, 12345, ones)], [0, 1, 12345, ones])

for _ in range(iters(1500)):
  width = rnd(1, 5) if rnd(0, 4) else 60
  n = rnd(1, 8 if width == 60 else 12)
  a = []
  pb = PrefixXorBasis()
  for r in range(n):
    a.append(random_value(width, a))
    pb.insert(a[r], r)
    s = {0}  # span of a[l..r], growing as l decreases
    for l in range(r, -1, -1):
      s |= {v ^ a[l] for v in s}
      check_eq(pb.max_xor(l), max(s))
      sub = XorBasis()  # the slots with pos >= l are a basis
      for b in range(60):
        if pb.p[b] and pb.pos[b] >= l:
          check(sub.insert(pb.p[b]))
      check_eq([sub.kth(k) for k in range(len(s))], sorted(s))


def speed(n):
  a = [rnd(0, (1 << 60) - 1) for _ in range(n)]
  ls = [rnd(0, r) for r in range(n)]
  with Timer("2e5 inserts + 2e5 queries (XorBasis)"):
    xb = XorBasis()
    total = 0
    for r in range(n):
      xb.insert(a[r])
      total += xb.max_xor(a[r]) & 1
  with Timer("2e5 inserts + 2e5 range queries (PrefixXorBasis)"):
    pb = PrefixXorBasis()
    for r in range(n):
      pb.insert(a[r], r)
      total += pb.max_xor(ls[r]) & 1
  check_eq(pb.max_xor(0), ones)
  check_eq(xb.max_xor(), ones)


speed(200_000)
print("OK linalg/xor-basis")
