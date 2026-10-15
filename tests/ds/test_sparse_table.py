import resource

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("ds/sparse-table")


def peak_mb():
  return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


# The Notes' size, run first so the peak memory reading is the table's:
# n = 1e6 takes about 150 MB (a slice copy per level made it 250 MB).
n = 10**6
a = [rnd(0, 10**9) for _ in range(n)]
queries = [(rnd(0, n // 2), rnd(n // 2 + 1, n)) for _ in range(n)]
before = peak_mb()
with Timer("build n = 1e6 + 1e6 queries"):
  st = lib.SparseTable(a)
  got = [st.query(l, r) for l, r in queries]
grown = peak_mb() - before
check(grown < 200, f"table took {grown} MB")
mid, part = n // 2, a[:]  # every query range contains mid
for i in range(mid - 1, -1, -1):
  part[i] = min(part[i], part[i + 1])  # min of a[i..mid]
for i in range(mid + 1, n):
  part[i] = min(part[i], part[i - 1])  # min of a[mid..i]
want = [min(part[l], part[r - 1]) for l, r in queries]
check(got == want, "n = 1e6: some query is wrong")
del st, a, queries, got, want, part

for _ in range(iters(2000)):
  n = rnd(0, 40)
  big = rnd(0, 1)
  a = [rnd(-(10**30), 10**30) if big else rnd(-5, 5) for _ in range(n)]
  st = lib.SparseTable(a)
  packed = lib.SparseTable([(x << 20) + i for i, x in enumerate(a)])
  for l in range(n):
    m = a[l]
    for r in range(l + 1, n + 1):
      m = min(m, a[r - 1])
      check_eq(st.query(l, r), m)
      check_eq(packed.query(l, r) & (2**20 - 1), a.index(m, l, r))

n, B = 200_000, 256
a = [rnd(0, 10**9) for _ in range(n)]
queries = []
for _ in range(n):
  l = rnd(0, n - 1)
  queries.append((l, rnd(l + 1, n)))
with Timer("build n = 2e5 + 2e5 queries"):
  st = lib.SparseTable(a)
  got = [st.query(l, r) for l, r in queries]
block = [min(a[i : i + B]) for i in range(0, n, B)]
for (l, r), g in zip(queries, got):
  m, i = 10**18, l
  while i < r and i % B:
    m = min(m, a[i])
    i += 1
  while i + B <= r:
    m = min(m, block[i // B])
    i += B
  while i < r:
    m = min(m, a[i])
    i += 1
  check_eq(g, m)
print("OK ds/sparse-table")
