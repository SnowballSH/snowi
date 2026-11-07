from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("ds/dsu")

for _ in range(iters(3000)):
  n = rnd(1, 15)
  d = lib.DSU(n)
  label = list(range(n))
  for _ in range(30):
    a, b = rnd(0, n - 1), rnd(0, n - 1)
    if rnd(0, 1):
      check_eq(d.join(a, b), label[a] != label[b])
      old, new = label[b], label[a]
      label = [new if x == old else x for x in label]
    check_eq(d.same(a, b), label[a] == label[b])
    check_eq(label[d.find(a)], label[a])
    check_eq(d.size(a), label.count(label[a]))

n = 200_000
pairs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("2e5 joins + 2e5 finds"):
  d = lib.DSU(n)
  merges = sum(d.join(a, b) for a, b in pairs)
  roots = [d.find(i) for i in range(n)]
check_eq(len(set(roots)), n - merges)
check_eq(sum(d.size(i) for i in set(roots)), n)
# Union by size keeps depth <= log2 n; one find on the deepest
# vertex of a doubling tree at least halves its depth.
n = 1 << 12
d = lib.DSU(n)
for i in range(1, n):
  if i % 2:
    d.join(i, 0)
  else:
    d.join(0, i)
for x in range(n):
  depth = 0
  while d.parent[x] != x:
    x, depth = d.parent[x], depth + 1
  check(depth <= 12, f"depth {depth}")
d = lib.DSU(n)
length = 1
while length < n:
  for i in range(0, n - length, 2 * length):
    d.join(i, i + length)
  length *= 2
check_eq(d.find(n - 1), 0)
x, depth = n - 1, 0
while d.parent[x] != x:
  x, depth = d.parent[x], depth + 1
check(depth <= 6, f"depth {depth} after one find")
print("OK ds/dsu")
