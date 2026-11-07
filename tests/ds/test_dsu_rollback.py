from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("ds/dsu-rollback")

for _ in range(iters(2000)):
  n = rnd(1, 12)
  d = lib.RollbackDSU(n)
  label = list(range(n))
  saved = []  # (checkpoint, labels)
  for _ in range(40):
    kind = rnd(0, 4)
    a, b = rnd(0, n - 1), rnd(0, n - 1)
    if kind <= 1:
      check_eq(d.join(a, b), label[a] != label[b])
      old, new = label[b], label[a]
      label = [new if x == old else x for x in label]
    elif kind == 2:
      saved.append((d.time(), label))
    elif saved:
      k = rnd(0, len(saved) - 1)
      d.rollback(saved[k][0])
      label = saved[k][1]
      del saved[k + 1 :]
    check_eq(d.find(a) == d.find(b), label[a] == label[b])
    check_eq(label[d.find(a)], label[a])
    check_eq(d.size(a), label.count(label[a]))
    check_eq(n - d.time(), len(set(label)))  # the Tweak's count

n = 200_000
pairs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("2e5 joins + 2e5 finds + rollbacks"):
  d = lib.RollbackDSU(n)
  merges = 0
  for i, (a, b) in enumerate(pairs):
    t = d.time()
    merged = d.join(a, b)
    if i % 3 == 0:
      d.rollback(t)
    else:
      merges += merged
  roots = set(d.find(i) for i in range(n))
check_eq(len(roots), n - merges)
check_eq(sum(d.size(r) for r in roots), n)
# Union by size: depth <= log2 n. A flipped comparison would turn
# one of these two join orders into a path.
d = lib.RollbackDSU(1 << 12)
for i in range(1, 1 << 12):
  if i % 2:
    d.join(i, 0)
  else:
    d.join(0, i)
for x in range(1 << 12):
  depth = 0
  while d.parent[x] != x:
    x, depth = d.parent[x], depth + 1
  check(depth <= 12, f"depth {depth}")
print("OK ds/dsu-rollback")
