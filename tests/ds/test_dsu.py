from testlib import Timer, check_eq, iters, load, rnd

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
print("OK ds/dsu")
