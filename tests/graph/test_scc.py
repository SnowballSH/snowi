from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("graph/scc")


def reach(adj):
  n = len(adj)
  r = [[False] * n for _ in range(n)]
  for s in range(n):
    r[s][s] = True
    stack = [s]
    while stack:
      for v in adj[stack.pop()]:
        if not r[s][v]:
          r[s][v] = True
          stack.append(v)
  return r


def check_order(adj, comp, k):
  check_eq(len(comp), len(adj))
  check_eq(sorted(set(comp)), list(range(k)))
  for u, out in enumerate(adj):
    for v in out:
      check(comp[u] >= comp[v])


for _ in range(iters(3000)):
  n = rnd(0, 9)
  adj = [[] for _ in range(n)]
  for _ in range(rnd(0, 2 * n + 2) if n else 0):
    adj[rnd(0, n - 1)].append(rnd(0, n - 1))
  comp, k = lib.scc(adj)
  check_order(adj, comp, k)
  r = reach(adj)
  for u in range(n):
    for v in range(n):
      check_eq(comp[u] == comp[v], r[u][v] and r[v][u])

# Blocks of cycles chained into a path, labels shuffled. The
# condensation is a path, so the numbering is forced.
n, b = 500_000, 1000
k = n // b
label = list(range(n))
rng.shuffle(label)
adj = [[] for _ in range(n)]
for v in range(n):
  adj[label[v]].append(label[v - b + 1 if v % b == b - 1 else v + 1])
  if v + b < n:
    adj[label[v]].append(label[v + rnd(1, b)])
with Timer("n = 5e5, m = 1e6"):
  comp, cnt = lib.scc(adj)
check_eq(cnt, k)
check(all(comp[label[v]] == k - 1 - v // b for v in range(n)))
cyc = [[(v + 1) % n] for v in range(n)]
with Timer("one cycle of 5e5 vertices"):
  check_eq(lib.scc(cyc)[1], 1)
print("OK graph/scc")
