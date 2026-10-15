from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("offline/dynamic-connectivity")


def labels(n, find):
  """Component label of every vertex: the smallest vertex in it."""
  root = [find(v) for v in range(n)]
  low = {}
  for v in range(n):
    low.setdefault(root[v], v)
  return [low[r] for r in root]


def brute_labels(n, edges):
  p = list(range(n))

  def find(x):
    while p[x] != x:
      p[x] = p[p[x]]
      x = p[x]
    return x

  for u, v in edges:
    p[find(u)] = find(v)
  return labels(n, find)


def stress(n, ops, p_add, p_query, samples=None):
  """Random ops: add (loops and parallel edges allowed), remove a
  present edge, or query. Checks every query's partition and
  component count, and that answers come once each, in order."""
  dc = lib.DynConn(n)
  present, want = [], []
  for _ in range(ops):
    x = rng.random()
    if x < p_query:
      check_eq(dc.query(), len(want))
      sampled = samples is None or rng.random() * ops * p_query < samples
      want.append(brute_labels(n, present) if sampled else None)
    elif x < p_query + p_add or not present:
      u, v = rnd(0, n - 1), rnd(0, n - 1)
      dc.add(u, v)
      present.append((u, v))
    else:
      k = rnd(0, len(present) - 1)
      present[k], present[-1] = present[-1], present[k]
      u, v = present.pop()
      if rnd(0, 1):
        u, v = v, u
      dc.remove(u, v)
  seen = []

  def answer(i):
    check_eq(i, len(seen))
    seen.append(i)
    if want[i] is not None:
      check_eq(labels(n, dc.d.find), want[i])
      check_eq(n - dc.d.time(), len(set(want[i])))

  dc.solve(answer)
  check_eq(len(seen), len(want))
  check_eq(dc.d.time(), 0)  # everything rolled back


# the packet's example, 0-indexed
dc = lib.DynConn(4)
dc.add(0, 1)
dc.add(1, 2)
dc.query()
dc.remove(0, 1)
dc.add(2, 3)
dc.query()
out = []
dc.solve(lambda i: out.append((4 - dc.d.time(), dc.d.find(0) == dc.d.find(2))))
check_eq(out, [(2, True), (2, False)])
# no queries; queries on an empty graph; one vertex
a = lib.DynConn(3)
a.add(0, 1)
a.remove(1, 0)
a.add(1, 2)
a.solve(lambda i: check(False))
b = lib.DynConn(1)
b.query()
b.add(0, 0)
b.query()
b.remove(0, 0)
b.query()
calls = []
b.solve(lambda i: calls.append(b.d.time()))
check_eq(calls, [0, 0, 0])

for _ in range(iters(2000)):
  stress(rnd(1, 7), rnd(0, 40), rng.uniform(0.2, 0.6), rng.uniform(0.1, 0.5))
for _ in range(iters(100)):
  stress(30, 300, 0.4, 0.3)
stress(1000, 20_000, 0.4, 0.3, 30)


def components_after_each(n, ops):  # the Usage, timed
  dc = lib.DynConn(n)
  for op, u, v in ops:
    if op == 0:
      dc.add(u, v)
    elif op == 1:
      dc.remove(u, v)
    else:
      dc.query()
  res = []
  d = dc.d
  dc.solve(lambda i: res.append(n - d.time()))
  return res


n = 100_000
ops, present = [], []
for _ in range(200_000):
  x = rng.random()
  if x < 0.25:
    ops.append((2, 0, 0))
  elif x < 0.65 or not present:
    u, v = rnd(0, n - 1), rnd(0, n - 1)
    ops.append((0, u, v))
    present.append((u, v))
  else:
    k = rnd(0, len(present) - 1)
    present[k], present[-1] = present[-1], present[k]
    ops.append((1,) + present.pop())
with Timer("2e5 operations, n = 1e5"):
  res = components_after_each(n, ops)
# the last query against a plain union-find
last = max(i for i, op in enumerate(ops) if op[0] == 2)
alive = {}
for op, u, v in ops[:last]:
  if op < 2:
    alive[(u, v)] = alive.get((u, v), 0) + (1 - 2 * op)
edges = [e for e, c in alive.items() for _ in range(c)]
check_eq(res[-1], len(set(brute_labels(n, edges))))
print("OK offline/dynamic-connectivity")
