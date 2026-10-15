from testlib import CONTENT, Timer, check_eq, iters, load, rnd

lib = load("offline/cdq")


def brute(p):
  return [
    sum(
      j != i and all(p[j][d] <= p[i][d] for d in range(3))
      for j in range(len(p))
    )
    for i in range(len(p))
  ]


def random_points(n, lo, hi):
  return [(rnd(lo, hi), rnd(lo, hi), rnd(lo, hi)) for _ in range(n)]


def cdq_time(X, Y, ask):
  """The Example's CDQ on time: op t inserts (X[t], Y[t]) or, if
  ask[t], counts inserted points <= (X[t], Y[t])."""
  q = len(X)
  rank = {v: i for i, v in enumerate(sorted(set(X)))}
  cx = [rank[x] for x in X]
  idx, res, fw = list(range(q)), [0] * q, lib.Fenwick(q)
  add, query = fw.add, fw.query

  def rec(l, r):
    if r - l < 2:
      return
    m = (l + r) // 2
    rec(l, m)
    rec(m, r)
    k, merged = l, []
    for i in range(m, r):
      t = idx[i]
      while k < m and Y[idx[k]] <= Y[t]:
        if not ask[idx[k]]:
          add(cx[idx[k]], 1)
        merged.append(idx[k])
        k += 1
      if ask[t]:
        res[t] += query(cx[t] + 1)
      merged.append(t)
    for j in range(l, k):
      if not ask[idx[j]]:
        add(cx[idx[j]], -1)
    idx[l:r] = merged + idx[k:m]

  rec(0, q)
  return res


def weighted_from_tweak():
  """dominance with the Tweak's weights, applied to the snippet's own
  source, so the names the Tweak mentions must be the code's."""
  src = (CONTENT / "offline/cdq.py").read_text()
  for old, new in (
    ("def dominance(p):", "def dominance(p, wt):"),
    ("add(c[y], 1)", "add(c[y], wt[y])"),
    ("add(c[idx[j]], -1)", "add(c[idx[j]], -wt[idx[j]])"),
    (  # cnt[last] + wt[last] - wt[i], telescoped along the copies
      "cnt[order[i]] = cnt[order[i + 1]]",
      "cnt[order[i]] = cnt[order[i + 1]] + wt[order[i + 1]] - wt[order[i]]",
    ),
  ):
    check_eq(src.count(old), 1, old)
    src = src.replace(old, new)
  namespace = {"Fenwick": lib.Fenwick}
  exec(compile(src, "cdq.py with weights", "exec"), namespace)
  return namespace["dominance"]


def stress_weighted(rounds, maxn, V):
  dominance = weighted_from_tweak()
  for _ in range(iters(rounds)):
    n = rnd(0, maxn)
    p, wt = random_points(n, -V, V), [rnd(-9, 9) for _ in range(n)]
    want = [
      sum(
        wt[j]
        for j in range(n)
        if j != i and all(p[j][d] <= p[i][d] for d in range(3))
      )
      for i in range(n)
    ]
    check_eq(dominance(p, wt), want, (p, wt))


def stress_time(rounds, maxq, V):
  for _ in range(iters(rounds)):
    q = rnd(0, maxq)
    X = [rnd(0, V) for _ in range(q)]
    Y = [rnd(0, V) for _ in range(q)]
    ask = [rnd(0, 1) for _ in range(q)]
    res = cdq_time(X, Y, ask)
    for t in range(q):
      if ask[t]:
        want = sum(
          not ask[s] and X[s] <= X[t] and Y[s] <= Y[t] for s in range(t)
        )
        check_eq(res[t], want)


# Fixed cases: the packet's example, copies, n = 0 and 1.
check_eq(
  lib.dominance([(1, 3, 2), (2, 1, 1), (3, 2, 3), (4, 4, 2)]),
  [0, 0, 1, 2],
)
check_eq(lib.dominance([]), [])
check_eq(lib.dominance([(7, -3, 2)]), [0])
check_eq(lib.dominance([(4, 4, 4)] * 5), [4] * 5)
check_eq(lib.dominance([[1, 2, 3], [1, 2, 3], [0, 9, 3]]), [1, 1, 0])
M = 10**18
check_eq(
  lib.dominance([(M, M, M), (-M, -M, -M), (M, -M, M), (-M, M, -M)]),
  [3, 0, 1, 1],
)
# Small ranges give many ties in every coordinate.
for rounds, maxn, V in ((2000, 12, 2), (2000, 25, 4), (300, 60, 50)):
  for _ in range(iters(rounds)):
    p = random_points(rnd(0, maxn), -V, V)
    check_eq(lib.dominance(p), brute(p))
stress_weighted(1500, 20, 2)
stress_time(1500, 30, 5)
stress_time(200, 120, 1000)

for n, V, label in ((100_000, 10**9, "n = 1e5"), (200_000, 1000, "n = 2e5")):
  p = random_points(n, 0, V)
  with Timer(label):
    cnt = lib.dominance(p)
  for _ in range(50):
    i = rnd(0, n - 1)
    a, b, c = p[i]
    want = sum(x <= a and y <= b and z <= c for x, y, z in p) - 1
    check_eq(cnt[i], want)
print("OK offline/cdq")
