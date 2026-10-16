from bisect import bisect_left, bisect_right

from testlib import Timer, check, check_eq, iters, load, rnd

# parallel_bs has no Uses; Fenwick is the Usage's and Example's state.
lib = load("ds/fenwick", "offline/parallel-binary-search")


def run(T, Q, reset, apply, pred):
  """Wraps the callbacks: counts calls and checks that apply(t)
  comes in order 0, 1, ... after each reset."""
  calls = {"reset": 0, "apply": 0, "pred": 0}
  nxt = [-1]

  def reset2():
    reset()
    calls["reset"] += 1
    nxt[0] = 0

  def apply2(t):
    check_eq(t, nxt[0])
    apply(t)
    nxt[0] += 1
    calls["apply"] += 1

  def pred2(i):
    check(nxt[0] >= 0 and 0 <= i < Q)
    calls["pred"] += 1
    return pred(i)

  res = lib.parallel_bs(T, Q, reset2, apply2, pred2)
  check_eq(len(res), Q)
  rounds = (T + 1).bit_length()  # ceil(log2(T + 2))
  check(calls["reset"] <= rounds, calls)
  check(calls["apply"] <= calls["reset"] * T)
  check(calls["pred"] <= calls["reset"] * Q)
  return res, calls


def meteors(m, own, l, r, x, need, solve):
  """The Example's METEORS code: a ring of m sectors, sector j owned
  by state own[j]; shower t adds x[t] to sectors l[t]..r[t]."""
  T, Q = len(l), len(need)
  sec = [[] for _ in range(Q)]
  for j in range(m):
    sec[own[j]].append(j)
  fw = lib.Fenwick(m)

  def reset():
    fw.t = [0] * m

  def apply(t):
    fw.add(l[t], x[t])
    fw.add(r[t] + 1, -x[t])
    if l[t] > r[t]:
      fw.add(0, x[t])

  def pred(i):
    s = 0
    for j in sec[i]:
      s += fw.query(j + 1)
      if s >= need[i]:
        break
    return s >= need[i]

  return solve(T, Q, reset, apply, pred)


def random_meteors(m, T, Q, maxx, maxneed):
  own = [rnd(0, Q - 1) for _ in range(m)]
  l = [rnd(0, m - 1) for _ in range(T)]
  r = [rnd(0, m - 1) for _ in range(T)]
  x = [rnd(0, maxx) for _ in range(T)]
  return own, l, r, x, [rnd(0, maxneed) for _ in range(Q)]


def meteors_one(m, own, l, r, x, need, i):
  """Brute force for one state: walk the showers, counting its
  sectors inside each one by bisection."""
  sec = [j for j in range(m) if own[j] == i]
  s, t = 0, 0
  while s < need[i] and t < len(l):
    a, b = l[t], r[t]
    if a <= b:
      s += x[t] * (bisect_right(sec, b) - bisect_left(sec, a))
    else:
      s += x[t] * (len(sec) - bisect_left(sec, a) + bisect_right(sec, b))
    t += 1
  return t if s >= need[i] else len(l) + 1


def stress_meteors(maxm, maxT, maxQ, maxx):
  m, T, Q = rnd(1, maxm), rnd(0, maxT), rnd(1, maxQ)  # every sector owned
  own, l, r, x, need = random_meteors(m, T, Q, maxx, 0)
  sec = [[j for j in range(m) if own[j] == i] for i in range(Q)]
  need = [rnd(0, maxx * T * len(sec[i]) // 4 + 1) for i in range(Q)]
  a, want = [0] * m, [T + 1] * Q
  for t in range(T + 1):
    for i in range(Q):
      if want[i] > T and sum(a[j] for j in sec[i]) >= need[i]:
        want[i] = t
    if t < T:
      j = l[t]
      while True:
        a[j] += x[t]
        if j == r[t]:
          break
        j = (j + 1) % m
  res, _ = meteors(m, own, l, r, x, need, run)
  check_eq(res, want)


def stress_dsu(maxn, maxT, maxQ):
  """DSU state: the first prefix of edges joining u and v."""
  n, T, Q = rnd(1, maxn), rnd(0, maxT), rnd(0, maxQ)
  e = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(T)]
  qs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(Q)]
  p = list(range(n))

  def find(v):
    while p[v] != v:
      v = p[v]
    return v

  want = []
  for u, v in qs:
    p[:] = range(n)
    t = 0
    while t <= T and find(u) != find(v):
      if t < T:
        p[find(e[t][0])] = find(e[t][1])
      t += 1
    want.append(t)  # T + 1 if never

  def reset():
    p[:] = range(n)

  def apply(t):
    p[find(e[t][0])] = find(e[t][1])

  res, _ = run(T, Q, reset, apply, lambda i: find(qs[i][0]) == find(qs[i][1]))
  check_eq(res, want)


# Q = 0; T = 0
res, calls = run(5, 0, lambda: None, lambda t: None, lambda i: True)
check_eq((res, calls["reset"]), ([], 0))
res, _ = run(0, 3, lambda: None, lambda t: None, lambda i: i != 1)
check_eq(res, [0, 1, 0])
for _ in range(iters(2000)):
  stress_meteors(6, 10, 5, 4)
for _ in range(iters(200)):
  stress_meteors(30, 100, 30, 10**9)
for _ in range(iters(1500)):
  stress_dsu(7, 12, 8)
for _ in range(iters(100)):
  stress_dsu(30, 200, 50)


def usage(n, T, Q, p, w, r, need):
  """The Usage verbatim, inside a function as in main()."""
  Fenwick, parallel_bs = lib.Fenwick, lib.parallel_bs
  fw = Fenwick(n)

  def reset():
    nonlocal fw
    fw = Fenwick(n)

  return parallel_bs(
    T,
    Q,
    reset,
    lambda t: fw.add(p[t], w[t]),
    lambda i: fw.query(r[i]) >= need[i],
  )


for _ in range(iters(300)):  # small: against a direct simulation
  n, T, Q = rnd(1, 8), rnd(0, 10), rnd(0, 8)
  p, w = [rnd(0, n - 1) for _ in range(T)], [rnd(0, 3) for _ in range(T)]
  r, need = [rnd(0, n) for _ in range(Q)], [rnd(0, 8) for _ in range(Q)]
  want = []
  for i in range(Q):
    s, t = 0, 0
    while t < T and s < need[i]:
      s += w[t] if p[t] < r[i] else 0
      t += 1
    want.append(t if s >= need[i] else T + 1)
  check_eq(usage(n, T, Q, p, w, r, need), want)

# Large: the Usage, T = Q = 2e5 point adds and prefix checks.
n = T = Q = 200_000
p = [rnd(0, n - 1) for _ in range(T)]
w = [rnd(0, 10**9) for _ in range(T)]
r = [rnd(0, n) for _ in range(Q)]
need = [rnd(0, 10**14) for _ in range(Q)]
with Timer("T = Q = 2e5, Usage"):
  res = usage(n, T, Q, p, w, r, need)
for _ in range(50):
  i, s, t = rnd(0, Q - 1), 0, 0
  while t < T and s < need[i]:
    s += w[t] if p[t] < r[i] else 0
    t += 1
  check_eq(res[i], t if s >= need[i] else T + 1)

# Large: the METEORS Example, n = m = k = 1e5 (POI limits are 3e5).
N = 100_000
case = random_meteors(N, N, N, 10**9, 10**9)
with Timer("METEORS n = m = k = 1e5, Example"):
  res = meteors(N, *case, lib.parallel_bs)
for _ in range(50):
  i = rnd(0, N - 1)
  check_eq(res[i], meteors_one(N, *case, i))
print("OK offline/parallel-binary-search")
