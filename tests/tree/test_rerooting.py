from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("tree/rerooting")
MOD = 998244353


def random_tree(n, shape):
  """Shapes: 0 random, 1 path, 2 star, 3 long thin, 4 binary-ish."""
  perm = list(range(n))
  rng.shuffle(perm)
  adj = [[] for _ in range(n)]
  for i in range(1, n):
    if shape == 1:
      p = i - 1
    elif shape == 2:
      p = 0
    elif shape == 3:
      p = rnd(max(0, i - 3), i - 1)
    elif shape == 4:
      p = (i - 1) // 2
    else:
      p = rnd(0, i - 1)
    adj[perm[i]].append(perm[p])
    adj[perm[p]].append(perm[i])
  for a in adj:
    rng.shuffle(a)
  return adj


def bfs_dist(adj, s):
  d, q = [-1] * len(adj), [s]
  d[s] = 0
  for u in q:
    for v in adj[u]:
      if d[v] < 0:
        d[v] = d[u] + 1
        q.append(v)
  return d


def rooted_dp(adj, r, e, merge, lift, fin):
  """The rooted DP from scratch: sub[v] = value of v's subtree."""
  n = len(adj)
  par, order = [-1] * n, [r]
  for u in order:
    for v in adj[u]:
      if v != par[u]:
        par[v] = u
        order.append(v)
  sub = [None] * n
  for v in reversed(order):
    acc = e
    for c in adj[v]:
      if c != par[v]:
        acc = merge(acc, lift(sub[c], c, v))
    sub[v] = fin(acc, v)
  return sub


def check_all(adj, e, merge, lift, fin):
  """Every field of Reroot against the DP from every root."""
  n = len(adj)
  rr = lib.Reroot(adj, e, merge, lift, fin)
  sub = [rooted_dp(adj, r, e, merge, lift, fin) for r in range(n)]
  d0 = bfs_dist(adj, 0)
  check_eq(rr.par[0], -1)
  check_eq(rr.dn[0], e)
  check_eq(rr.up[0], e)
  for v in range(n):
    check_eq(rr.dp[v], sub[v][v])
    if v:
      p = rr.par[v]
      check(p in adj[v] and d0[p] == d0[v] - 1)
      check_eq(rr.dn[v], lift(sub[0][v], v, p))
      check_eq(rr.up[v], lift(sub[v][p], p, v))
  return rr


def merge_pair(a, b):
  return a[0] + b[0], a[1] + b[1]


def lift_pair(x, u, v):
  return x[0], x[1] + x[0]


def fin_pair(x, v):
  return x[0] + 1, x[1]


def check_modular(adj, m):
  """No inverse; depends on both endpoints and the direction."""
  n = len(adj)
  A = [rnd(0, 3) for _ in range(n)]
  B = [rnd(0, 3) for _ in range(n)]
  C = [rnd(0, m - 1) for _ in range(n)]
  D = [rnd(0, 2) for _ in range(n)]
  check_all(
    adj,
    1,
    lambda a, b: a * b % m,
    lambda x, u, v: (x * A[u] + B[v] + D[u]) % m,
    lambda acc, v: (acc * D[v] + C[v]) % m,
  )


def check_tree(adj):
  n = len(adj)
  dist = [bfs_dist(adj, v) for v in range(n)]
  # Usage: eccentricity.
  ecc = check_all(adj, 0, max, lambda x, u, v: x + 1, lambda a, v: a)
  check_eq(ecc.dp, [max(d) for d in dist])
  # Example: sum of distances, subtree sizes.
  rr = check_all(adj, (0, 0), merge_pair, lift_pair, fin_pair)
  check_eq(rr.dp, [(n, sum(d)) for d in dist])
  for v in range(1, n):
    below = sum(dist[0][x] == dist[0][v] + dist[v][x] for x in range(n))
    check_eq(rr.dn[v][0], below)
  # Example: connected vertex sets containing v (all subsets).
  cs = check_all(
    adj,
    1,
    lambda a, b: a * b % MOD,
    lambda x, u, v: (x + 1) % MOD,
    lambda x, v: x,
  )
  if n <= 9:
    want = [0] * n
    for mask in range(1, 1 << n):
      s = (mask & -mask).bit_length() - 1
      seen, st = 1 << s, [s]
      while st:
        u = st.pop()
        for v in adj[u]:
          if mask >> v & 1 and not seen >> v & 1:
            seen |= 1 << v
            st.append(v)
      if seen == mask:
        for v in range(n):
          want[v] += mask >> v & 1
    check_eq(cs.dp, want)
  check_modular(adj, 7)
  check_modular(adj, MOD)
  # Sorted tuples (multisets) catch a child merged twice or skipped.
  check_all(
    adj,
    (),
    lambda a, b: tuple(sorted(a + b)),
    lambda x, u, v: tuple(sorted(x + (1000 * u + v,))),
    lambda x, v: tuple(sorted(x + (-v - 1,))),
  )


check_tree([[]])
check_tree([[1], [0]])
for it in range(iters(800)):
  n = rnd(1, 6 if it < 200 else 13)
  check_tree(random_tree(n, it % 5))

# Large: n = 2e5 (random, path, star), independent formulas.
n = 200_000
for shape in range(3):
  adj = random_tree(n, shape)
  with Timer(f"shape {shape}, n = 2e5: eccentricity (ints)"):
    ecc = lib.Reroot(adj, 0, max, lambda x, u, v: x + 1, lambda a, v: a)
  with Timer(f"shape {shape}, n = 2e5: sum of distances (tuples)"):
    rr = lib.Reroot(adj, (0, 0), merge_pair, lift_pair, fin_pair)
  d0 = bfs_dist(adj, 0)
  order = sorted(range(n), key=d0.__getitem__)
  par, size = [-1] * n, [1] * n
  for v in order:
    for u in adj[v]:
      if d0[u] == d0[v] + 1:
        par[u] = v
  for v in reversed(order[1:]):
    size[par[v]] += size[v]
  ans = [0] * n
  ans[0] = sum(d0)
  for v in order[1:]:
    ans[v] = ans[par[v]] + n - 2 * size[v]
  check_eq(rr.dp, [(n, a) for a in ans])
  a = max(range(n), key=d0.__getitem__)
  da = bfs_dist(adj, a)
  b = max(range(n), key=da.__getitem__)
  db = bfs_dist(adj, b)
  check_eq(ecc.dp, [max(x, y) for x, y in zip(da, db)])
print("OK tree/rerooting")
