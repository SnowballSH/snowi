from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("flow/hopcroft-karp")
HopcroftKarp = lib.HopcroftKarp


def brute(n, m, adj):
  """Maximum matching by DP over left vertices and used right sets."""
  f = [0] * (1 << m)
  for u in reversed(range(n)):
    nf = f[:]
    for mask in range(1 << m):
      for v in adj[u]:
        if not mask >> v & 1:
          nf[mask] = max(nf[mask], 1 + f[mask | 1 << v])
    f = nf
  return f[0]


def certify(hk, n, m, edges, k):
  """A valid matching of size k plus a vertex cover of size k."""
  es = set(edges)
  matched = 0
  for u in range(n):
    if hk.ml[u] >= 0:
      check((u, hk.ml[u]) in es)
      check_eq(hk.mr[hk.ml[u]], u)
      matched += 1
  for v in range(m):
    if hk.mr[v] >= 0:
      check_eq(hk.ml[hk.mr[v]], v)
  check_eq(matched, k)
  cl, cr = hk.cover()
  check_eq(len(cl) + len(cr), k)
  inl, inr = set(cl), set(cr)
  for u, v in edges:
    check(u in inl or v in inr)


for _ in range(iters(3000)):
  n, m = rnd(0, 7), rnd(0, 7)
  e = rnd(0, 18) if n and m else 0
  edges = [(rnd(0, n - 1), rnd(0, m - 1)) for _ in range(e)]
  adj = [[] for _ in range(n)]
  hk = HopcroftKarp(n, m)
  for u, v in edges:
    adj[u].append(v)
    hk.add_edge(u, v)
  k = hk.max_matching()
  check_eq(k, brute(n, m, adj))
  certify(hk, n, m, edges, k)
  check_eq(hk.max_matching(), k)  # nothing left to augment
  certify(hk, n, m, edges, k)
  if n and m:  # Tweak: add edges later, continue
    for _ in range(3):
      u, v = rnd(0, n - 1), rnd(0, m - 1)
      edges.append((u, v))
      adj[u].append(v)
      hk.add_edge(u, v)
    k = hk.max_matching()
    check_eq(k, brute(n, m, adj))
    certify(hk, n, m, edges, k)

# Diamond: left A_i, B_i matched to a_i, b_i, both joined to
# a_(i+1) and b_(i+1); a free root reaches a_0, b_0 and the last
# layer is a dead end. A chain C_j - c_j of length L + 1 sets the
# phase length. Without visit marks and the persistent it[]
# pointers the search would repeat 2^L times.
L = 60
n, m, r = 3 * L + 3, 3 * L + 2, 3 * L + 1
edges = []
for i in range(L):
  for u in (2 * i, 2 * i + 1):
    edges.append((u, u))
    if i + 1 < L:
      edges += [(u, 2 * i + 2), (u, 2 * i + 3)]
for j in range(2 * L, 3 * L + 1):
  edges += [(j, j), (j, j + 1)]
edges += [(r, 0), (r, 1), (r + 1, 2 * L)]
hk = HopcroftKarp(n, m)
for u, v in edges:
  hk.add_edge(u, v)
check_eq(hk.max_matching(), m)
certify(hk, n, m, edges, m)

n, m = 200_000, 150_000
edges = [(rnd(0, n - 1), rnd(0, m - 1)) for _ in range(400_000)]
hk = HopcroftKarp(n, m)
for u, v in edges:
  hk.add_edge(u, v)
with Timer("n = 2e5, m = 1.5e5, 4e5 random edges"):
  k = hk.max_matching()
certify(hk, n, m, edges, k)

n = 200_000  # chain: the last augmenting path has length 2n
edges = []
hk = HopcroftKarp(n, n)
for u in range(n):
  if u + 1 < n:
    edges.append((u, u + 1))
  edges.append((u, u))
for u, v in edges:
  hk.add_edge(u, v)
with Timer("chain of 2e5, one augmenting path of length 4e5"):
  check_eq(hk.max_matching(), n)
certify(hk, n, n, edges, n)
print("OK flow/hopcroft-karp")
