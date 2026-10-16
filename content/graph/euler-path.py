"""
Title: Euler path
Description: Hierholzer's algorithm: a walk from `src` that uses
  every edge exactly once. Returns the vertices `vs` ($m + 1$ of
  them, ending where the walk ends) and the edge ids `es`, where edge
  `es[i]` goes from `vs[i]` to `vs[i + 1]`. Returns `None` if no
  such walk exists.
Precondition: Edges `(u, v)` with ids `0..m-1`, directed `u -> v` or
  undirected; vertices 0-indexed. Self-loops and parallel edges are
  fine; isolated vertices are ignored.
Time: $O(n + m)$, iterative.
Usage:
  res = euler_path(n, edges, directed, src)
  if res: vs, es = res  # None: impossible from src
Tweak: Circuit only: also require `vs[-1] == src`.
Example:
  >>> edges = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 2)]
  >>> euler_path(4, edges, False, 0)  # all degrees even
  ([0, 1, 2, 3, 2, 0], [0, 1, 3, 4, 2])
  >>> euler_path(4, edges, True, 0)  # directed: also a circuit
  ([0, 1, 2, 3, 2, 0], [0, 1, 3, 4, 2])
  >>> edges.append((0, 3))  # now 0 and 3 have odd degree
  >>> euler_path(4, edges, False, 1) is None
  True
  >>> vs, es = euler_path(4, edges, False, 3); vs
  [3, 2, 1, 0, 2, 3, 0]
  >>> # De Bruijn sequence: every length-3 binary word once
  >>> E = [(w, (2 * w + c) % 4) for w in range(4) for c in (0, 1)]
  >>> vs, es = euler_path(4, E, True, 0)
  >>> "00" + "".join(str(e % 2) for e in es)
  '0001011100'
Notes: Undirected: a path exists iff the edges are connected and 0
  or 2 vertices have odd degree (a circuit iff 0); it must start at
  an odd vertex if there are two. Directed: the edges are connected
  (ignoring directions) and every vertex has $"in" = "out"$
  (circuit), or the start has $"out" = "in" + 1$, the end
  $"in" = "out" + 1$ and all others balance. Think of it for "use
  every edge/domino/word once", "reconstruct a string from its
  overlapping pieces" (pieces are edges), De Bruijn sequences
  (Example: vertices are words of length $k - 1$, edge `w * a + c`
  appends letter `c`), and "draw without lifting the pen". Covering
  all edges of an undirected component with $2k$ odd vertices needs
  $max(1, k)$ paths: pair odd vertices with dummy edges, run a
  circuit, cut at the dummies. Orienting an undirected graph so that
  every vertex has $|"in" - "out"| <= 1$: in each component, pair
  its odd vertices with dummy edges and orient along an Euler
  circuit (one call per component: a walk never leaves its
  component).
Source: kactl EulerWalk (CC0), with explicit degree checks.
Status: stress-tested against exhaustive search
"""


def euler_path(n, edges, directed, src):
  m = len(edges)
  adj = [[] for _ in range(n)]
  other, deg = [0] * m, [0] * n  # far end of e: other[e] ^ u
  for e, (u, v) in enumerate(edges):
    other[e] = u ^ v
    adj[u].append(e)
    deg[u] += 1
    if directed:
      deg[v] -= 1
    else:
      adj[v].append(e)
      deg[v] += 1
  # Pretend the walk enters src once more: then at most one
  # vertex (the end) may be unbalanced.
  deg[src] += -1 if directed else 1
  if sum(d != 0 if directed else d & 1 for d in deg) > 1:
    return None
  used, it = [False] * m, [0] * n
  vs, es, sv, se = [], [], [src], [-1]
  while sv:
    u = sv[-1]
    out, i = adj[u], it[u]
    while i < len(out) and used[out[i]]:
      i += 1
    if i == len(out):
      it[u] = i
      vs.append(sv.pop())
      es.append(se.pop())
    else:
      e = out[i]
      it[u] = i + 1
      used[e] = True
      sv.append(other[e] ^ u)
      se.append(e)
  if len(vs) != m + 1:  # some edge is not connected to src
    return None
  es.pop()
  return vs[::-1], es[::-1]
