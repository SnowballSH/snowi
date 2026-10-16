from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("geometry/manhattan-mst")


def dist(a, b):
  return abs(a[0] - b[0]) + abs(a[1] - b[1])


def brute(ps):
  """Prim on the complete graph, O(n^2)."""
  n = len(ps)
  d, done, total = [float("inf")] * n, [False] * n, 0
  d[0] = 0
  for _ in range(n):
    u = min((v for v in range(n) if not done[v]), key=d.__getitem__)
    done[u] = True
    total += d[u]
    for v in range(n):
      d[v] = min(d[v], dist(ps[u], ps[v]))
  return total


def check_points(ps):
  n = len(ps)
  total, tree = lib.manhattan_mst(ps)
  check_eq(len(tree), max(n - 1, 0))
  par = list(range(n))

  def find(x):
    while par[x] != x:
      x = par[x]
    return x

  w = 0
  for i, j in tree:  # a spanning tree whose weight is total
    check(find(i) != find(j), f"cycle in {ps}")
    par[find(i)] = find(j)
    w += dist(ps[i], ps[j])
  check_eq(w, total)
  check_eq(total, brute(ps), f"{ps}")


check_points([(5, 5)])
check_points([(0, 0), (0, 0)])
check_points([(0, 0), (1, 1), (2, 2), (3, 3)])
check_points([(0, 0), (1, -1), (2, -2), (3, -3)])
for it in range(iters(2000)):
  n = rnd(1, 25)
  c = (1, 3, 10, 1000, 1 << 30)[it % 5]
  ps = [(rnd(-c, c), rnd(-c, c)) for _ in range(n)]
  if it % 7 == 0:  # points on diagonals and axes: octant borders
    ps = [(x, (x, -x, 0)[rnd(0, 2)]) for x, _ in ps]
  check_points(ps)

n = 200_000
c = 1 << 30
ps = [(rnd(-c, c), rnd(-c, c)) for _ in range(n)]
with Timer("manhattan_mst of 2e5 points"):
  total, tree = lib.manhattan_mst(ps)
check_eq(len(tree), n - 1)
print("OK geometry/manhattan-mst")
