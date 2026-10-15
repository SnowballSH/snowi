from itertools import product

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("graph/two-sat")


def holds(v, clauses, groups):
  def lit(x):
    return not v[~x] if x < 0 else v[x]

  if any(not lit(a) and not lit(b) for a, b in clauses):
    return False
  return all(sum(map(lit, g)) <= 1 for g in groups)


def rnd_lit(n):
  x = rnd(0, n - 1)
  return x if rnd(0, 1) else ~x


sat, total = 0, iters(3000)
for _ in range(total):
  n = rnd(1, 9)
  ts = lib.TwoSat(n)
  clauses, groups = [], []
  for _ in range(rnd(0, 3 * n)):
    a, b = rnd_lit(n), rnd_lit(n)
    if rnd(0, 3) == 0:
      ts.implies(a, b)
      clauses.append((~a, b))
    else:
      ts.add_clause(a, b)
      clauses.append((a, b))
  for _ in range(rnd(0, 2)):
    g = [rnd_lit(n) for _ in range(rnd(0, 5))]
    ts.at_most_one(g)
    groups.append(g)
  want = any(
    holds(v, clauses, groups) for v in product([False, True], repeat=n)
  )
  val = ts.solve()
  check_eq(val is not None, want)
  if val is not None:
    sat += 1
    check(holds(val[:n], clauses, groups))
check(total // 5 < sat < total * 4 // 5, f"{sat} of {total} satisfiable")

# A chain of implications forced at its start; one big group.
n = 200_000
with Timer("2e5 variables, 2e5 implications"):
  ts = lib.TwoSat(n)
  for i in range(n - 1):
    ts.implies(i, i + 1)
  ts.add_clause(0, 0)
  check_eq(ts.solve(), [True] * n)
ts.add_clause(~(n - 1), ~(n - 1))
check(ts.solve() is None)
with Timer("at_most_one over 2e5 literals"):
  amo = lib.TwoSat(n)
  amo.at_most_one(list(range(n)))
  amo.add_clause(n // 2, n // 2)
  val = amo.solve()
check_eq(val[:n], [i == n // 2 for i in range(n)])
amo.add_clause(n // 3, n - 1)
check(amo.solve() is None)
print("OK graph/two-sat")
