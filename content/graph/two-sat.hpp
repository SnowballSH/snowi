/**
 * Title: 2-SAT
 * Description: Decides whether boolean variables `0..n-1` can
 *   satisfy all clauses $(a or b)$, and finds such values. A literal
 *   is `x` for "x is true" or `~x` for "x is false". `solve()`
 *   returns false if impossible; otherwise `val[x]` is the value of
 *   `x`.
 * Precondition: Literals use variables below `n` (or ones made by
 *   `add_var`).
 * Time: $O(n + m)$ for $m$ clauses; `at_most_one` of $k$ literals
 *   adds $k - 2$ variables and $3k - 5$ clauses.
 * Usage:
 *   TwoSat ts(n);
 *   ts.add_clause(a, ~b); // a or not b
 *   ts.implies(a, b); ts.add_clause(a, a); // force a
 *   ts.at_most_one({a, ~b, c});
 *   if (ts.solve()) ... ts.val[x] ...
 * Example:
 *   // n flags, flag i at x[i][0] or x[i][1]; any two flags at
 *   // least D apart. Variable i: flag i uses x[i][1].
 *   TwoSat ts(n);
 *   rep(i, 0, n) rep(j, 0, i) rep(p, 0, 2) rep(q, 0, 2)
 *     if (abs(x[i][p] - x[j][q]) < D) // not both choices
 *       ts.add_clause(p ? ~i : i, q ? ~j : j);
 *   if (ts.solve()) rep(i, 0, n) cout << x[i][ts.val[i]] << '\n';
 * Notes: Think 2-SAT when every object has exactly two options and
 *   every constraint involves at most two objects ("if A picks X
 *   then B must not pick Y"). "At least one of three" is not 2-SAT:
 *   look for a reason it is automatic or guess one element. A
 *   global constraint (one clique, one special vertex) often becomes
 *   pairwise after fixing one element: try each. Minimize a
 *   threshold: binary search plus 2-SAT. Equal/different
 *   constraints only: a parity DSU also counts solutions
 *   ($2^"components"$). Counting 2-SAT solutions in general is
 *   #P-hard. Forced literal `x`: `~x` reaches `x`. "At most one per
 *   group" with pairwise clauses is quadratic: use `at_most_one`.
 *   With $10^6$ clauses, flat edge arrays beat `vector<vi>`.
 * Uses: graph/scc
 * Source: kactl 2-SAT (CC0); UCup S4-S5 Training Packet 18.
 * Status: stress-tested against exhaustive search
 */
#pragma once
#include "scc.hpp"

struct TwoSat {
  int n;
  vector<vi> g;
  vector<bool> val;
  TwoSat(int vars = 0) : n(vars), g(2 * vars) {}
  static int node(int x) { return x < 0 ? 2 * ~x + 1 : 2 * x; }
  int add_var() { return g.resize(2 * ++n), n - 1; }
  void add_clause(int a, int b) {
    g[node(~a)].push_back(node(b));
    g[node(~b)].push_back(node(a));
  }
  void implies(int a, int b) { add_clause(~a, b); }
  void at_most_one(const vi& lits) {
    if (sz(lits) < 2) return;
    int pre = lits[0]; // "some of lits[0..i-1] is true"
    rep(i, 1, sz(lits)) {
      add_clause(~pre, ~lits[i]);
      if (i + 1 == sz(lits)) break;
      int next = add_var();
      implies(pre, next), implies(lits[i], next), pre = next;
    }
  }
  bool solve() {
    auto [comp, k] = scc(g);
    val.assign(n, false);
    rep(x, 0, n) {
      if (comp[2 * x] == comp[2 * x + 1]) return false;
      val[x] = comp[2 * x] < comp[2 * x + 1];
    }
    return true;
  }
};
