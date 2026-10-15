#include "../testlib.hpp"
#include "../../content/graph/two-sat.hpp"

struct Instance {
  int n;
  vector<pii> clauses; // (a or b); implications are stored as such
  vector<vi> groups;   // at most one true
  bool holds(const vector<bool>& v) const {
    auto lit = [&](int x) { return x < 0 ? !v[~x] : (bool)v[x]; };
    for (auto [a, b] : clauses)
      if (!lit(a) && !lit(b)) return false;
    for (auto& g : groups) {
      int cnt = 0;
      for (int x : g) cnt += lit(x);
      if (cnt > 1) return false;
    }
    return true;
  }
};

int rnd_lit(int n) {
  int x = (int)rnd(0, n - 1);
  return rnd(0, 1) ? x : ~x;
}

int main() {
  int sat = 0, total = iters(4000);
  rep(it, 0, total) {
    Instance in{(int)rnd(1, 10), {}, {}};
    TwoSat ts(in.n);
    int m = (int)rnd(0, 3 * in.n), groups = (int)rnd(0, 2);
    rep(i, 0, m) {
      int a = rnd_lit(in.n), b = rnd_lit(in.n);
      if (rnd(0, 3) == 0) {
        ts.implies(a, b);
        in.clauses.push_back({~a, b});
      } else {
        ts.add_clause(a, b);
        in.clauses.push_back({a, b});
      }
    }
    rep(i, 0, groups) {
      vi g;
      int len = (int)rnd(0, 5);
      rep(j, 0, len) g.push_back(rnd_lit(in.n));
      ts.at_most_one(g);
      in.groups.push_back(g);
    }
    bool want = false;
    rep(mask, 0, 1 << in.n) {
      vector<bool> v(in.n);
      rep(x, 0, in.n) v[x] = mask >> x & 1;
      if (in.holds(v)) {
        want = true;
        break;
      }
    }
    bool got = ts.solve();
    CHECK_EQ(got, want);
    if (got) {
      sat++;
      CHECK(in.holds(vector<bool>(ts.val.begin(), ts.val.begin() + in.n)));
    }
  }
  CHECK(sat > total / 5 && sat < total * 4 / 5);
  // Large: a chain of implications forced at its start, and one
  // group of 1e5 literals with exactly one forced true.
  int n = 100'000;
  TwoSat ts(n);
  rep(i, 0, n - 1) ts.implies(i, i + 1);
  ts.add_clause(0, 0);
  CHECK(ts.solve());
  CHECK(count(ts.val.begin(), ts.val.end(), true) == n);
  ts.add_clause(~(n - 1), ~(n - 1));
  CHECK(!ts.solve());
  TwoSat amo(n);
  vi all_lits(n);
  iota(all(all_lits), 0);
  amo.at_most_one(all_lits);
  amo.add_clause(n / 2, n / 2);
  CHECK(amo.solve());
  rep(i, 0, n) CHECK_EQ((bool)amo.val[i], i == n / 2);
  amo.add_clause(n / 3, n - 1);
  CHECK(!amo.solve());
  puts("OK graph/two-sat");
}
