#include "../testlib.hpp"
#include "../../content/offline/dynamic-connectivity.hpp"

// Component label of every vertex: the smallest vertex in it.
vi labels(int n, auto find) {
  vi root(n), low(n, INT_MAX), res(n);
  rep(v, 0, n) root[v] = find(v), low[root[v]] = min(low[root[v]], v);
  rep(v, 0, n) res[v] = low[root[v]];
  return res;
}

vi brute_labels(int n, const vector<pii>& edges) {
  vi p(n);
  iota(all(p), 0);
  auto find = [&](int x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  for (auto [u, v] : edges) p[find(u)] = find(v);
  return labels(n, find);
}

// Random ops: add (loops and parallel edges allowed), remove a
// present edge, or query. Checks every query's partition and
// component count, and that answers come once each, in order.
void stress(int n, int ops, double p_add, double p_query, int samples = -1) {
  DynConn dc(n);
  vector<pii> present;
  vector<vi> want;
  rep(t, 0, ops) {
    double x = rnd_real(0, 1);
    if (x < p_query) {
      CHECK_EQ(dc.query(), sz(want));
      bool check = samples < 0 || rnd_real(0, ops * p_query) < samples;
      want.push_back(check ? brute_labels(n, present) : vi());
    } else if (x < p_query + p_add || present.empty()) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      dc.add(u, v), present.push_back({u, v});
    } else {
      int k = (int)rnd(0, sz(present) - 1);
      swap(present[k], present.back());
      auto [u, v] = present.back();
      present.pop_back();
      rnd(0, 1) ? dc.remove(u, v) : dc.remove(v, u);
    }
  }
  int next = 0;
  dc.solve([&](int i) {
    CHECK_EQ(i, next++);
    if (want[i].empty()) return;
    vi got = labels(n, [&](int v) { return dc.d.find(v); });
    CHECK(got == want[i]);
    CHECK_EQ(n - dc.d.time() / 2, sz(set<int>(all(want[i]))));
  });
  CHECK_EQ(next, sz(want));
  CHECK_EQ(dc.d.time(), 0); // everything rolled back
}

int main() {
  { // the packet's example, 0-indexed
    DynConn dc(4);
    dc.add(0, 1), dc.add(1, 2), dc.query();
    dc.remove(0, 1), dc.add(2, 3), dc.query();
    vi comps, same;
    dc.solve([&](int) {
      comps.push_back(4 - dc.d.time() / 2);
      same.push_back(dc.d.find(0) == dc.d.find(2));
    });
    CHECK(comps == vi({2, 2}) && same == vi({1, 0}));
  }
  { // no queries; queries on an empty graph; one vertex
    DynConn a(3);
    a.add(0, 1), a.remove(1, 0), a.add(1, 2);
    a.solve([&](int) { CHECK(false); });
    DynConn b(1);
    b.query(), b.add(0, 0), b.query(), b.remove(0, 0), b.query();
    int calls = 0;
    b.solve([&](int) {
      calls++;
      CHECK_EQ(b.d.time(), 0);
    });
    CHECK_EQ(calls, 3);
  }
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 7), ops = (int)rnd(0, 40);
    stress(n, ops, rnd_real(0.2, 0.6), rnd_real(0.1, 0.5));
  }
  rep(it, 0, iters(200)) stress(30, 300, 0.4, 0.3);
  // Large: 1e5 operations, a few queries checked (each check
  // costs O(n log n)).
  elapsed();
  stress(100'000, 100'000, 0.45, 0.2, 5);
  stress(1000, 100'000, 0.4, 0.3, 50);
  fprintf(stderr, "  2 x 1e5 operations: %.2fs\n", elapsed());
  puts("OK offline/dynamic-connectivity");
}
