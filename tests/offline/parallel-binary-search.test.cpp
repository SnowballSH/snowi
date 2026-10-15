#include "../testlib.hpp"
#include "../../content/ds/fenwick.hpp"
#include "../../content/offline/parallel-binary-search.hpp"

struct Calls {
  int resets = 0;
  ll applies = 0, preds = 0;
};

// Wraps the callbacks: counts calls and checks that apply(t)
// comes in order 0, 1, ... after each reset and pred reads a
// state of 0..T updates.
vi run(int T, int Q, auto reset, auto apply, auto pred, Calls& c) {
  int next = -1;
  vi res = parallel_bs(
      T, Q, [&] { reset(), c.resets++, next = 0; },
      [&](int t) {
        CHECK_EQ(t, next);
        apply(t), next++, c.applies++;
      },
      [&](int i) {
        CHECK(next >= 0 && 0 <= i && i < Q);
        c.preds++;
        return pred(i);
      });
  CHECK_EQ(sz(res), Q);
  int rounds = 0;
  while ((1 << rounds) < T + 2) rounds++;
  CHECK(c.resets <= rounds); // ceil(log2(T + 2))
  CHECK(c.applies <= (ll)c.resets * T);
  CHECK(c.preds <= (ll)c.resets * Q);
  return res;
}

// METEORS as in the Example: a ring of m sectors, sector j owned
// by state own[j]; shower t adds x[t] to sectors l[t]..r[t].
void stress_meteors(int maxm, int maxT, int maxQ, ll maxx) {
  int m = (int)rnd(1, maxm), T = (int)rnd(0, maxT), Q = (int)rnd(0, maxQ);
  vi own(m), l(T), r(T);
  vector<ll> x(T), need(Q);
  for (int& o : own) o = (int)rnd(0, max(Q - 1, 0));
  rep(t, 0, T)
    l[t] = (int)rnd(0, m - 1), r[t] = (int)rnd(0, m - 1), x[t] = rnd(0, maxx);
  vector<vi> sec(Q);
  if (Q)
    rep(j, 0, m) sec[own[j]].push_back(j);
  rep(i, 0, Q) need[i] = rnd(0, maxx * T * sz(sec[i]) / 4 + 1);
  // brute force: simulate, record each state's first success
  vector<ll> a(m);
  vi want(Q, T + 1);
  rep(t, 0, T + 1) {
    rep(i, 0, Q) {
      ll s = 0;
      for (int j : sec[i]) s += a[j];
      if (want[i] > T && s >= need[i]) want[i] = t;
    }
    if (t == T) break;
    for (int j = l[t];; j = (j + 1) % m) {
      a[j] += x[t];
      if (j == r[t]) break;
    }
  }
  Fenwick fw(m);
  Calls c;
  vi res = run(
      T, Q, [&] { fw = Fenwick(m); },
      [&](int t) {
        fw.add(l[t], x[t]), fw.add(r[t] + 1, -x[t]);
        if (l[t] > r[t]) fw.add(0, x[t]);
      },
      [&](int i) {
        ll s = 0;
        for (int j : sec[i])
          if ((s += fw.query(j + 1)) >= need[i]) break;
        return s >= need[i];
      },
      c);
  CHECK(res == want);
}

// DSU state: the first prefix of edges joining u[i] and v[i].
void stress_dsu(int maxn, int maxT, int maxQ) {
  int n = (int)rnd(1, maxn), T = (int)rnd(0, maxT), Q = (int)rnd(0, maxQ);
  vector<pii> e(T), qs(Q);
  for (auto& [a, b] : e) a = (int)rnd(0, n - 1), b = (int)rnd(0, n - 1);
  for (auto& [a, b] : qs) a = (int)rnd(0, n - 1), b = (int)rnd(0, n - 1);
  vi p(n);
  auto find = [&](int v) {
    while (p[v] != v) v = p[v];
    return v;
  };
  vi want(Q);
  rep(i, 0, Q) {
    iota(all(p), 0);
    auto [u, v] = qs[i];
    int t = 0;
    for (; t <= T; t++) {
      if (find(u) == find(v)) break;
      if (t < T) p[find(e[t].first)] = find(e[t].second);
    }
    want[i] = t; // T + 1 if never
  }
  Calls c;
  vi res = run(
      T, Q, [&] { iota(all(p), 0); },
      [&](int t) { p[find(e[t].first)] = find(e[t].second); },
      [&](int i) { return find(qs[i].first) == find(qs[i].second); }, c);
  CHECK(res == want);
}

int main() {
  { // The Python doctest's METEORS case; Q = 0; T = 0.
    Calls c;
    vi empty = run(5, 0, [] {}, [](int) {}, [](int) { return true; }, c);
    CHECK(empty.empty() && c.resets == 0);
    Calls d;
    vi r = run(0, 3, [] {}, [](int) {}, [](int i) { return i != 1; }, d);
    CHECK(r == vi({0, 1, 0}));
  }
  rep(it, 0, iters(3000)) stress_meteors(6, 10, 5, 4);
  rep(it, 0, iters(300)) stress_meteors(30, 100, 30, 1'000'000'000);
  rep(it, 0, iters(2000)) stress_dsu(7, 12, 8);
  rep(it, 0, iters(200)) stress_dsu(30, 200, 50);
  // Large: the Usage, T = Q = 1e5 point adds and prefix checks.
  elapsed();
  int n = 100'000, T = 100'000, Q = 100'000;
  vi p(T), r(Q);
  vector<ll> w(T), need(Q);
  rep(t, 0, T) p[t] = (int)rnd(0, n - 1), w[t] = rnd(0, 1'000'000'000);
  rep(i, 0, Q) r[i] = (int)rnd(0, n), need[i] = rnd(0, 50'000'000'000'000);
  Fenwick fw(n);
  Calls c;
  vi res = run(
      T, Q, [&] { fw = Fenwick(n); }, [&](int t) { fw.add(p[t], w[t]); },
      [&](int i) { return fw.query(r[i]) >= need[i]; }, c);
  fprintf(stderr, "  T = Q = 1e5: %.2fs, %d rounds\n", elapsed(), c.resets);
  rep(k, 0, 100) {
    int i = (int)rnd(0, Q - 1), t = 0;
    ll s = 0;
    while (t < T && s < need[i]) s += p[t] < r[i] ? w[t] : 0, t++;
    CHECK_EQ(res[i], s >= need[i] ? t : T + 1);
  }
  puts("OK offline/parallel-binary-search");
}
