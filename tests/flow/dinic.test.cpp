#include "../testlib.hpp"
#include "../../content/flow/dinic.hpp"

struct Arc {
  int u, v;
  ll cap, rev;
};

// Minimum s-t cut by enumerating every source side.
ll brute_cut(int n, const vector<Arc>& arcs, int s, int t) {
  ll best = LLONG_MAX;
  rep(mask, 0, 1 << n) {
    if (!(mask >> s & 1) || (mask >> t & 1)) continue;
    ll c = 0;
    for (auto& a : arcs) {
      bool in_u = mask >> a.u & 1, in_v = mask >> a.v & 1;
      if (in_u && !in_v) c += a.cap;
      if (in_v && !in_u) c += a.rev;
    }
    best = min(best, c);
  }
  return best;
}

// Edmonds-Karp on a capacity matrix (BFS shortest augmenting paths).
ll edmonds_karp(int n, const vector<Arc>& arcs, int s, int t) {
  vector<vector<ll>> c(n, vector<ll>(n));
  for (auto& a : arcs)
    if (a.u != a.v) c[a.u][a.v] += a.cap, c[a.v][a.u] += a.rev;
  ll total = 0;
  while (true) {
    vi par(n, -1);
    par[s] = s;
    queue<int> q;
    q.push(s);
    while (sz(q) && par[t] < 0) {
      int u = q.front();
      q.pop();
      rep(v, 0, n)
        if (par[v] < 0 && c[u][v] > 0) par[v] = u, q.push(v);
    }
    if (par[t] < 0) return total;
    ll f = LLONG_MAX;
    for (int v = t; v != s; v = par[v]) f = min(f, c[par[v]][v]);
    for (int v = t; v != s; v = par[v]) c[par[v]][v] -= f, c[v][par[v]] += f;
    total += f;
  }
}

// Flow bounds, conservation, value, and a cut of equal capacity.
void certify(Dinic& d, int n, const vector<Arc>& arcs, const vi& ids, int s,
             int t, ll value) {
  vector<ll> net(n);
  ll cut = 0;
  rep(i, 0, sz(arcs)) {
    auto& a = arcs[i];
    ll f = d.flow(ids[i]);
    CHECK(-a.rev <= f && f <= a.cap);
    net[a.u] -= f, net[a.v] += f;
    bool su = d.source_side(a.u), sv = d.source_side(a.v);
    if (su && !sv) {
      cut += a.cap;
      CHECK_EQ(f, a.cap);
    }
    if (sv && !su) {
      cut += a.rev;
      CHECK_EQ(f, -a.rev);
    }
  }
  rep(v, 0, n)
    if (v != s && v != t) CHECK_EQ(net[v], 0LL);
  CHECK_EQ(net[t], value);
  CHECK_EQ(-net[s], value);
  CHECK(d.source_side(s));
  CHECK(!d.source_side(t));
  CHECK_EQ(cut, value);
}

vector<Arc> random_arcs(int n, int m, ll maxc) {
  vector<Arc> arcs;
  rep(i, 0, m) {
    Arc a{(int)rnd(0, n - 1), (int)rnd(0, n - 1), rnd(0, maxc), 0};
    int kind = (int)rnd(0, 3);
    if (kind == 1) a.rev = a.cap;
    if (kind == 2) a.rev = rnd(0, maxc);
    arcs.push_back(a);
  }
  return arcs;
}

pair<Dinic, vi> build(int n, const vector<Arc>& arcs) {
  Dinic d(n);
  vi ids;
  for (auto& a : arcs) ids.push_back(d.add_edge(a.u, a.v, a.cap, a.rev));
  return {d, ids};
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(2, 8), m = (int)rnd(0, 16);
    ll maxc = rnd(0, 3) ? rnd(1, 10) : (ll)1e15;
    auto arcs = random_arcs(n, m, maxc);
    int s = (int)rnd(0, n - 1), t = (int)rnd(0, n - 2);
    if (t >= s) t++;
    auto [d, ids] = build(n, arcs);
    ll f = d.max_flow(s, t);
    CHECK_EQ(f, brute_cut(n, arcs, s, t));
    certify(d, n, arcs, ids, s, t, f);
    CHECK_EQ(d.max_flow(s, t), 0LL); // already maximum
    certify(d, n, arcs, ids, s, t, f);
    // Tweak: reset and run again on another pair.
    for (auto& adj : d.g)
      for (auto& x : adj) x.cap = x.orig;
    swap(s, t);
    f = d.max_flow(s, t);
    CHECK_EQ(f, brute_cut(n, arcs, s, t));
    certify(d, n, arcs, ids, s, t, f);
  }
  rep(it, 0, iters(300)) {
    int n = (int)rnd(2, 30), m = (int)rnd(0, 120);
    auto arcs = random_arcs(n, m, rnd(1, 1000));
    auto [d, ids] = build(n, arcs);
    ll f = d.max_flow(0, n - 1);
    CHECK_EQ(f, edmonds_karp(n, arcs, 0, n - 1));
    certify(d, n, arcs, ids, 0, n - 1, f);
  }
  // Extremes: capacities near 2^62 on parallel paths.
  {
    ll big = 1LL << 61;
    vector<Arc> arcs = {
        {0, 1, big, 0}, {1, 3, big, big}, {0, 2, big, 0}, {2, 3, big, 0}};
    auto [d, ids] = build(4, arcs);
    CHECK_EQ(d.max_flow(0, 3), 2 * big);
    certify(d, 4, arcs, ids, 0, 3, 2 * big);
  }
  // Diamond chain of 2 x 60 vertices ending in a dead end, plus a
  // separate path to t of the same length. Without the current-edge
  // pointer the dead end is explored once per path: 2^60 times.
  {
    int L = 60, n = 3 * L + 2, s = 3 * L, t = 3 * L + 1;
    vector<Arc> arcs = {{s, 0, 10, 0}, {s, 1, 10, 0}, {s, 2 * L, 1, 0}};
    rep(i, 0, L - 1) {
      for (int a : {2 * i, 2 * i + 1})
        arcs.push_back({a, 2 * i + 2, 10, 0}),
            arcs.push_back({a, 2 * i + 3, 10, 0});
      arcs.push_back({2 * L + i, 2 * L + i + 1, 1, 0});
    }
    arcs.push_back({3 * L - 1, t, 1, 0});
    auto [d, ids] = build(n, arcs);
    CHECK_EQ(d.max_flow(s, t), 1LL);
    certify(d, n, arcs, ids, s, t, 1);
  }
  // Large: unit bipartite matching network, then a random
  // network with large capacities, then a long path.
  elapsed();
  {
    int side = 20000, n = 2 * side + 2, s = n - 2, t = n - 1;
    vector<Arc> arcs;
    rep(i, 0, side) arcs.push_back({s, i, 1, 0});
    rep(i, 0, side) arcs.push_back({side + i, t, 1, 0});
    rep(i, 0, 60000)
      arcs.push_back(
          {(int)rnd(0, side - 1), side + (int)rnd(0, side - 1), 1, 0});
    auto [d, ids] = build(n, arcs);
    ll f = d.max_flow(s, t);
    CHECK(f > side * 9 / 10);
    certify(d, n, arcs, ids, s, t, f);
  }
  {
    int n = 1000;
    auto arcs = random_arcs(n, 20000, (ll)1e12);
    auto [d, ids] = build(n, arcs);
    ll f = d.max_flow(0, n - 1);
    certify(d, n, arcs, ids, 0, n - 1, f);
  }
  {
    int n = 20000;
    vector<Arc> arcs;
    rep(i, 0, n - 1) arcs.push_back({i, i + 1, rnd(5, 9), 0});
    rep(i, 0, n - 3) arcs.push_back({i, i + 3, 1, 0});
    auto [d, ids] = build(n, arcs);
    ll f = d.max_flow(0, n - 1);
    certify(d, n, arcs, ids, 0, n - 1, f);
  }
  fprintf(stderr, "  large cases: %.2fs\n", elapsed());
  puts("OK flow/dinic");
}
