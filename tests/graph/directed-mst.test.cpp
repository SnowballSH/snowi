#include "../testlib.hpp"
#include "../../content/graph/directed-mst.hpp"

typedef array<ll, 3> E; // a, b, w
const ll BIG = LLONG_MAX / 4;

// Every way to give each non-root vertex one incoming edge; keep
// those where following parents always ends at the root.
optional<ll> brute(int n, int root, const vector<E>& es) {
  vector<vector<E>> in(n);
  for (auto& e : es)
    if (e[1] != root && e[0] != e[1]) in[e[1]].push_back(e);
  optional<ll> best;
  vi choice(n), par(n, -1);
  auto rec = [&](auto&& self, int v, ll cost) -> void {
    if (v == n) {
      rep(s, 0, n) {
        int x = s, steps = 0;
        while (x != root && steps++ <= n) x = par[x];
        if (x != root) return;
      }
      if (!best || cost < *best) best = cost;
      return;
    }
    if (v == root) return self(self, v + 1, cost);
    for (auto& e : in[v]) par[v] = (int)e[0], self(self, v + 1, cost + e[2]);
  };
  rec(rec, 0, 0);
  return best;
}

// Textbook O(n m) Chu-Liu/Edmonds with explicit contraction.
optional<ll> slow(int n, int root, vector<E> es) {
  ll res = 0;
  while (true) {
    vector<ll> in(n, BIG);
    vi pre(n, -1), id(n, -1), vis(n, -1);
    for (auto [a, b, w] : es)
      if (a != b && w < in[b]) in[b] = w, pre[b] = (int)a;
    in[root] = 0;
    rep(v, 0, n)
      if (in[v] == BIG) return nullopt;
    int cnt = 0;
    rep(v, 0, n) {
      res += in[v];
      int x = v;
      while (vis[x] != v && id[x] < 0 && x != root) vis[x] = v, x = pre[x];
      if (x != root && id[x] < 0) {
        for (int y = pre[x]; y != x; y = pre[y]) id[y] = cnt;
        id[x] = cnt++;
      }
    }
    if (!cnt) return res;
    rep(v, 0, n)
      if (id[v] < 0) id[v] = cnt++;
    for (auto& [a, b, w] : es) {
      if (id[a] != id[b]) w -= in[b];
      a = id[a], b = id[b];
    }
    n = cnt, root = id[root];
  }
}

// par must be an arborescence of existing edges costing `cost`.
void check_tree(int n, int root, const vector<E>& es, ll cost, const vi& par) {
  CHECK_EQ(sz(par), n);
  CHECK_EQ(par[root], -1);
  map<pii, ll> cheapest;
  for (auto [a, b, w] : es) {
    pii key{(int)a, (int)b};
    if (!cheapest.count(key) || w < cheapest[key]) cheapest[key] = w;
  }
  ll sum = 0;
  rep(v, 0, n)
    if (v != root) {
      CHECK(cheapest.count({par[v], v}));
      sum += cheapest[{par[v], v}];
    }
  CHECK_EQ(sum, cost);
  vi state(n); // 0 new, 1 on stack, 2 reaches root
  state[root] = 2;
  rep(s, 0, n) {
    vi st;
    int x = s;
    while (state[x] == 0) state[x] = 1, st.push_back(x), x = par[x];
    CHECK_EQ(state[x], 2);
    for (int y : st) state[y] = 2;
  }
}

void run(int n, int root, const vector<E>& es, optional<ll> want) {
  DirectedMST g(n);
  for (auto [a, b, w] : es) g.add_edge((int)a, (int)b, w);
  rep(rep_, 0, 2) { // solve() twice: it must not change the graph
    auto [cost, par] = g.solve(root);
    if (!want) {
      CHECK_EQ(cost, -1LL);
      CHECK(par.empty());
    } else {
      CHECK_EQ(cost, *want);
      check_tree(n, root, es, cost, par);
    }
  }
}

vector<E> random_edges(int n, int m, ll lo, ll hi) {
  vector<E> es(m);
  for (auto& e : es) e = {rnd(0, n - 1), rnd(0, n - 1), rnd(lo, hi)};
  return es;
}

// Weights at the Precondition's limit, n * max|w| = 4e18, mostly
// near the extremes; UBSan flags any overflow.
vector<E> extreme_edges(int n, int m) {
  ll big = 4'000'000'000'000'000'000 / n;
  vector<E> es(m);
  for (auto& e : es) {
    ll w = rnd(0, 1) ? big - rnd(0, 2) : rnd(0, big);
    e = {rnd(0, n - 1), rnd(0, n - 1), rnd(0, 1) ? w : -w};
  }
  return es;
}

// Free root (Tweak): super vertex n with edges n -> v of weight
// B = 1 + sum |w|; a single root works iff exactly one par[v] == n.
void free_root(int n, const vector<E>& es) {
  ll B = 1;
  for (auto& e : es) B += abs(e[2]);
  DirectedMST g(n + 1);
  for (auto [a, b, w] : es) g.add_edge((int)a, (int)b, w);
  rep(v, 0, n) g.add_edge(n, v, B);
  auto [cost, par] = g.solve(n);
  optional<ll> want;
  rep(r, 0, n) {
    auto c = brute(n, r, es);
    if (c && (!want || *c < *want)) want = c;
  }
  if (count(all(par), n) == 1) CHECK(want == cost - B);
  else CHECK(!want);
}

int main() {
  run(1, 0, {}, 0);
  run(1, 0, {{0, 0, -5}}, 0);
  run(2, 0, {}, nullopt);
  run(2, 1, {{0, 1, 3}}, nullopt);
  run(2, 0, {{0, 1, -1}}, -1);                // cost -1 with a valid answer
  run(3, 0, {{1, 2, 1}, {2, 1, 1}}, nullopt); // unreachable cycle
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(1, 7), m = (int)rnd(0, 3 * n);
    bool neg = rnd(0, 1);
    auto es = random_edges(n, m, neg ? -6 : 0, 10);
    int root = (int)rnd(0, n - 1);
    auto want = brute(n, root, es);
    CHECK(want == slow(n, root, es));
    run(n, root, es, want);
  }
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(1, 7), m = (int)rnd(0, 3 * n);
    auto es = extreme_edges(n, m);
    int root = (int)rnd(0, n - 1);
    run(n, root, es, brute(n, root, es));
  }
  rep(it, 0, iters(1000)) {
    int n = (int)rnd(1, 6), m = (int)rnd(0, 2 * n);
    free_root(n, random_edges(n, m, rnd(0, 1) ? -6 : 0, 10));
  }
  {
    // Free root at its limit (n + 1) B <= 4e18: n vertices with only
    // self-loops, so every vertex takes the super edge.
    int n = 1000;
    ll w = 3'000'000'000'000, B = 1 + n * w;
    DirectedMST g(n + 1);
    rep(v, 0, n) g.add_edge(v, v, w), g.add_edge(n, v, B);
    auto [cost, par] = g.solve(n);
    CHECK_EQ(cost, n * B);
    CHECK_EQ(count(all(par), n), n);
  }
  rep(it, 0, iters(200)) {
    int n = (int)rnd(2, 60), m = (int)rnd(n, 6 * n);
    auto es = random_edges(n, m, -1000, 1000);
    rep(v, 1, n) es.push_back({rnd(0, v - 1), v, rnd(0, 2000)});
    int root = 0;
    run(n, root, es, slow(n, root, es));
  }
  // One big cycle of cheap edges: the whole ring contracts.
  int n = 100'000;
  vector<E> es;
  ll lo = BIG;
  rep(v, 1, n) {
    es.push_back({v, v % (n - 1) + 1, 1});
    ll w = rnd(1000, 1'000'000);
    es.push_back({0, v, w}), lo = min(lo, w);
  }
  run(n, 0, es, lo + n - 2);
  // Weights k, k+1, k+2, k-1, k+3, k-2, ... into one vertex build
  // a skew heap whose right spine has m/2 nodes; a recursive merge
  // overflows the stack on the first pop.
  ll k = 100'000;
  es = {{0, 1, k}, {0, 1, k + 1}};
  rep(i, 1, k) es.push_back({0, 1, k + 1 + i}), es.push_back({0, 1, k - i});
  es.push_back({0, 1, 2 * k + 1}); // one high weight after the minimum
  run(2, 0, es, 1);
  // Random graph with a spanning path from 0, against O(n m).
  n = 2000;
  es = random_edges(n, 20 * n, 0, 1'000'000'000);
  rep(v, 1, n) es.push_back({v - 1, v, rnd(0, 1'000'000'000)});
  run(n, 0, es, slow(n, 0, es));
  n = 100'000;
  es = random_edges(n, 2 * n, 0, 1'000'000'000);
  rep(v, 1, n) es.push_back({rnd(0, v - 1), v, rnd(0, 1'000'000'000)});
  double t0 = elapsed();
  DirectedMST g(n);
  for (auto [a, b, w] : es) g.add_edge((int)a, (int)b, w);
  auto [cost, par] = g.solve(0);
  fprintf(stderr, "  n = 1e5, m = 3e5: %.2fs\n", elapsed() - t0);
  check_tree(n, 0, es, cost, par);
  puts("OK graph/directed-mst");
}
