#include "../testlib.hpp"
#include "../../content/ds/lazy-segtree.hpp"
#include "../../content/tree/hld.hpp"

const ll MOD = 998244353;
// The header's Example: maps y -> b y + c, then(f, g) applies f
// first, and ds/segtree nodes that fold both ways. Both is the
// Example's T, with e = Both(0) and op(a, b) = a + b.
struct M {
  ll b, c;
};
M then(M f, M g) { return {g.b * f.b % MOD, (g.b * f.c + g.c) % MOD}; }
struct Both {
  M f{1, 0}, r{1, 0};
  constexpr Both(ll = 0) {}
  Both(M m) : f(m), r(m) {}
  friend Both operator+(Both a, Both b) {
    Both s;
    s.f = then(a.f, b.f), s.r = then(b.r, a.r);
    return s;
  }
};

// The snippet fixes T = ll; inside a template whose parameter is
// named ll the same code also runs over Both.
template<class ll> struct Generic {
#include "../../content/ds/segtree.hpp"
};
using SegTree = Generic<long long>::SegTree;
using BothTree = Generic<Both>::SegTree;

// Shapes: 0 random, 1 path, 2 star, 3 long thin, 4 binary-ish;
// labels shuffled.
vector<vi> random_tree(int n, int shape) {
  vi perm(n);
  iota(all(perm), 0), shuffle(all(perm), test_rng);
  vector<vi> adj(n);
  rep(i, 1, n) {
    int p = shape == 1   ? i - 1
            : shape == 2 ? 0
            : shape == 3 ? (int)rnd(max(0, i - 3), i - 1)
            : shape == 4 ? (i - 1) / 2
                         : (int)rnd(0, i - 1);
    adj[perm[i]].push_back(perm[p]), adj[perm[p]].push_back(perm[i]);
  }
  for (auto& a : adj) shuffle(all(a), test_rng);
  return adj;
}

struct Brute {
  vi par, dep, size;
  Brute(const vector<vi>& adj, int root)
      : par(sz(adj), -1), dep(sz(adj)), size(sz(adj), 1) {
    vi q{root}, seen(sz(adj));
    seen[root] = 1;
    rep(i, 0, sz(q))
      for (int v : adj[q[i]])
        if (!seen[v]) {
          seen[v] = 1, par[v] = q[i], dep[v] = dep[q[i]] + 1;
          q.push_back(v);
        }
    for (int i = sz(q) - 1; i > 0; i--) size[par[q[i]]] += size[q[i]];
  }
  int lca(int u, int v) {
    while (dep[u] > dep[v]) u = par[u];
    while (dep[v] > dep[u]) v = par[v];
    while (u != v) u = par[u], v = par[v];
    return u;
  }
  bool is_anc(int a, int v) {
    while (v >= 0 && v != a) v = par[v];
    return v == a;
  }
  // Path u .. v; side[i] = -1 strictly on u's side, 0 lca, 1 v's.
  vi path(int u, int v, vi& side) {
    int a = lca(u, v);
    vi left, right;
    for (; u != a; u = par[u]) left.push_back(u);
    for (; v != a; v = par[v]) right.push_back(v);
    side.assign(sz(left), -1), side.push_back(0);
    left.push_back(a);
    for (int i = sz(right); i--;) left.push_back(right[i]), side.push_back(1);
    return left;
  }
};

void check_tree(const vector<vi>& adj, int root) {
  int n = sz(adj);
  HLD h(adj, root);
  Brute b(adj, root);
  vi at(n, -1); // vertex at each position
  rep(v, 0, n) {
    CHECK_EQ(h.par[v], b.par[v]);
    CHECK_EQ(h.depth[v], b.dep[v]);
    CHECK_EQ(h.size[v], b.size[v]);
    CHECK(0 <= h.pos[v] && h.pos[v] < n && at[h.pos[v]] < 0);
    at[h.pos[v]] = v;
    int best = 0;
    for (int c : adj[v])
      if (c != b.par[v]) best = max(best, b.size[c]);
    if (best == 0) CHECK_EQ(h.heavy[v], -1);
    else {
      CHECK_EQ(b.par[h.heavy[v]], v);
      CHECK_EQ(b.size[h.heavy[v]], best);
      CHECK_EQ(h.pos[h.heavy[v]], h.pos[v] + 1);
    }
    bool top = v == root || h.heavy[b.par[v]] != v;
    CHECK_EQ(h.head[v], top ? v : h.head[b.par[v]]);
    auto [l, r] = h.subtree(v);
    rep(x, 0, n) {
      bool in = l <= h.pos[x] && h.pos[x] < r;
      CHECK_EQ(in, b.is_anc(v, x));
    }
    int light = 0;
    for (int x = v; x != root; x = b.par[x]) light += h.head[x] == x;
    CHECK(n == 1 || (1 << light) <= n);
  }
  vector<M> maps(n);
  for (auto& m : maps) m = {rnd(0, MOD - 1), rnd(0, MOD - 1)};
  vector<Both> leaves(n);
  rep(x, 0, n) leaves[h.pos[x]] = Both(maps[x]);
  BothTree st(leaves);
  rep(u, 0, n)
    rep(v, 0, n) {
      CHECK_EQ(h.lca(u, v), b.lca(u, v));
      vi side, p = b.path(u, v, side);
      for (bool edges : {false, true}) {
        // Walk the ranges: u's side right to left, v's side left to
        // right; the result must be the path, LCA dropped for edges.
        vi want, got;
        vector<pii> rs = h.path(u, v, edges);
        rep(i, 0, sz(p))
          if (!edges || side[i]) want.push_back(p[i]);
        for (auto [l, r] : rs) {
          CHECK(l < r);
          int s = 0;
          rep(x, l, r) s = side[find(all(p), at[x]) - p.begin()] ?: s;
          if (s < 0)
            for (int x = r; x-- > l;) got.push_back(at[x]);
          else
            rep(x, l, r) got.push_back(at[x]);
        }
        CHECK(got == want);
      }
      // Notes: a u - v path crosses at most 2 log2 n light edges
      // (log2 n each side), so path() gives at most that plus 2.
      int light = 0, a = b.lca(u, v);
      for (int x : p) light += x != a && h.head[x] == x;
      CHECK(light <= 2 * __lg(n));
      CHECK(sz(h.path(u, v)) <= light + 2);
      // The header's Example from `int w` on (st holds Both nodes).
      int w = h.lca(u, v);
      M f{1, 0};
      for (auto [l, r] : h.path(u, w)) f = then(f, st.query(l, r).r);
      for (auto [l, r] : h.path(w, v, true)) f = then(f, st.query(l, r).f);
      M want{1, 0};
      for (int x : p) want = then(want, maps[x]);
      CHECK(f.b == want.b && f.c == want.c);
    }
}

// The header's Usage: ds/segtree for path and subtree sums with
// point updates, ds/lazy-segtree for path adds and minimums.
void check_usage(const vector<vi>& adj, int root, int ops) {
  int n = sz(adj);
  HLD h(adj, root);
  Brute br(adj, root);
  vector<ll> a = rnd_vec(n, -1e9, 1e9), b(n);
  rep(v, 0, n) b[h.pos[v]] = a[v];
  SegTree st(b);
  LazySegTree lz(b);
  auto on_path = [&](int u, int v) {
    vi side;
    return br.path(u, v, side);
  };
  rep(it, 0, ops) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    ll x = rnd(-1e9, 1e9);
    int type = (int)rnd(0, 4);
    if (type == 0) { // point update
      a[v] = x, st.set(h.pos[v], x), lz.set(h.pos[v], x);
    } else if (type == 1) { // path add (lazy) and path sum
      for (auto [l, r] : h.path(u, v)) lz.apply(l, r, x);
      for (int y : on_path(u, v)) a[y] += x, st.set(h.pos[y], a[y]);
    } else if (type == 2) { // path sum, path min, edge sum
      ll s = 0, m = LLONG_MAX, se = 0, ws = 0, wm = LLONG_MAX;
      for (auto [l, r] : h.path(u, v))
        s += st.query(l, r), m = min(m, lz.query(l, r));
      for (auto [l, r] : h.path(u, v, true)) se += st.query(l, r);
      for (int y : on_path(u, v)) ws += a[y], wm = min(wm, a[y]);
      CHECK_EQ(s, ws);
      CHECK_EQ(m, wm);
      CHECK_EQ(se, ws - a[br.lca(u, v)]);
    } else { // subtree sum and min
      auto [l, r] = h.subtree(v);
      ll ws = 0, wm = LLONG_MAX;
      rep(y, 0, n)
        if (br.is_anc(v, y)) ws += a[y], wm = min(wm, a[y]);
      CHECK_EQ(st.query(l, r), ws);
      CHECK_EQ(lz.query(l, r), wm);
    }
  }
}

int main() {
  { // Reviewer's case: walked left to right, path(3, 2) is 3 0 1 2.
    vector<vi> adj{{1, 2}, {0, 3, 4}, {0}, {1}, {1, 5}, {4}};
    CHECK(HLD(adj).path(3, 2) == (vector<pii>{{4, 5}, {0, 2}, {5, 6}}));
    check_tree(adj, 0);
  }
  rep(it, 0, iters(1200)) {
    int n = (int)rnd(1, it < 100 ? 4 : 20);
    auto adj = random_tree(n, it % 5);
    int root = (int)rnd(0, n - 1);
    check_tree(adj, root);
    check_usage(adj, root, 30);
  }
  // Large: random recursive tree and a 2e5 path (deep), 1e5 path
  // sums checked against prefix sums from the root.
  int n = 200'000;
  rep(shape, 0, 2) {
    auto adj = random_tree(n, shape);
    int root = (int)rnd(0, n - 1);
    HLD h(adj, root);
    vector<ll> a = rnd_vec(n, -1e9, 1e9), b(n), pre(n);
    rep(v, 0, n) b[h.pos[v]] = a[v];
    SegTree st(b);
    vi q{root}; // prefix sums from the root, in BFS order
    pre[root] = a[root];
    rep(i, 0, sz(q))
      for (int v : adj[q[i]])
        if (v != h.par[q[i]]) pre[v] = pre[q[i]] + a[v], q.push_back(v);
    rep(it, 0, 100'000) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      int w = h.lca(u, v);
      CHECK(h.depth[w] <= min(h.depth[u], h.depth[v]));
      ll s = 0;
      for (auto [l, r] : h.path(u, v)) s += st.query(l, r);
      CHECK_EQ(s, pre[u] + pre[v] - 2 * pre[w] + a[w]);
    }
  }
  puts("OK tree/hld");
}
