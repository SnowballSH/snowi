#include "../testlib.hpp"

// (sum, max) at once; SM(0) is the identity, so the snippet's
// `e = 0` and `op = a + b` work unchanged for this type.
struct SM {
  ll s = 0, m = LLONG_MIN;
  constexpr SM(int = 0) {}
  constexpr SM(ll s_, ll m_) : s(s_), m(m_) {}
  friend SM operator+(SM a, SM b) { return {a.s + b.s, max(a.m, b.m)}; }
  bool operator==(const SM&) const = default;
};
// The header's Example: T = (weight, edge id), op = max of pairs,
// identity {LLONG_MIN, -1}.
struct ME {
  ll w = LLONG_MIN;
  int id = -1;
  constexpr ME(int = 0) {}
  constexpr ME(ll w_, int id_) : w(w_), id(id_) {}
  friend ME operator+(ME a, ME b) {
    return pair(a.w, a.id) < pair(b.w, b.id) ? b : a;
  }
};

// The snippet fixes T = ll; inside a template whose parameter is
// named ll the same code runs over other types.
template<class ll> struct Generic {
#include "../../content/tree/link-cut.hpp"
};
using LinkCut = Generic<long long>::LinkCut;
using LinkCutSM = Generic<SM>::LinkCut;
using LinkCutME = Generic<ME>::LinkCut;

// An explicit forest: edge sets, one marked root per tree.
struct Brute {
  int n;
  vector<set<int>> adj;
  vi is_root;
  vector<ll> val;
  Brute(vector<ll> a) : n(sz(a)), adj(n), is_root(n, 1), val(a) {}
  vi comp(int s) { // BFS order of s's tree
    vi q{s}, seen(n);
    seen[s] = 1;
    rep(i, 0, sz(q))
      for (int v : adj[q[i]])
        if (!seen[v]) seen[v] = 1, q.push_back(v);
    return q;
  }
  int root(int v) {
    for (int x : comp(v))
      if (is_root[x]) return x;
    abort();
  }
  vi path(int u, int v) { // vertices from u to v
    vi par(n, -1), q{u};
    par[u] = u;
    rep(i, 0, sz(q))
      for (int x : adj[q[i]])
        if (par[x] < 0) par[x] = q[i], q.push_back(x);
    vi p{v};
    while (p.back() != u) p.push_back(par[p.back()]);
    reverse(all(p));
    return p;
  }
  void make_root(int v) { is_root[root(v)] = 0, is_root[v] = 1; }
  void link(int u, int v) {
    is_root[root(u)] = 0;
    adj[u].insert(v), adj[v].insert(u);
  }
  void cut(int u) {
    int r = root(u);
    if (r == u) return;
    int p = path(u, r)[1];
    adj[u].erase(p), adj[p].erase(u), is_root[u] = 1;
  }
  int lca(int u, int v) {
    vi a = path(u, root(u)), b = path(v, root(v));
    for (int x : a)
      if (count(all(b), x)) return x;
    abort();
  }
};

// Random operations on n vertices against the brute force, for the
// default (sum) and the (sum, max) instantiation.
void stress(int n, int ops) {
  vector<ll> a = rnd_vec(n, -1e9, 1e9);
  vector<SM> b(n);
  rep(i, 0, n) b[i] = SM(a[i], a[i]);
  LinkCut lc(a);
  LinkCutSM ls(b);
  Brute br(a);
  auto fold = [&](const vi& p) {
    SM r;
    for (int x : p) r = r + SM(br.val[x], br.val[x]);
    return r;
  };
  rep(it, 0, ops) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    bool same = br.root(u) == br.root(v);
    int type = (int)rnd(0, 9);
    if (type <= 1) { // link, or cut an edge when already connected
      if (!same) {
        lc.link(u, v), ls.link(u, v), br.link(u, v);
      } else if (u != v) {
        vi p = br.path(u, v);
        int k = (int)rnd(0, sz(p) - 2), x = p[k], y = p[k + 1];
        lc.cut(x, y), ls.cut(x, y), br.make_root(y), br.cut(x);
      }
    } else if (type == 2) { // cut u from its parent
      lc.cut(u), ls.cut(u), br.cut(u);
    } else if (type == 3) {
      lc.make_root(u), ls.make_root(u), br.make_root(u);
    } else if (type == 4) {
      ll x = rnd(-1e9, 1e9);
      lc.set(u, x), ls.set(u, SM(x, x)), br.val[u] = x;
    } else if (type == 5) {
      CHECK_EQ(lc.connected(u, v), same);
      CHECK_EQ(ls.connected(u, v), same);
      CHECK_EQ(lc.find_root(u), br.root(u));
      CHECK_EQ(ls.find_root(v), br.root(v));
    } else if (type == 6 && same) {
      CHECK_EQ(lc.lca(u, v), br.lca(u, v));
      CHECK_EQ(ls.lca(u, v), br.lca(u, v));
    } else if (type == 7 && same) { // path fold, root becomes u
      SM want = fold(br.path(u, v));
      CHECK_EQ(lc.query(u, v), want.s);
      CHECK(ls.query(u, v) == want);
      br.make_root(u);
    } else if (type == 8) { // root path, root unchanged
      SM want = fold(br.path(br.root(u), u));
      lc.access(u), ls.access(u);
      CHECK_EQ(lc.t[u].agg, want.s);
      CHECK(ls.t[u].agg == want);
    } else if (type == 9) { // values survive splaying
      CHECK_EQ(lc.t[u].val, br.val[u]);
      CHECK(ls.t[u].val == SM(br.val[u], br.val[u]));
    }
  }
}

// The header's Example against Kruskal on the edges seen so far.
void check_msf(int n, vi a, vi b, vector<ll> w) {
  int m = sz(w);
  LinkCutME lc(vector<LinkCutME::T>(n + m, LinkCutME::e));
  ll total = 0;
  auto add = [&](int i) {
    int u = a[i], v = b[i];
    if (lc.connected(u, v)) {
      auto [mw, j] = lc.query(u, v); // heaviest edge on path
      if (mw <= w[i]) return;
      lc.cut(n + j, a[j]), lc.cut(n + j, b[j]), total -= mw;
    }
    lc.set(n + i, {w[i], i});
    lc.link(n + i, u), lc.link(n + i, v), total += w[i];
  };
  rep(i, 0, m) {
    add(i);
    vi ord(i + 1), dsu(n);
    iota(all(ord), 0), iota(all(dsu), 0);
    sort(all(ord), [&](int x, int y) { return w[x] < w[y]; });
    auto find = [&](int x) {
      while (dsu[x] != x) x = dsu[x] = dsu[dsu[x]];
      return x;
    };
    ll want = 0;
    for (int j : ord)
      if (find(a[j]) != find(b[j])) want += w[j], dsu[find(a[j])] = find(b[j]);
    CHECK_EQ(total, want);
    rep(t, 0, 3) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      CHECK_EQ(lc.connected(u, v), find(u) == find(v));
    }
  }
}

int main() {
  { // n = 1: a lone vertex
    LinkCut lc(vector<ll>{7});
    CHECK_EQ(lc.query(0, 0), 7LL);
    CHECK_EQ(lc.lca(0, 0), 0);
    lc.cut(0), lc.set(0, -3);
    CHECK_EQ(lc.find_root(0), 0);
    CHECK_EQ(lc.query(0, 0), -3LL);
  }
  rep(it, 0, iters(1500)) stress((int)rnd(1, it < 500 ? 5 : 16), 60);
  // Negative weights: the identity must lie below every edge.
  check_msf(3, {0, 1, 2}, {1, 2, 0}, {-5, -3, -4});
  rep(it, 0, iters(200)) {
    int n = (int)rnd(1, 8), m = 25;
    ll maxw = rnd(1, 20);
    vi a(m), b(m);
    vector<ll> w(m);
    rep(i, 0, m) {
      a[i] = (int)rnd(0, n - 1), b[i] = (int)rnd(0, n - 1);
      w[i] = rnd(-maxw, maxw);
    }
    check_msf(n, a, b, w);
  }

  // Large 1: a path 0 - 1 - ... - n-1 (deep), values in a Fenwick
  // tree; random everts, cuts and relinks of path edges; path sums
  // are interval sums and the lca with root r is the clamp of r.
  int n = 100'000;
  double t0 = elapsed();
  {
    vector<ll> a = rnd_vec(n, -1e9, 1e9);
    LinkCut lc(a);
    rep(i, 1, n) lc.link(i - 1, i);
    // In order along the path: quadratic if splay rotates wrongly.
    rep(i, 0, n) CHECK_EQ(lc.find_root(i), n - 1);
    vector<ll> fw(n + 1); // Fenwick over a
    auto add = [&](int i, ll x) {
      for (i++; i <= n; i += i & -i) fw[i] += x;
    };
    auto pre = [&](int i) { // sum of a[0..i)
      ll s = 0;
      for (; i > 0; i -= i & -i) s += fw[i];
      return s;
    };
    rep(i, 0, n) add(i, a[i]);
    int root = n - 1;
    rep(it, 0, 100'000) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      int type = (int)rnd(0, 4);
      if (type == 0) {
        ll x = rnd(-1e9, 1e9);
        lc.set(u, x), add(u, x - a[u]), a[u] = x;
      } else if (type == 1) {
        lc.make_root(u), root = u;
      } else if (type == 2 && u + 1 < n) { // remove and restore
        lc.cut(u, u + 1), root = u + 1;
        CHECK(!lc.connected(u, u + 1));
        lc.link(u, u + 1);
      } else if (type == 3) {
        int lo = min(u, v), hi = max(u, v);
        CHECK_EQ(lc.lca(u, v), clamp(root, lo, hi));
      } else {
        CHECK_EQ(lc.query(u, v), pre(max(u, v) + 1) - pre(min(u, v)));
        root = u;
      }
    }
  }
  // Large 2: a random tree rooted at 0 whose subtrees move under
  // new parents (fixed root); lca and root path sums by climbing.
  {
    vector<ll> a = rnd_vec(n, -1e9, 1e9);
    LinkCut lc(a);
    vi par(n, -1), mark(n, -1);
    rep(i, 1, n) par[i] = (int)rnd(0, i - 1), lc.link(i, par[i]);
    rep(it, 0, 50'000) {
      int u = (int)rnd(1, n - 1), v = (int)rnd(0, n - 1);
      for (int y = u; y >= 0; y = par[y]) mark[y] = it;
      if (it % 2 == 0) { // move u's subtree under v
        int x = v;
        while (x >= 0 && x != u) x = par[x];
        if (x == u) continue; // v is in u's subtree
        lc.cut(u), lc.link(u, v), par[u] = v;
        CHECK_EQ(lc.find_root(u), 0);
      } else {
        int x = v;
        while (mark[x] != it) x = par[x];
        CHECK_EQ(lc.lca(u, v), x);
        ll s = 0;
        for (int y = v; y >= 0; y = par[y]) s += a[y];
        lc.access(v);
        CHECK_EQ(lc.t[v].agg, s);
      }
    }
  }
  fprintf(stderr, "  two large runs, n = 1e5: %.2fs\n", elapsed() - t0);
  puts("OK tree/link-cut");
}
