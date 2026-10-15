#include "../testlib.hpp"
#include "../../content/flow/gomory-hu.hpp"

typedef array<ll, 3> E3;
typedef vector<vector<ll>> Mat;

ll cut_weight(const vector<E3>& edges, const vector<bool>& in) {
  ll c = 0;
  for (auto [u, v, w] : edges)
    if (in[u] != in[v]) c += w;
  return c;
}

// All-pairs minimum cut by enumerating every vertex subset.
Mat brute(int n, const vector<E3>& edges) {
  Mat lam(n, vector<ll>(n, LLONG_MAX));
  rep(mask, 0, 1 << n) {
    vector<bool> in(n);
    rep(v, 0, n) in[v] = mask >> v & 1;
    ll c = cut_weight(edges, in);
    rep(u, 0, n)
      rep(v, 0, n)
        if (in[u] && !in[v]) lam[u][v] = min(lam[u][v], c);
  }
  return lam;
}

// Edmonds-Karp on a capacity matrix (BFS shortest augmenting paths).
ll edmonds_karp(int n, const vector<E3>& edges, int s, int t) {
  Mat c(n, vector<ll>(n));
  for (auto [u, v, w] : edges)
    if (u != v) c[u][v] += w, c[v][u] += w;
  ll total = 0;
  while (true) {
    vi par(n, -1);
    par[s] = s;
    vi q{s};
    rep(i, 0, sz(q))
      rep(v, 0, n)
        if (par[v] < 0 && c[q[i]][v] > 0) par[v] = q[i], q.push_back(v);
    if (par[t] < 0) return total;
    ll f = LLONG_MAX;
    for (int v = t; v != s; v = par[v]) f = min(f, c[par[v]][v]);
    for (int v = t; v != s; v = par[v]) c[par[v]][v] -= f, c[v][par[v]] += f;
    total += f;
  }
}

// Checks that tree is a spanning tree whose every edge splits the
// vertices into a cut of exactly its weight. Returns the lightest
// edge on each tree path (LLONG_MAX on the diagonal).
Mat certify(int n, const vector<E3>& edges, const vector<E3>& tree) {
  CHECK_EQ(sz(tree), n ? n - 1 : 0);
  vector<vector<pii>> adj(n); // (neighbour, tree edge id)
  rep(i, 0, sz(tree)) {
    auto [u, v, w] = tree[i];
    CHECK(0 <= u && u < n && 0 <= v && v < n && u != v && w >= 0);
    adj[u].push_back({(int)v, i}), adj[v].push_back({(int)u, i});
  }
  Mat low(n, vector<ll>(n, -1));
  rep(s, 0, n) {
    low[s][s] = LLONG_MAX;
    vi st{s};
    while (sz(st)) {
      int u = st.back();
      st.pop_back();
      for (auto [v, i] : adj[u])
        if (low[s][v] < 0)
          low[s][v] = min(low[s][u], tree[i][2]), st.push_back(v);
    }
    rep(v, 0, n) CHECK(low[s][v] >= 0); // connected
  }
  rep(i, 0, sz(tree)) {
    vector<bool> in(n);
    vi st{(int)tree[i][0]};
    in[st[0]] = true;
    while (sz(st)) {
      int u = st.back();
      st.pop_back();
      for (auto [v, j] : adj[u])
        if (j != i && !in[v]) in[v] = true, st.push_back(v);
    }
    CHECK(!in[tree[i][1]]);
    CHECK_EQ(cut_weight(edges, in), tree[i][2]);
  }
  return low;
}

vector<E3> random_edges(int n, int m, ll maxc) {
  vector<E3> e;
  rep(i, 0, m) e.push_back({rnd(0, n - 1), rnd(0, n - 1), rnd(0, maxc)});
  return e;
}

// Heavy random clusters joined by light edges: deep, uneven trees.
vector<E3> clustered(int n, int m) {
  vi cl(n);
  rep(v, 0, n) cl[v] = (int)rnd(0, 4);
  vector<E3> e;
  rep(i, 0, m) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    e.push_back({u, v, cl[u] == cl[v] ? rnd(5, 20) : rnd(0, 2)});
  }
  return e;
}

int main() {
  CHECK_EQ(sz(gomory_hu(0, {})), 0);
  CHECK_EQ(sz(gomory_hu(1, {{0, 0, 5}})), 0);
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(1, 8), m = (int)rnd(0, 14);
    ll maxc = rnd(0, 3) ? rnd(1, 10) : (ll)1e15;
    auto e = random_edges(n, m, maxc);
    Mat low = certify(n, e, gomory_hu(n, e)), lam = brute(n, e);
    rep(u, 0, n)
      rep(v, 0, n)
        if (u != v) CHECK_EQ(low[u][v], lam[u][v]);
  }
  rep(it, 0, iters(40)) {
    int n = (int)rnd(10, 30);
    int m = (int)rnd(n, 4 * n);
    auto e = rnd(0, 1) ? random_edges(n, m, rnd(1, 20)) : clustered(n, m);
    Mat low = certify(n, e, gomory_hu(n, e));
    rep(k, 0, 20) {
      int s = (int)rnd(0, n - 1), t = (int)rnd(0, n - 2);
      t += t >= s;
      CHECK_EQ(low[s][t], edmonds_karp(n, e, s, t));
    }
  }
  // Large: certify the cut property, then compare pairs with
  // fresh max flows.
  elapsed();
  {
    int n = 800;
    auto e = random_edges(n, 8000, (ll)1e9);
    auto tree = gomory_hu(n, e);
    fprintf(stderr, "  n = 800, m = 8000: %.2fs\n", elapsed());
    Mat low = certify(n, e, tree);
    rep(k, 0, 30) {
      int s = (int)rnd(0, n - 1), t = (int)rnd(0, n - 2);
      t += t >= s;
      Dinic d(n);
      for (auto [u, v, c] : e) d.add_edge(u, v, c, c);
      CHECK_EQ(low[s][t], d.max_flow(s, t));
    }
  }
  puts("OK flow/gomory-hu");
}
