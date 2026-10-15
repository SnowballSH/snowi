#include "../testlib.hpp"
#include "../../content/tree/rerooting.hpp"

const ll MOD = 998244353;

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

vi bfs_dist(const vector<vi>& adj, int s) {
  vi d(sz(adj), -1), q{s};
  d[s] = 0;
  rep(i, 0, sz(q))
    for (int v : adj[q[i]])
      if (d[v] < 0) d[v] = d[q[i]] + 1, q.push_back(v);
  return d;
}

// The rooted DP from scratch with root r: sub[v] is the value of the
// subtree of v. Children are merged in adjacency order.
template<class T> vector<T> rooted_dp(const vector<vi>& adj, int r, T e,
                                      auto merge, auto lift, auto fin) {
  int n = sz(adj);
  vi par(n, -1), order{r};
  rep(i, 0, n)
    for (int v : adj[order[i]])
      if (v != par[order[i]]) par[v] = order[i], order.push_back(v);
  vector<T> sub(n);
  for (int i = n; i--;) {
    int v = order[i];
    T acc = e;
    for (int c : adj[v])
      if (c != par[v]) acc = merge(acc, lift(sub[c], c, v));
    sub[v] = fin(acc, v);
  }
  return sub;
}

// Compares every field of Reroot with DPs from every root.
template<class T>
void check_all(const vector<vi>& adj, T e, auto merge, auto lift, auto fin) {
  int n = sz(adj);
  Reroot<T> rr(adj, e, merge, lift, fin);
  vector<vector<T>> sub(n);
  rep(r, 0, n) sub[r] = rooted_dp(adj, r, e, merge, lift, fin);
  vi d0 = bfs_dist(adj, 0);
  CHECK_EQ(rr.par[0], -1);
  rep(v, 0, n) {
    CHECK(rr.dp[v] == sub[v][v]);
    if (v == 0) {
      CHECK(rr.dn[0] == e && rr.up[0] == e);
      continue;
    }
    int p = rr.par[v];
    CHECK(0 <= p && p < n && d0[p] == d0[v] - 1);
    CHECK(count(all(adj[v]), p) == 1);
    CHECK(rr.dn[v] == lift(sub[0][v], v, p));
    CHECK(rr.up[v] == lift(sub[v][p], p, v));
  }
}

typedef pair<ll, ll> P;

void check_tree(const vector<vi>& adj) {
  int n = sz(adj);
  // Usage: eccentricity.
  auto mx = [](int a, int b) { return max(a, b); };
  auto inc = [](int x, int, int) { return x + 1; };
  auto same = [](int acc, int) { return acc; };
  check_all(adj, 0, mx, inc, same);
  Reroot<int> ecc(adj, 0, mx, inc, same);
  rep(v, 0, n) {
    vi d = bfs_dist(adj, v);
    CHECK_EQ(ecc.dp[v], *max_element(all(d)));
  }
  // Example: sum of distances and subtree sizes.
  auto merge = [](P a, P b) {
    return P{a.first + b.first, a.second + b.second};
  };
  auto lift = [](P x, int, int) { return P{x.first, x.second + x.first}; };
  auto fin = [](P x, int) { return P{x.first + 1, x.second}; };
  check_all(adj, P{0, 0}, merge, lift, fin);
  Reroot<P> rr(adj, P{0, 0}, merge, lift, fin);
  vi d0 = bfs_dist(adj, 0);
  rep(v, 0, n) {
    vi d = bfs_dist(adj, v);
    CHECK(rr.dp[v] == P(n, accumulate(all(d), 0LL)));
    if (v) {
      int below = 0; // vertices whose path to 0 passes v
      rep(x, 0, n) below += d0[x] == d0[v] + d[x];
      CHECK_EQ(rr.dn[v].first, (ll)below);
    }
  }
  // Example: connected vertex sets containing v, by brute force
  // over all subsets for small n.
  Reroot<ll> cs(
      adj, 1LL, [](ll a, ll b) { return a * b % MOD; },
      [](ll x, int, int) { return (x + 1) % MOD; },
      [](ll acc, int) { return acc; });
  if (n <= 10) {
    vector<ll> want(n);
    rep(mask, 1, 1 << n) {
      int s = __builtin_ctz(mask), seen = 1 << s;
      vi st{s};
      while (sz(st)) {
        int u = st.back();
        st.pop_back();
        for (int v : adj[u])
          if ((mask >> v & 1) && !(seen >> v & 1))
            seen |= 1 << v, st.push_back(v);
      }
      if (seen == mask)
        rep(v, 0, n) want[v] += mask >> v & 1;
    }
    rep(v, 0, n) CHECK_EQ(cs.dp[v], want[v]);
  }
  // No inverse, depends on both endpoints and the direction: a
  // product with zeros (mod 7 and mod MOD).
  for (ll m : {7LL, MOD}) {
    vector<ll> A = rnd_vec(n, 0, 3), B = rnd_vec(n, 0, 3),
               C = rnd_vec(n, 0, m - 1), D = rnd_vec(n, 0, 2);
    check_all(
        adj, 1LL, [&](ll a, ll b) { return a * b % m; },
        [&](ll x, int u, int v) { return (x * A[u] + B[v] + D[u]) % m; },
        [&](ll acc, int v) { return (acc * D[v] + C[v]) % m; });
  }
  // Order-sensitive values (sorted multisets) catch a child merged
  // twice or skipped even when sums would cancel.
  typedef vector<ll> V;
  check_all(
      adj, V{},
      [](V a, V b) {
        a.insert(a.end(), all(b));
        sort(all(a));
        return a;
      },
      [](V x, int u, int v) {
        x.push_back(1000 * u + v);
        sort(all(x));
        return x;
      },
      [](V x, int v) {
        x.push_back(-v - 1);
        sort(all(x));
        return x;
      });
}

int main() {
  check_tree(vector<vi>(1));
  check_tree({{1}, {0}});
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, it < 300 ? 6 : 14);
    check_tree(random_tree(n, it % 5));
  }
  // Large: n = 2e5, a random tree, a path (deep) and a star (one
  // vertex with 2e5 children), checked with independent formulas.
  int n = 200'000;
  rep(shape, 0, 3) {
    auto adj = random_tree(n, shape);
    double t0 = elapsed();
    Reroot<P> rr(
        adj, P{0, 0},
        [](P a, P b) { return P{a.first + b.first, a.second + b.second}; },
        [](P x, int, int) { return P{x.first, x.second + x.first}; },
        [](P x, int) { return P{x.first + 1, x.second}; });
    Reroot<int> ecc(
        adj, 0, [](int a, int b) { return max(a, b); },
        [](int x, int, int) { return x + 1; },
        [](int acc, int) { return acc; });
    fprintf(stderr, "  shape %d, n = 2e5, two reroots: %.2fs\n", shape,
            elapsed() - t0);
    // Sum of distances: ans[c] = ans[p] + n - 2 size[c].
    vi d0 = bfs_dist(adj, 0), order(n), par(n, -1), size(n, 1);
    iota(all(order), 0);
    sort(all(order), [&](int a, int b) { return d0[a] < d0[b]; });
    for (int v : order)
      for (int u : adj[v])
        if (d0[u] == d0[v] + 1) par[u] = v;
    for (int i = n; --i;) size[par[order[i]]] += size[order[i]];
    vector<ll> ans(n);
    ans[0] = accumulate(all(d0), 0LL);
    for (int i = 1; i < n; i++) {
      int v = order[i];
      ans[v] = ans[par[v]] + n - 2LL * size[v];
    }
    // Eccentricity: max distance to the two ends of a diameter.
    int a = (int)(max_element(all(d0)) - d0.begin());
    vi da = bfs_dist(adj, a);
    int b = (int)(max_element(all(da)) - da.begin());
    vi db = bfs_dist(adj, b);
    rep(v, 0, n) {
      CHECK(rr.dp[v] == P(n, ans[v]));
      CHECK_EQ(ecc.dp[v], max(da[v], db[v]));
    }
  }
  puts("OK tree/rerooting");
}
