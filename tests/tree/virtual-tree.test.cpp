#include "../testlib.hpp"
#include "../../content/tree/lca-rmq.hpp"
#include "../../content/tree/virtual-tree.hpp"

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
  vi par, dep;
  vector<vi> anc; // all-pairs LCA by parent walking
  Brute(const vector<vi>& adj, int root) : par(sz(adj), -1), dep(sz(adj)) {
    int n = sz(adj);
    vi q{root}, seen(n);
    seen[root] = 1;
    rep(i, 0, sz(q))
      for (int v : adj[q[i]])
        if (!seen[v]) {
          seen[v] = 1, par[v] = q[i], dep[v] = dep[q[i]] + 1;
          q.push_back(v);
        }
    anc.assign(n, vi(n));
    rep(u, 0, n)
      rep(v, 0, n) {
        int x = u, y = v;
        while (dep[x] > dep[y]) x = par[x];
        while (dep[y] > dep[x]) y = par[y];
        while (x != y) x = par[x], y = par[y];
        anc[u][v] = x;
      }
  }
  int lca(int u, int v) { return anc[u][v]; }
};

// Closure under pairwise LCA, sorted by tin; parent = nearest proper
// ancestor inside the set.
void check_vt(Brute& b, const vi& tin, const vi& keys, const vi& vs,
              const vi& par) {
  set<int> K(all(keys)), S = K;
  for (int x : K)
    for (int y : K) S.insert(b.lca(x, y));
  CHECK_EQ(sz(vs), sz(S));
  CHECK(set<int>(all(vs)) == S);
  CHECK(sz(S) <= max(0, 2 * sz(set<int>(all(keys))) - 1));
  CHECK_EQ(sz(par), sz(vs));
  rep(i, 0, sz(vs)) {
    if (i) CHECK(tin[vs[i - 1]] < tin[vs[i]]);
    int x = b.par[vs[i]];
    while (x >= 0 && !S.count(x)) x = b.par[x];
    CHECK_EQ(par[i] < 0 ? -1 : vs[par[i]], x);
  }
}

void check_tree(const vector<vi>& adj, int root) {
  int n = sz(adj);
  LCA lc(adj, root);
  RmqLCA rl(adj, root);
  Brute b(adj, root);
  vi is_key(n);
  rep(q, 0, 8) {
    int k = q == 0 ? 0 : q == 1 ? 1 : (int)rnd(0, 2 * n);
    vi keys(k);
    for (int& v : keys) v = (int)rnd(0, n - 1);
    auto [vs, par] = virtual_tree(lc, keys);
    check_vt(b, lc.tin, keys, vs, par);
    auto [vs2, par2] = virtual_tree(rl, keys);
    check_vt(b, rl.tin, keys, vs2, par2);
    // The header's Example: sum of distances over pairs of keys.
    int m = sz(vs);
    k = 0;
    vector<ll> cnt(m); // keys in the subtree of vs[i]
    for (int v : keys) k += !is_key[v], is_key[v] = 1;
    rep(i, 0, m) cnt[i] = is_key[vs[i]];
    ll total = 0;
    for (int i = m - 1; i > 0; i--) {
      ll len = lc.depth[vs[i]] - lc.depth[vs[par[i]]];
      total += len * cnt[i] * (k - cnt[i]);
      cnt[par[i]] += cnt[i];
    }
    for (int v : keys) is_key[v] = 0;
    set<int> S(all(keys));
    ll want = 0;
    for (int x : S)
      for (int y : S)
        if (x < y) want += b.dep[x] + b.dep[y] - 2 * b.dep[b.lca(x, y)];
    CHECK_EQ(total, want);
    // Usage: children before parents.
    vi dp(m, 1);
    for (int i = sz(vs) - 1; i > 0; i--) dp[par[i]] += dp[i];
    if (m) CHECK_EQ(dp[0], m);
  }
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, it < 100 ? 4 : 30);
    auto adj = random_tree(n, it % 5);
    check_tree(adj, (int)rnd(0, n - 1));
  }
  // Large: queries with sum k = 1e5 on a 2e5 path rooted at one end
  // (nearest ancestor in S = previous index in S) and a random tree.
  int n = 200'000;
  vi perm(n), idx(n);
  iota(all(perm), 0), shuffle(all(perm), test_rng);
  rep(i, 0, n) idx[perm[i]] = i;
  vector<vi> path(n);
  rep(i, 1, n)
    path[perm[i]].push_back(perm[i - 1]), path[perm[i - 1]].push_back(perm[i]);
  LCA lp(path, perm[0]);
  for (int done = 0; done < 100'000;) {
    int k = (int)rnd(1, 2000);
    vi keys(k);
    for (int& v : keys) v = (int)rnd(0, n - 1);
    done += k;
    auto [vs, par] = virtual_tree(lp, keys);
    vi pos;
    for (int v : keys) pos.push_back(idx[v]);
    sort(all(pos)), pos.erase(unique(all(pos)), pos.end());
    CHECK_EQ(sz(vs), sz(pos));
    rep(i, 0, sz(vs)) {
      CHECK_EQ(idx[vs[i]], pos[i]);
      CHECK_EQ(par[i], i - 1);
    }
  }
  auto adj = random_tree(n, 0);
  LCA lc(adj, 0);
  RmqLCA rl(adj, 0);
  for (int done = 0; done < 100'000;) {
    int k = (int)rnd(1, 2000);
    vi keys(k);
    for (int& v : keys) v = (int)rnd(0, n - 1);
    done += k;
    auto [vs, par] = virtual_tree(lc, keys);
    auto [vs2, par2] = virtual_tree(rl, keys);
    CHECK(vs == vs2 && par == par2);
    CHECK(sz(vs) <= 2 * k - 1);
    rep(i, 1, sz(vs)) {
      int p = vs[par[i]];
      CHECK(lc.tin[p] < lc.tin[vs[i]] && lc.lca(p, vs[i]) == p);
    }
  }
  puts("OK tree/virtual-tree");
}
