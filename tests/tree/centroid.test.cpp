#include "../testlib.hpp"
#include "../../content/tree/centroid.hpp"

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

vector<vi> all_dists(const vector<vi>& adj) {
  int n = sz(adj);
  vector<vi> d(n, vi(n, -1));
  rep(s, 0, n) {
    vi q{s};
    d[s][s] = 0;
    rep(i, 0, sz(q))
      for (int v : adj[q[i]])
        if (d[s][v] < 0) { d[s][v] = d[s][q[i]] + 1, q.push_back(v); }
  }
  return d;
}

// Component of s among vertices with ok[v], by BFS.
vi component(const vector<vi>& adj, int s, const vector<bool>& ok) {
  vi q{s};
  vector<bool> seen(sz(adj));
  seen[s] = true;
  rep(i, 0, sz(q))
    for (int v : adj[q[i]])
      if (ok[v] && !seen[v]) { seen[v] = true, q.push_back(v); }
  sort(all(q));
  return q;
}

void check_tree(const vector<vi>& adj) {
  int n = sz(adj);
  Centroid cd(adj);
  CHECK_EQ(sz(cd.par), n);
  auto d = all_dists(adj);
  int roots = 0, lg = __lg(n);
  vector<vi> piece(n); // piece[c]: v whose centroid chain has c
  rep(v, 0, n) {
    roots += cd.par[v] < 0;
    int p = cd.par[v];
    CHECK_EQ(cd.level[v], p < 0 ? 0 : cd.level[p] + 1);
    CHECK(cd.level[v] <= lg);
    for (int c = v; c >= 0; c = cd.par[c]) {
      piece[c].push_back(v);
      CHECK_EQ(cd.dist[cd.level[c]][v], d[c][v]);
    }
  }
  CHECK_EQ(roots, 1);
  rep(c, 0, n) {
    // The piece is c's component once higher centroids are gone.
    vector<bool> ok(n);
    rep(v, 0, n) ok[v] = cd.level[v] >= cd.level[c];
    sort(all(piece[c]));
    CHECK(component(adj, c, ok) == piece[c]);
    // c is a centroid: removing it leaves parts of size <= half.
    ok[c] = false;
    for (int v : adj[c])
      if (ok[v]) CHECK(2 * sz(component(adj, v, ok)) <= sz(piece[c]));
  }
  // Usage: nearest marked vertex.
  int INF = INT_MAX / 2;
  vi best(n, INT_MAX / 2);
  auto mark = [&](int v) {
    for (int c = v; c >= 0; c = cd.par[c])
      best[c] = min(best[c], cd.dist[cd.level[c]][v]);
  };
  auto nearest = [&](int v) { // INT_MAX / 2 if none
    int r = INT_MAX / 2;
    for (int c = v; c >= 0; c = cd.par[c])
      r = min(r, best[c] + cd.dist[cd.level[c]][v]);
    return r;
  };
  vi marked;
  rep(op, 0, 2 * n) {
    int v = (int)rnd(0, n - 1);
    if (rnd(0, 2) == 0) mark(v), marked.push_back(v);
    int want = INF;
    for (int x : marked) want = min(want, d[v][x]);
    CHECK_EQ(nearest(v), want);
  }
  // Tweak: sum of distances to a changing marked set.
  {
    vector<ll> cnt(n), sum(n), sub(n);
    vi on(n);
    rep(op, 0, 2 * n) {
      int v = (int)rnd(0, n - 1);
      if (rnd(0, 1) && !on[v]) {
        on[v] = 1;
        for (int c = v; c >= 0; c = cd.par[c]) {
          cnt[c]++, sum[c] += cd.dist[cd.level[c]][v];
          if (cd.par[c] >= 0) sub[c] += cd.dist[cd.level[c] - 1][v];
        }
      }
      ll got = 0, want = 0;
      for (int c = v, ch = -1; c >= 0; ch = c, c = cd.par[c]) {
        ll dv = cd.dist[cd.level[c]][v];
        got += sum[c] + cnt[c] * dv;
        if (ch >= 0) got -= sub[ch] + cnt[ch] * dv;
      }
      rep(x, 0, n)
        if (on[x]) want += d[v][x];
      CHECK_EQ(got, want);
    }
  }
  // Example: pairs at distance exactly K, for every K.
  rep(K, 1, n + 1) {
    vi size(n); // piece sizes
    rep(v, 0, n)
      for (int c = v; c >= 0; c = cd.par[c]) size[c]++;
    vector<vi> cnt(n), sub(n);
    rep(c, 0, n) cnt[c].resize(size[c]);
    rep(c, 0, n) sub[c].resize(size[c] + 1);
    ll pairs = 0;
    rep(v, 0, n) {
      int ch = -1; // the previous centroid on v's way up
      for (int c = v; c >= 0; ch = c, c = cd.par[c]) {
        int dv = cd.dist[cd.level[c]][v], e = K - dv;
        if (0 <= e && e < size[c]) pairs += cnt[c][e];
        if (ch >= 0 && 0 <= e && e <= size[ch]) pairs -= sub[ch][e];
        cnt[c][dv]++;
        if (ch >= 0) sub[ch][dv]++;
      }
    }
    ll want = 0;
    rep(u, 0, n)
      rep(v, 0, u) want += d[u][v] == K;
    CHECK_EQ(pairs, want);
  }
}

int main() {
  rep(it, 0, iters(1200)) {
    int n = (int)rnd(1, it < 100 ? 4 : 24);
    check_tree(random_tree(n, it % 5));
  }
  // Large: a 2e5 path (deep), nearest marked vertex against a set of
  // marked positions; a random tree and a star for the level bound.
  int n = 200'000;
  vi perm(n);
  iota(all(perm), 0), shuffle(all(perm), test_rng);
  vector<vi> path(n);
  rep(i, 1, n)
    path[perm[i]].push_back(perm[i - 1]), path[perm[i - 1]].push_back(perm[i]);
  Centroid cd(path);
  vi best(n, INT_MAX / 2);
  set<int> marked;
  rep(it, 0, 100'000) {
    int i = (int)rnd(0, n - 1), v = perm[i];
    if (it % 2) {
      marked.insert(i);
      for (int c = v; c >= 0; c = cd.par[c])
        best[c] = min(best[c], cd.dist[cd.level[c]][v]);
    } else {
      int r = INT_MAX / 2, want = INT_MAX / 2;
      for (int c = v; c >= 0; c = cd.par[c])
        r = min(r, best[c] + cd.dist[cd.level[c]][v]);
      auto it2 = marked.lower_bound(i);
      if (it2 != marked.end()) want = *it2 - i;
      if (it2 != marked.begin()) want = min(want, i - *prev(it2));
      CHECK_EQ(r, want);
    }
  }
  for (int shape : {0, 2, 4}) {
    auto adj = random_tree(n, shape);
    Centroid c2(adj);
    vi size(n);
    rep(v, 0, n)
      for (int c = v; c >= 0; c = c2.par[c]) size[c]++;
    rep(c, 0, n) CHECK(size[c] <= n >> c2.level[c]);
  }
  puts("OK tree/centroid");
}
