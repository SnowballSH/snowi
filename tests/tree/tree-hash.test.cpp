#include "../testlib.hpp"
#include "../../content/tree/tree-hash.hpp"

const ll MOD = 998244353;

vector<vi> from_edges(int n, const vector<pii>& es) {
  vector<vi> adj(n);
  for (auto [u, v] : es) adj[u].push_back(v), adj[v].push_back(u);
  for (auto& a : adj) shuffle(all(a), test_rng);
  return adj;
}
vector<pii> edges_of(const vector<vi>& adj) {
  vector<pii> es;
  rep(u, 0, sz(adj))
    for (int v : adj[u])
      if (u < v) es.push_back({u, v});
  return es;
}
// Shapes: 0 random, 1 path, 2 star, 3 long thin, 4 binary-ish;
// labels shuffled.
vector<vi> random_tree(int n, int shape) {
  vi perm(n);
  iota(all(perm), 0), shuffle(all(perm), test_rng);
  vector<pii> es;
  rep(i, 1, n) {
    int p = shape == 1   ? i - 1
            : shape == 2 ? 0
            : shape == 3 ? (int)rnd(max(0, i - 3), i - 1)
            : shape == 4 ? (i - 1) / 2
                         : (int)rnd(0, i - 1);
    es.push_back({perm[i], perm[p]});
  }
  return from_edges(n, es);
}
vector<vi> relabel(const vector<vi>& adj, const vi& perm) {
  vector<pii> es;
  for (auto [u, v] : edges_of(adj)) es.push_back({perm[u], perm[v]});
  return from_edges(sz(adj), es);
}
vi random_perm(int n) {
  vi p(n);
  iota(all(p), 0), shuffle(all(p), test_rng);
  return p;
}

// Brute force over all relabellings (n <= 8): the smallest edge
// bitmask, with the root sent to 0 when root >= 0. Equal iff
// isomorphic. Also counts the relabellings that keep the edges
// (and the root): the automorphisms.
pair<ll, ll> canon(const vector<vi>& adj, int root) {
  int n = sz(adj);
  vector<pii> es = edges_of(adj);
  vi p(n);
  iota(all(p), 0);
  ll best = LLONG_MAX, self = 0, mine = 0;
  for (auto [u, v] : es) mine |= 1LL << (min(u, v) * 8 + max(u, v));
  do {
    ll m = 0;
    for (auto [u, v] : es) {
      int a = min(p[u], p[v]), b = max(p[u], p[v]);
      m |= 1LL << (a * 8 + b);
    }
    if (root < 0 || p[root] == 0) best = min(best, m);
    if (root < 0 || p[root] == root) self += m == mine;
  } while (next_permutation(all(p)));
  return {best, self};
}
// The subtree of v in adj rooted at r, relabelled 0..k-1, v -> 0.
vector<vi> subtree(const vector<vi>& adj, int r, int v) {
  int n = sz(adj);
  vi par(n, -1), q{r};
  par[r] = r;
  rep(i, 0, sz(q))
    for (int x : adj[q[i]])
      if (par[x] < 0) par[x] = q[i], q.push_back(x);
  vi idx(n, -1), st{v};
  vector<pii> es;
  idx[v] = 0;
  int k = 1;
  while (sz(st)) {
    int x = st.back();
    st.pop_back();
    for (int y : adj[x])
      if (y != par[x])
        idx[y] = k++, es.push_back({idx[x], idx[y]}), st.push_back(y);
  }
  return from_edges(k, es);
}
vi brute_centers(const vector<vi>& adj) {
  int n = sz(adj), best = INT_MAX;
  vi ecc(n), res;
  rep(s, 0, n) {
    vi d(n, -1), q{s};
    d[s] = 0;
    rep(i, 0, sz(q))
      for (int x : adj[q[i]])
        if (d[x] < 0) d[x] = d[q[i]] + 1, q.push_back(x);
    ecc[s] = *max_element(all(d)), best = min(best, ecc[s]);
  }
  rep(v, 0, n)
    if (ecc[v] == best) res.push_back(v);
  return res;
}

// One shared object and the canonical form behind every id seen:
// ids must be a bijection with isomorphism classes.
TreeHash th;
map<ll, int> rooted_id, unrooted_id;
map<int, ll> rooted_canon, unrooted_canon;
void consistent(map<ll, int>& c2i, map<int, ll>& i2c, ll c, int id) {
  CHECK(id >= 1 && id <= sz(th.ids));
  CHECK_EQ(c2i.emplace(c, id).first->second, id);
  CHECK_EQ(i2c.emplace(id, c).first->second, c);
}

// The header's Example: automorphisms fixing r, then unrooted.
ll automorphisms(const vector<vi>& adj, int r) {
  vi id = th.rooted(adj, r);
  map<pii, int> cnt;
  ll aut = 1;
  rep(v, 0, sz(adj))
    if (v != r) aut = aut * ++cnt[{th.par[v], id[v]}] % MOD;
  return aut;
}

void check_small(const vector<vi>& adj) {
  int n = sz(adj);
  vi cs = th.centers(adj);
  sort(all(cs));
  CHECK(cs == brute_centers(adj));
  auto [uc, uaut] = canon(adj, -1);
  consistent(unrooted_id, unrooted_canon, uc, th.unrooted(adj));
  ll aut = automorphisms(adj, cs[0]);
  if (sz(cs) == 2 &&
      th.rooted(adj, cs[0])[cs[0]] == th.rooted(adj, cs[1])[cs[1]])
    aut *= 2;
  CHECK_EQ(aut, uaut);
  int r = (int)rnd(0, n - 1);
  vi id = th.rooted(adj, r);
  vi par(n, -2), q{r}; // BFS parents from r
  par[r] = -1;
  rep(i, 0, sz(q))
    for (int x : adj[q[i]])
      if (par[x] == -2) par[x] = q[i], q.push_back(x);
  CHECK(th.par == par);
  auto [rc, raut] = canon(adj, r);
  consistent(rooted_id, rooted_canon, rc, id[r]);
  CHECK_EQ(automorphisms(adj, r), raut);
  rep(v, 0, n) {
    auto sub = subtree(adj, r, v);
    consistent(rooted_id, rooted_canon, canon(sub, 0).first, id[v]);
  }
}

vector<vi> from_prufer(const vi& code, int n) {
  vi deg(n, 1);
  for (int x : code) deg[x]++;
  vector<pii> es;
  for (int x : code)
    rep(v, 0, n)
      if (deg[v] == 1) {
        es.push_back({v, x}), deg[v]--, deg[x]--;
        break;
      }
  int a = -1;
  rep(v, 0, n)
    if (deg[v] == 1) {
      if (a < 0) a = v;
      else es.push_back({a, v});
    }
  return from_edges(n, es);
}

int main() {
  elapsed();
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, it < 500 ? 5 : 7);
    auto adj = random_tree(n, it % 5);
    check_small(adj);
    check_small(relabel(adj, random_perm(n)));
  }
  fprintf(stderr, "  small random trees: %.2fs\n", elapsed());
  // All labelled trees: non-isomorphic trees on n vertices are
  // 1 1 1 2 3 6 11 (OEIS A000055), rooted 1 1 2 4 9 20 48.
  vi unrooted_want{1, 1, 1, 2, 3, 6, 11}, rooted_want{1, 1, 2, 4, 9, 20, 48};
  rep(n, 1, 8) {
    set<int> un, ro;
    vi code(max(0, n - 2));
    while (true) {
      auto adj = from_prufer(code, n);
      un.insert(th.unrooted(adj));
      rep(r, 0, n) ro.insert(th.rooted(adj, r)[r]);
      int i = 0;
      while (i < sz(code) && ++code[i] == n) code[i++] = 0;
      if (i == sz(code)) break;
    }
    CHECK_EQ(sz(un), unrooted_want[n - 1]);
    CHECK_EQ(sz(ro), rooted_want[n - 1]);
  }
  fprintf(stderr, "  all labelled trees: %.2fs\n", elapsed());
  // Large: n = 1e5 random, path (two centers) and complete binary
  // tree, each against a relabelled copy; then a tree with fewer
  // leaves must differ.
  int n = 100'000;
  double t0 = elapsed();
  for (int shape : {0, 1, 4}) {
    auto a = random_tree(n, shape);
    vi perm = random_perm(n);
    auto b = relabel(a, perm);
    TreeHash h;
    int r = (int)rnd(0, n - 1);
    vi ia = h.rooted(a, r), ib = h.rooted(b, perm[r]);
    rep(v, 0, n) CHECK_EQ(ia[v], ib[perm[v]]);
    CHECK_EQ(h.unrooted(a), h.unrooted(b));
    vi ca = h.centers(a), cb = h.centers(b);
    for (int& c : ca) c = perm[c];
    sort(all(ca)), sort(all(cb));
    CHECK(ca == cb);
    if (shape == 1) CHECK_EQ(sz(ca), 2);
    if (shape == 0) { // move a leaf of x (degree >= 4) to a leaf y
      int x = -1, leaf = -1, y = -1;
      rep(v, 0, n)
        if (sz(a[v]) >= 4)
          for (int c : a[v])
            if (sz(a[c]) == 1) x = v, leaf = c;
      CHECK(x >= 0);
      rep(v, 0, n)
        if (sz(a[v]) == 1 && v != leaf) y = v;
      vector<pii> es;
      for (auto [u, v] : edges_of(a))
        if (minmax(u, v) != minmax(x, leaf)) es.push_back({u, v});
      es.push_back({y, leaf});
      CHECK(h.unrooted(a) != h.unrooted(from_edges(n, es)));
    }
  }
  fprintf(stderr, "  three trees of 1e5, each twice: %.2fs\n", elapsed() - t0);
  puts("OK tree/tree-hash");
}
