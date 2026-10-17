#include "../testlib.hpp"
#include "../../content/contest/gen.hpp"

struct Dsu {
  vi p;
  Dsu(int n) : p(n, -1) {}
  int find(int x) { return p[x] < 0 ? x : p[x] = find(p[x]); }
  bool join(int a, int b) {
    a = find(a), b = find(b);
    if (a == b) return false;
    p[a] = b;
    return true;
  }
};

vi degrees(int n, const vector<pii>& e) {
  vi deg(n);
  for (auto [u, v] : e) deg[u]++, deg[v]++;
  return deg;
}

// Endpoints in range, no loops, no unordered pair twice.
bool simple(int n, const vector<pii>& e) {
  set<pii> seen;
  for (auto [u, v] : e)
    if (u < 0 || v < 0 || u >= n || v >= n || u == v ||
        !seen.insert(minmax(u, v)).second)
      return false;
  return true;
}

int components(int n, const vector<pii>& e) {
  Dsu d(n);
  int c = n;
  for (auto [u, v] : e) c -= d.join(u, v);
  return c;
}

bool is_tree(int n, const vector<pii>& e) {
  return sz(e) == n - 1 && simple(n, e) && components(n, e) == 1;
}

bool acyclic(int n, const vector<pii>& e) {
  vector<vi> adj(n);
  vi indeg(n), order;
  for (auto [u, v] : e) adj[u].push_back(v), indeg[v]++;
  rep(v, 0, n)
    if (!indeg[v]) order.push_back(v);
  rep(i, 0, sz(order))
    for (int w : adj[order[i]])
      if (!--indeg[w]) order.push_back(w);
  return sz(order) == n;
}

// Pearson statistic of counts against `cells` equally likely cells;
// cells never seen count too.
template<class K> double chi2(const map<K, int>& cnt, int cells) {
  ll total = 0;
  for (auto& kv : cnt) total += kv.second;
  double want = (double)total / cells, s = want * (cells - sz(cnt));
  for (auto& kv : cnt) s += (kv.second - want) * (kv.second - want) / want;
  return s;
}
// A correct generator fails with probability below about 1e-6.
bool chi2_ok(double stat, int cells) {
  int dof = cells - 1;
  return stat < dof + 7 * sqrt(2.0 * dof) + 10;
}

vector<pii> normalized(vector<pii> e) {
  for (auto& [u, v] : e)
    if (u > v) swap(u, v);
  sort(all(e));
  return e;
}

const int SHAPES[] = {Gen::UNIFORM, Gen::PARENT, Gen::PATH, Gen::STAR,
                      Gen::CATERPILLAR};

void check_shape(int n, int shape, const vector<pii>& e) {
  CHECK(is_tree(n, e));
  vi deg = degrees(n, e);
  int max_deg = n == 1 ? 0 : *max_element(all(deg));
  if (shape == Gen::PATH) CHECK(max_deg <= 2);
  if (shape == Gen::STAR) CHECK(max_deg == n - 1);
  if (shape == Gen::CATERPILLAR) { // the non-leaves form a path
    vi inner(n);
    for (auto [u, v] : e)
      if (deg[u] > 1 && deg[v] > 1) inner[u]++, inner[v]++;
    CHECK(*max_element(all(inner)) <= 2);
  }
}

void test_numbers() {
  rep(it, 0, iters(300)) {
    Gen g(rnd(0, 1e18));
    ll lo = rnd(-5, 5), hi = lo + rnd(0, 6);
    set<ll> hit;
    rep(i, 0, 200) {
      ll x = g.rnd(lo, hi);
      CHECK(lo <= x && x <= hi);
      hit.insert(x);
    }
    CHECK(hit.count(lo) && hit.count(hi));
    ll big = g.rnd(-(ll)1e18, (ll)1e18);
    CHECK(-(ll)1e18 <= big && big <= (ll)1e18);

    int n = (int)rnd(0, 12);
    bool tight = rnd(0, 1);
    hi = lo + (tight ? n - 1 : rnd(n, 3 * n + 3));
    for (bool distinct : {false, true}) {
      auto a = g.array(n, lo, hi, distinct);
      CHECK_EQ(sz(a), n);
      for (ll x : a) CHECK(lo <= x && x <= hi);
      if (distinct) CHECK_EQ(sz(set<ll>(all(a))), n);
    }
    vi p = g.perm(n), id(n);
    iota(all(id), 0);
    sort(all(p));
    CHECK(p == id);

    string alphabet = rnd(0, 1) ? "ab" : "xyz";
    string s = g.str(n, alphabet);
    CHECK_EQ(sz(s), n);
    for (char c : s) CHECK(alphabet.find(c) != string::npos);
  }
  // Huge range, distinct.
  Gen g(rnd(0, 1e18));
  auto a = g.array(1000, -(ll)1e18, (ll)1e18, true);
  CHECK_EQ(sz(set<ll>(all(a))), 1000);
  for (ll x : a) CHECK(-(ll)1e18 <= x && x <= (ll)1e18);
  string s = g.str(1000);
  CHECK(count(all(s), 'a') > 400 && count(all(s), 'b') > 400);
  CHECK_EQ(count(all(s), 'a') + count(all(s), 'b'), 1000L);
  // Every letter of a longer alphabet, equally often.
  map<char, int> letters;
  for (char c : g.str(5000, "vwxyz")) letters[c]++;
  CHECK_EQ(sz(letters), 5);
  CHECK(chi2_ok(chi2(letters, 5), 5));
  // Orders: distinct arrays as ordered tuples, and permutations.
  map<vector<ll>, int> tuples, pairs;
  map<vi, int> perms;
  rep(it, 0, 24000) {
    tuples[g.array(2, 3, 6, true)]++, pairs[g.array(2, 0, 1)]++;
    perms[g.perm(4)]++;
  }
  CHECK(chi2_ok(chi2(tuples, 12), 12));
  CHECK(chi2_ok(chi2(pairs, 4), 4));
  CHECK(chi2_ok(chi2(perms, 24), 24));
}

// Ranges at the 64-bit limits; the sanitizer aborts on overflow.
void test_extremes() {
  const ll lo = LLONG_MIN, hi = LLONG_MAX;
  Gen g(rnd(0, 1e18));
  CHECK(g.rnd(lo, lo) == lo && g.rnd(hi, hi) == hi);
  using Case = tuple<int, ll, ll>;
  for (auto [n, a, b] :
       {Case(3, 0, hi), Case(1, lo, hi), Case(4, lo, hi), Case(2, lo, lo + 1),
        Case(5, lo, lo + 4), Case(2, hi - 1, hi), Case(3, hi - 2, hi)})
    for (bool distinct : {false, true}) {
      auto v = g.array(n, a, b, distinct);
      CHECK_EQ(sz(v), n);
      for (ll x : v) CHECK(a <= x && x <= b);
      if (distinct) CHECK_EQ(sz(set<ll>(all(v))), n);
    }
}

void test_seeds() {
  rep(it, 0, iters(50)) {
    ll seed = rnd(0, 1e18);
    Gen a(seed), b(seed);
    int n = (int)rnd(20, 40), m = (int)rnd(n - 1, 2 * n);
    CHECK(a.tree(n) == b.tree(n));
    CHECK(a.graph(n, m) == b.graph(n, m));
    CHECK(a.array(n, 0, 1e9, true) == b.array(n, 0, 1e9, true));
    CHECK_EQ(a.str(n), b.str(n));
    CHECK(Gen(seed).perm(n) != Gen(seed + 1).perm(n));
  }
}

void test_trees() {
  rep(it, 0, iters(3000)) {
    Gen g(rnd(0, 1e18));
    int n = it < 50 ? 1 + it % 3 : (int)rnd(1, 30);
    for (int shape : SHAPES) { // with m = n - 1, graph is its tree
      check_shape(n, shape, g.tree(n, shape));
      check_shape(n, shape, g.graph(n, n - 1, true, shape));
    }
  }
  // Uniform labelled trees: n^(n-2) of them, all equally likely.
  Gen g(rnd(0, 1e18));
  for (auto [n, samples] : {pii(3, 3000), pii(4, 16000), pii(5, 37500)}) {
    map<vector<pii>, int> cnt;
    rep(it, 0, samples) cnt[normalized(g.tree(n))]++;
    int cells = (int)pow(n, n - 2);
    CHECK_EQ(sz(cnt), cells);
    CHECK(chi2_ok(chi2(cnt, cells), cells));
  }
  // Relabelled, shuffled and flipped: star centres and path
  // ends are uniform, the centre is first in half of the edges,
  // and consecutive path edges rarely share a vertex.
  map<int, int> centre, ends;
  int first = 0, total = 0, adjacent = 0;
  rep(it, 0, 6000) {
    auto star = g.tree(6, Gen::STAR);
    vi deg = degrees(6, star);
    int c = int(max_element(all(deg)) - deg.begin());
    centre[c]++;
    for (auto [u, v] : star) first += u == c, total++;
    auto path = g.tree(6, Gen::PATH);
    deg = degrees(6, path);
    rep(v, 0, 6)
      if (deg[v] == 1) ends[v]++;
    rep(i, 1, 5) {
      auto [a, b] = path[i - 1];
      auto [x, y] = path[i];
      adjacent += a == x || a == y || b == x || b == y;
    }
  }
  CHECK(chi2_ok(chi2(centre, 6), 6));
  CHECK(chi2_ok(chi2(ends, 6), 6));
  CHECK(abs(first - total / 2) < total / 50);
  // 4 of 10 edge pairs are adjacent in a shuffled path of 5 edges.
  CHECK(abs(adjacent - 6000 * 4 * 4 / 10) < 6000 * 4 / 20);
  // Random-parent trees on 4 vertices: 2 of the 6 equally likely
  // parent arrays are stars, so 1/3 of them (uniform: 4/16).
  int stars = 0;
  rep(it, 0, 12000) {
    vi deg = degrees(4, g.tree(4, Gen::PARENT));
    stars += *max_element(all(deg)) == 3;
  }
  CHECK(abs(stars - 4000) < 300);
  // Caterpillars of 8 vertices: every spine (non-leaf) size from 1
  // (a star) to 6 (a path) occurs, and most are neither (~70%).
  map<int, int> spine;
  rep(it, 0, 3000) {
    vi deg = degrees(8, g.tree(8, Gen::CATERPILLAR));
    spine[(int)count_if(all(deg), [](int d) { return d > 1; })]++;
  }
  CHECK_EQ(sz(spine), 6);
  CHECK(spine[2] + spine[3] + spine[4] + spine[5] > 1500);
}

void test_graphs() {
  rep(it, 0, iters(3000)) {
    Gen g(rnd(0, 1e18));
    int n = it < 60 ? 1 + it % 4 : (int)rnd(1, 14);
    ll pairs = (ll)n * (n - 1) / 2;
    bool connected = rnd(0, 1);
    ll lo = connected ? n - 1 : 0;
    int m = (int)(it % 3 == 0 ? pairs : it % 3 == 1 ? lo : rnd(lo, pairs));
    auto e = g.graph(n, m, connected);
    CHECK_EQ(sz(e), m);
    CHECK(simple(n, e));
    if (connected) CHECK_EQ(components(n, e), 1);
    auto d = g.dag(n, m);
    CHECK_EQ(sz(d), m);
    CHECK(simple(n, d));
    CHECK(acyclic(n, d));
  }
  // m random pairs of n = 5 vertices: all C(10, 3) sets equally
  // likely, so the pair decoding is a bijection.
  Gen g(rnd(0, 1e18));
  map<vector<pii>, int> sets;
  rep(it, 0, 24000) sets[normalized(g.graph(5, 3, false))]++;
  CHECK(chi2_ok(chi2(sets, 120), 120));
  // With m = n - 1 a connected graph is its spanning tree, which is
  // uniform by default: all 16 labelled trees on 4 vertices.
  map<vector<pii>, int> spanning;
  rep(it, 0, 16000) spanning[normalized(g.graph(4, 3))]++;
  CHECK_EQ(sz(spanning), 16);
  CHECK(chi2_ok(chi2(spanning, 16), 16));
  // Transitive tournaments on 4 vertices: the topological order is a
  // random one of the 24.
  map<vi, int> orders;
  rep(it, 0, 24000) {
    vi pos(4);
    for (auto [u, v] : g.dag(4, 6)) pos[u]++;
    orders[pos]++;
  }
  CHECK(chi2_ok(chi2(orders, 24), 24));
  // Edge orientation of graphs is random too.
  int forward = 0;
  rep(it, 0, 2000)
    for (auto [u, v] : g.graph(6, 9)) forward += u < v;
  CHECK(abs(forward - 9000) < 600);
}

// Generation time only, printed; the sanitized build is ~5x slower.
void test_large() {
  Gen g(rnd(0, 1e18));
  int n = 200000, k = 700; // k: complete graph, 244650 edges
  double start = elapsed();
  vector<vector<pii>> trees;
  for (int shape : SHAPES) trees.push_back(g.tree(n, shape));
  double mid = elapsed();
  auto e = g.graph(n, 2 * n), d = g.dag(n, 2 * n);
  auto complete = g.graph(k, k * (k - 1) / 2);
  auto a = g.array(n, 0, n - 1, true);
  double stop = elapsed();
  fprintf(stderr, "  5 trees: %.2fs, graph + dag + K_700 + array: %.2fs\n",
          mid - start, stop - mid);
  CHECK(stop - start < 10); // catches quadratic behaviour
  for (auto& t : trees) CHECK(is_tree(n, t));
  CHECK(simple(n, e) && components(n, e) == 1 && acyclic(n, d));
  CHECK_EQ(sz(complete), k * (k - 1) / 2);
  sort(all(a));
  rep(i, 0, n) CHECK_EQ(a[i], (ll)i);
}

int main() {
  test_numbers();
  test_extremes();
  test_seeds();
  test_trees();
  test_graphs();
  test_large();
  puts("OK contest/gen");
}
