/**
 * Title: Random test generators
 * Description: Seeded random inputs for stress tests. Trees and
 *   graphs are 0-indexed edge lists, relabelled by a random
 *   permutation, shuffled and randomly oriented (`mix(e, n)` does
 *   this to any edge list); `dag` orients them along a random
 *   topological order. A connected `graph` contains a
 *   `tree(n, shape)`. Shapes: `UNIFORM` (all labelled trees
 *   equally likely), `PARENT` (each vertex below a random earlier
 *   one), `PATH`, `STAR`, `CATERPILLAR` (a path with leaves on it).
 * Precondition: `tree`, `graph`: `n >= 1`. `array(n, lo, hi,
 *   true)`: `hi - lo + 1 >= n`. `graph(n, m)`: `m <= n(n-1)/2`,
 *   and `m >= n - 1` if connected.
 * Time: $O((n + m) log n)$.
 * Usage:
 *   int main(int, char** argv) {
 *     Gen g(atoll(argv[1]));
 *     int n = g.rnd(1, 6), m = g.rnd(n - 1, n * (n - 1) / 2);
 *     printf("%d %d\n", n, m);
 *     for (auto [u, v] : g.graph(n, m)) // weight: drop if unused
 *       printf("%d %d %lld\n", u + 1, v + 1, g.rnd(1, 9));
 *   }
 * Example:
 *   auto t = g.tree(n, Gen::CATERPILLAR); // or PATH, STAR, ...
 *   auto deep = g.graph(n, m, true, Gen::PATH); // long paths
 *   auto d = g.dag(n, m);                 // edges u -> v
 *   vector<ll> a = g.array(n, 1, 3);      // small range: ties
 *   vector<ll> b = g.array(n, 1, 1e9, true); // distinct
 *   vi p = g.perm(n);
 *   string s = g.str(n, "abc");
 *   vector<pii> e; // own shape, e.g. a binary tree, then mix:
 *   rep(i, 1, n) e.emplace_back(i, (i - 1) / 2);
 *   e = g.mix(e, n);
 *   int T = g.rnd(1, 3); // multi-test: print T, then T cases
 * Notes: Rooted at vertex 0, `PARENT` and `UNIFORM` trees with
 *   $n = 2 dot 10^5$ have depth only about 40 and 1100: they miss
 *   deep-recursion crashes and the worst cases of HLD and
 *   small-to-large merging. Test `PATH` and `CATERPILLAR` too,
 *   also as spanning trees of `graph`. Small value ranges and
 *   alphabets force ties, hash collisions and periodic strings.
 *   Keep `n` tiny (about 8) so the brute force is fast and a
 *   failing input is short to read. Multi-test input: print `T`,
 *   then each case from the same `g`, keeping the sum of sizes
 *   within the limits.
 * Status: stress-tested: tree shapes, simple, connected and acyclic
 *   graphs, 64-bit ranges, chi-square uniformity of trees,
 *   spanning trees, letters, edge sets and orders
 */
#pragma once

struct Gen {
  mt19937_64 rng;
  Gen(ll seed) : rng(seed) {}
  ll rnd(ll lo, ll hi) {
    return uniform_int_distribution(lo, hi)(rng);
  }
  vector<ll> array(int n, ll lo, ll hi, bool distinct = false) {
    vector<ll> a;
    unordered_set<ll> s; // Floyd's sampling without replacement
    rep(i, 0, n) {
      ll j = distinct ? hi - (n - 1 - i) : hi, x = rnd(lo, j);
      if (distinct && !s.insert(x).second) s.insert(x = j);
      a.push_back(x);
    }
    shuffle(all(a), rng);
    return a;
  }
  vi perm(int n) {
    vi p(n);
    iota(all(p), 0), shuffle(all(p), rng);
    return p;
  }
  // Relabel, shuffle; orient randomly (acyclic: by old label).
  vector<pii> mix(vector<pii> e, int n, bool acyclic = false) {
    vi lab = perm(n);
    for (auto& [u, v] : e) {
      if (acyclic ? u > v : rng() & 1) swap(u, v);
      u = lab[u], v = lab[v];
    }
    shuffle(all(e), rng);
    return e;
  }
  enum { UNIFORM, PARENT, PATH, STAR, CATERPILLAR };
  vector<pii> tree(int n, int shape = UNIFORM) {
    vector<pii> e;
    if (shape == UNIFORM) { // random Pruefer code, then n - 1
      vi code(n - 1, n - 1), deg(n, 1);
      rep(i, 0, n - 2) deg[code[i] = rnd(0, n - 1)]++;
      priority_queue<int, vi, greater<>> q;
      rep(i, 0, n)
        if (deg[i] == 1) q.push(i);
      for (int v : code) {
        e.emplace_back(q.top(), v), q.pop();
        if (--deg[v] == 1) q.push(v);
      }
    } else { // path 0..k-1, the others below earlier vertices
      int k = shape == PATH ? n : 1;
      if (shape == CATERPILLAR) k = rnd(1, n);
      rep(i, 1, n) {
        int below = shape == PARENT ? i : k;
        e.emplace_back(i, i < k ? i - 1 : rnd(0, below - 1));
      }
    }
    return mix(e, n);
  }
  vector<pii> graph(int n, int m, bool connected = true,
                    int shape = UNIFORM) {
    auto e = connected ? tree(n, shape) : vector<pii>();
    unordered_set<ll> s(2 * n); // tree edges as u * n + v
    for (auto [u, v] : e)
      s.insert((ll)u * n + v), s.insert((ll)v * n + u);
    // pair t = {u, u + d mod n}, u = t % n, d = t / n + 1
    for (ll t : array(m, 0, (ll)n * (n - 1) / 2 - 1, true)) {
      int u = t % n, v = (u + t / n + 1) % n;
      if (sz(e) < m && !s.count((ll)u * n + v))
        e.emplace_back(u, v);
    }
    return mix(e, n);
  }
  vector<pii> dag(int n, int m) {
    return mix(graph(n, m, false), n, true);
  }
  string str(int n, string alphabet = "ab") {
    string s(n, ' ');
    for (char& c : s) c = alphabet[rnd(0, sz(alphabet) - 1)];
    return s;
  }
};
