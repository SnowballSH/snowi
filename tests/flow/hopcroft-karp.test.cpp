#include "../testlib.hpp"
#include "../../content/flow/hopcroft-karp.hpp"

// Maximum matching by DP over left vertices and used right sets.
int brute(int n, int m, const vector<vi>& adj) {
  vi f(1 << m, 0), nf;
  for (int u = n - 1; u >= 0; u--) {
    nf = f;
    rep(mask, 0, 1 << m)
      for (int v : adj[u])
        if (!(mask >> v & 1)) nf[mask] = max(nf[mask], 1 + f[mask | 1 << v]);
    f = nf;
  }
  return f[0];
}

// A valid matching of size k plus a vertex cover of size k.
void certify(HopcroftKarp& hk, int n, int m, const vector<pii>& edges, int k) {
  set<pii> es(all(edges));
  int matched = 0;
  rep(u, 0, n)
    if (hk.ml[u] >= 0) {
      CHECK(es.count({u, hk.ml[u]}));
      CHECK_EQ(hk.mr[hk.ml[u]], u);
      matched++;
    }
  rep(v, 0, m)
    if (hk.mr[v] >= 0) CHECK_EQ(hk.ml[hk.mr[v]], v);
  CHECK_EQ(matched, k);
  auto [cl, cr] = hk.cover();
  CHECK_EQ(sz(cl) + sz(cr), k);
  vector<bool> inl(n), inr(m);
  for (int u : cl) inl[u] = true;
  for (int v : cr) inr[v] = true;
  for (auto [u, v] : edges) CHECK(inl[u] || inr[v]);
}

int main() {
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(0, 8), m = (int)rnd(0, 8);
    int e = n && m ? (int)rnd(0, 20) : 0;
    vector<pii> edges;
    vector<vi> adj(n);
    HopcroftKarp hk(n, m);
    rep(i, 0, e) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, m - 1);
      edges.push_back({u, v}), adj[u].push_back(v);
      hk.add_edge(u, v);
    }
    int k = hk.max_matching();
    CHECK_EQ(k, brute(n, m, adj));
    certify(hk, n, m, edges, k);
    CHECK_EQ(hk.max_matching(), k); // nothing left to augment
    certify(hk, n, m, edges, k);
    if (n && m) { // Tweak: add edges later, continue
      rep(i, 0, 3) {
        int u = (int)rnd(0, n - 1), v = (int)rnd(0, m - 1);
        edges.push_back({u, v}), adj[u].push_back(v);
        hk.add_edge(u, v);
      }
      k = hk.max_matching();
      CHECK_EQ(k, brute(n, m, adj));
      certify(hk, n, m, edges, k);
    }
  }
  { // Diamond: left A_i, B_i matched to a_i, b_i, both joined to
    // a_(i+1) and b_(i+1); a free root reaches a_0, b_0 and the
    // last layer is a dead end. A chain C_j - c_j of length L + 1
    // sets the phase length. Without visit marks: 2^L searches.
    int L = 60, n = 3 * L + 3, m = 3 * L + 2, r = 3 * L + 1;
    vector<pii> edges;
    rep(i, 0, L)
      for (int u : {2 * i, 2 * i + 1}) {
        edges.push_back({u, u});
        if (i + 1 < L)
          edges.push_back({u, 2 * i + 2}), edges.push_back({u, 2 * i + 3});
      }
    rep(j, 0, L + 1)
      edges.push_back({2 * L + j, 2 * L + j}),
          edges.push_back({2 * L + j, 2 * L + j + 1});
    edges.push_back({r, 0}), edges.push_back({r, 1});
    edges.push_back({r + 1, 2 * L});
    HopcroftKarp hk(n, m);
    for (auto [u, v] : edges) hk.add_edge(u, v);
    CHECK_EQ(hk.max_matching(), m);
    certify(hk, n, m, edges, m);
  }
  elapsed();
  { // Large random graph; the cover certifies optimality.
    int n = 50000, m = 40000;
    vector<pii> edges;
    HopcroftKarp hk(n, m);
    rep(i, 0, 150000) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, m - 1);
      edges.push_back({u, v}), hk.add_edge(u, v);
    }
    int k = hk.max_matching();
    certify(hk, n, m, edges, k);
  }
  { // Chain: the last augmenting path has length 2n.
    int n = 20000;
    vector<pii> edges;
    HopcroftKarp hk(n, n);
    rep(u, 0, n) {
      if (u + 1 < n) edges.push_back({u, u + 1}), hk.add_edge(u, u + 1);
      edges.push_back({u, u}), hk.add_edge(u, u);
    }
    CHECK_EQ(hk.max_matching(), n);
    certify(hk, n, n, edges, n);
  }
  fprintf(stderr, "  large cases: %.2fs\n", elapsed());
  puts("OK flow/hopcroft-karp");
}
