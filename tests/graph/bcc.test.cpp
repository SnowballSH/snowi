#include "../testlib.hpp"
#include "../../content/graph/bcc.hpp"

struct Brute {
  int n;
  vector<pii> edges;
  // Connected components avoiding vertex `skip_v` and edge `skip_e`.
  vi label(int skip_v, int skip_e) {
    vi lab(n, -1);
    rep(s, 0, n)
      if (s != skip_v && lab[s] < 0) {
        vi stack = {s};
        lab[s] = s;
        while (sz(stack)) {
          int u = stack.back();
          stack.pop_back();
          rep(e, 0, sz(edges)) {
            auto [a, b] = edges[e];
            if (e == skip_e || a == skip_v || b == skip_v) continue;
            for (int v : {a == u ? b : -1, b == u ? a : -1})
              if (v >= 0 && lab[v] < 0) lab[v] = s, stack.push_back(v);
          }
        }
      }
    return lab;
  }
  bool bridge(int e) {
    auto [u, v] = edges[e];
    if (u == v) return false;
    vi lab = label(-1, e);
    return lab[u] != lab[v];
  }
  bool cut(int v) {
    vi lab = label(v, -1), nb;
    for (auto [a, b] : edges)
      if (a != b && (a == v || b == v)) nb.push_back(a ^ b ^ v);
    for (int a : nb)
      for (int b : nb)
        if (lab[a] != lab[b]) return true;
    return false;
  }
  // Union the edges of every simple cycle (length >= 2).
  vi cycle_classes() {
    int m = sz(edges);
    vi rep_(m);
    iota(all(rep_), 0);
    function<int(int)> find = [&](int x) {
      return rep_[x] == x ? x : rep_[x] = find(rep_[x]);
    };
    vi path_e, on(n);
    function<void(int, int)> dfs = [&](int s, int u) {
      rep(e, 0, m) {
        auto [a, b] = edges[e];
        if (a == b || (a != u && b != u)) continue;
        if (count(all(path_e), e)) continue;
        int v = a ^ b ^ u;
        if (v == s && sz(path_e) >= 1) {
          for (int f : path_e) rep_[find(f)] = find(e);
        } else if (!on[v]) {
          on[v] = 1, path_e.push_back(e);
          dfs(s, v);
          on[v] = 0, path_e.pop_back();
        }
      }
    };
    rep(s, 0, n) on[s] = 1, dfs(s, s), on[s] = 0;
    vi cls(m);
    rep(e, 0, m) cls[e] = find(e);
    return cls;
  }
};

void check_tree(int n, const vector<pii>& edges, BCC& b) {
  CHECK_EQ(sz(b.tree), n + b.nb);
  vector<set<int>> want(n + b.nb);
  rep(e, 0, sz(edges))
    if (b.comp[e] >= 0) {
      auto [u, v] = edges[e];
      for (int x : {u, v})
        want[x].insert(n + b.comp[e]), want[n + b.comp[e]].insert(x);
    }
  rep(x, 0, n + b.nb) {
    CHECK_EQ(sz(b.tree[x]), sz(want[x]));
    CHECK(set<int>(all(b.tree[x])) == want[x]);
  }
  rep(v, 0, n) CHECK_EQ((bool)b.cut[v], sz(b.tree[v]) >= 2);
}

// Notes: the vertices whose removal separates connected u != v are
// exactly the vertex nodes strictly inside their tree path.
void check_separators(int n, Brute& br, BCC& b) {
  vi lab = br.label(-1, -1);
  vector<vi> without(n);
  rep(x, 0, n) without[x] = br.label(x, -1);
  rep(u, 0, n) {
    vi par(sz(b.tree), -2), q = {u};
    par[u] = -1;
    rep(i, 0, sz(q))
      for (int y : b.tree[q[i]])
        if (par[y] == -2) par[y] = q[i], q.push_back(y);
    rep(v, 0, n) {
      if (v == u || lab[u] != lab[v]) continue;
      CHECK(par[v] >= 0);
      vector<bool> inner(n), sep(n);
      for (int y = par[v]; y != u; y = par[y])
        if (y < n) inner[y] = true;
      rep(x, 0, n) sep[x] = x != u && x != v && without[x][u] != without[x][v];
      CHECK(sep == inner);
    }
  }
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 7), m = (int)rnd(0, 9);
    vector<pii> edges;
    rep(i, 0, m) edges.push_back({(int)rnd(0, n - 1), (int)rnd(0, n - 1)});
    BCC b(n, edges);
    Brute br{n, edges};
    vi cls = br.cycle_classes();
    set<int> blocks;
    rep(e, 0, m) {
      CHECK_EQ((bool)b.bridge[e], br.bridge(e));
      if (edges[e].first == edges[e].second) {
        CHECK_EQ(b.comp[e], -1);
        continue;
      }
      CHECK(0 <= b.comp[e] && b.comp[e] < b.nb);
      blocks.insert(b.comp[e]);
      rep(f, 0, m)
        if (edges[f].first != edges[f].second)
          CHECK_EQ(b.comp[e] == b.comp[f], cls[e] == cls[f]);
    }
    CHECK_EQ(sz(blocks), b.nb);
    rep(v, 0, n) CHECK_EQ((bool)b.cut[v], br.cut(v));
    check_tree(n, edges, b);
    check_separators(n, br, b);
  }
  // Path 0-1-2-3: only 2 separates 1 and 3, though 1 is a cut vertex.
  {
    vector<pii> p4 = {{0, 1}, {1, 2}, {2, 3}};
    BCC b(4, p4);
    Brute br{4, p4};
    check_separators(4, br, b);
  }
  // Large: k cycles of length len chained by bridges, with an extra
  // parallel edge in each cycle and a self-loop; labels shuffled.
  int k = 2000, len = 50, n = k * len;
  vi lab(n);
  iota(all(lab), 0);
  shuffle(all(lab), test_rng);
  vector<pii> edges;
  vi kind; // -1 loop, -2 bridge, else cycle id
  rep(c, 0, k) {
    rep(i, 0, len) {
      int u = c * len + i, v = c * len + (i + 1) % len;
      edges.push_back({lab[u], lab[v]}), kind.push_back(c);
    }
    edges.push_back({lab[c * len], lab[c * len + 1]}), kind.push_back(c);
    edges.push_back({lab[c * len + 2], lab[c * len + 2]});
    kind.push_back(-1);
    if (c + 1 < k) {
      edges.push_back({lab[c * len + len - 1], lab[(c + 1) * len]});
      kind.push_back(-2);
    }
  }
  vi order(sz(edges));
  iota(all(order), 0);
  shuffle(all(order), test_rng);
  vector<pii> shuffled;
  vi skind;
  for (int e : order) shuffled.push_back(edges[e]), skind.push_back(kind[e]);
  BCC b(n, shuffled);
  CHECK_EQ(b.nb, 2 * k - 1);
  map<int, int> block_of_cycle;
  rep(e, 0, sz(shuffled)) {
    CHECK_EQ((bool)b.bridge[e], skind[e] == -2);
    if (skind[e] == -1) CHECK_EQ(b.comp[e], -1);
    if (skind[e] >= 0) {
      auto [p, fresh] = block_of_cycle.insert({skind[e], b.comp[e]});
      CHECK(fresh || p->second == b.comp[e]);
    }
  }
  CHECK_EQ(sz(block_of_cycle), k);
  int cuts = 0;
  rep(v, 0, n) cuts += b.cut[v];
  CHECK_EQ(cuts, 2 * (k - 1));
  check_tree(n, shuffled, b);
  // A path: every edge a bridge, every inner vertex a cut vertex.
  vector<pii> path;
  rep(v, 0, n - 1) path.push_back({v, v + 1});
  BCC p(n, path);
  CHECK_EQ(p.nb, n - 1);
  rep(v, 0, n) CHECK_EQ((bool)p.cut[v], v > 0 && v < n - 1);
  puts("OK graph/bcc");
}
