/**
 * Title: Link-cut tree
 * Description: A forest under edge insertions and deletions, with
 *   path queries, point updates. Each tree has a root, and
 *   `make_root(v)` moves it to `v` (evert). `link(u, v)` adds the
 *   edge `u - v` with `u` as a child of `v`; the root of `v`'s tree
 *   stays the root. `cut(u)` removes the edge from `u` to its parent
 *   (nothing at a root); `cut(u, v)` removes the edge `u - v` and
 *   leaves `v` as the root. `query(u, v)` folds the values on the
 *   path from `u` to `v` and leaves `u` as the root. `lca(u, v)` and
 *   `find_root(v)` use the current root. After `access(v)`,
 *   `t[v].agg` is the fold of the root path of `v`.
 * Precondition: 0-indexed. `link`: `u` and `v` in different trees.
 *   `cut(u, v)`: the edge exists. `lca`, `query`: `u` and `v`
 *   connected. `op` is associative and commutative (paths get
 *   reversed) with identity `e`.
 * Time: Amortised $O(log n)$ per operation.
 * Usage:
 *   LinkCut lc(a); // vertex values a[v], no edges
 *   lc.link(u, v); lc.cut(u, v); lc.connected(u, v);
 *   lc.query(u, v); lc.set(v, x); lc.lca(u, v);
 *   lc.access(v); ll s = lc.t[v].agg; // root path, same root
 * Tweak: Change only `T`, `e` and `op` (default: path sum). Edge
 *   values: add a node per edge holding its value (vertices hold
 *   `e`) and link `u - x - v`.
 *   Non-commutative `op`: keep the reverse-order fold `ragg` next
 *   to `agg` everywhere (`pull` combines `c[1]`, `val`, `c[0]`)
 *   and swap a node's `agg` and `ragg` the moment its `flip` is
 *   toggled (each child in `push`, `x` in `make_root`), not when
 *   it is pushed: `pull` reads children with pending flips.
 * Example:
 *   // Minimum spanning forest under edge insertions. T is
 *   // pair<ll, int> (weight, edge id), op = max, e = {LLONG_MIN,
 *   // -1}. Node n + i is edge i = (a[i], b[i], w[i]).
 *   LinkCut lc(vector<LinkCut::T>(n + m, LinkCut::e));
 *   ll total = 0;
 *   auto add = [&](int i) {
 *     int u = a[i], v = b[i];
 *     if (lc.connected(u, v)) {
 *       auto [mw, j] = lc.query(u, v); // heaviest edge on path
 *       if (mw <= w[i]) return;
 *       lc.cut(n + j, a[j]), lc.cut(n + j, b[j]), total -= mw;
 *     }
 *     lc.set(n + i, {w[i], i});
 *     lc.link(n + i, u), lc.link(n + i, v), total += w[i];
 *   };
 * Notes: Cues: edges added and removed online with path questions
 *   ("are `u` and `v` connected now", "max weight on the path"), a
 *   spanning forest that must evict the heaviest or oldest edge of
 *   the cycle a new edge closes, parent pointers that change (a
 *   rooted forest where a subtree is moved under a new parent:
 *   `cut(u)`, then `link(u, p)`). For a fixed root, call
 *   `link(u, p)` only on a root `u` (after `cut(u)`) and never call
 *   `make_root`, `query` or `cut(u, v)`; use `access(v)` and
 *   `t[v].agg` for root paths.
 *   Moving a subtree changes the depths of all its vertices at once:
 *   put that change on the new edge (a node between `u` and `p`, or
 *   the value of `u` itself) instead of on the vertices. To find the
 *   first vertex on a root path with a property of its fold (value
 *   at most `x`), `access(v)` and walk down from `v` using `agg` of
 *   the children, pushing as you go, then `splay` the vertex found.
 *   Offline, with all edge changes known in advance and only
 *   connectivity or per-component questions,
 *   `offline/dynamic-connectivity` (DSU with rollback) is simpler;
 *   it knows no paths, so path queries, `lca` or a moving root still
 *   need this.
 *   Subtree folds need each node to also keep the fold of its
 *   virtual (path-parent) children. Speed with `-O2`: $10^5$ mixed
 *   operations take about 0.15 s; $10^6$ on $10^6$ vertices about 3
 *   s, mostly cache misses.
 * Source: Sleator and Tarjan; array splay as in kactl LinkCutTree
 *   (CC0).
 * Status: stress-tested against an explicit forest
 */
#pragma once

struct LinkCut {
  typedef ll T;
  static constexpr T e = 0;
  T op(T a, T b) { return a + b; }
  struct Node {
    int c[2], p, flip;
    T val, agg;
  };
  int n; // node n is null and holds e
  vector<Node> t;
  vi st;
  LinkCut(const vector<T>& a)
      : n(sz(a)), t(n + 1, {{n, n}, n, 0, e, e}) {
    rep(i, 0, n) t[i].val = t[i].agg = a[i];
  }
  bool top(int x) {
    auto& c = t[t[x].p].c;
    return c[0] != x && c[1] != x;
  }
  void pull(int x) {
    Node& u = t[x];
    u.agg = op(op(t[u.c[0]].agg, u.val), t[u.c[1]].agg);
  }
  void push(int x) {
    if (t[x].flip) {
      swap(t[x].c[0], t[x].c[1]), t[x].flip = 0;
      for (int c : t[x].c) t[c].flip ^= 1;
    }
  }
  void rot(int x) {
    int y = t[x].p, z = t[y].p, d = t[y].c[1] == x,
        b = t[x].c[!d];
    if (!top(y)) t[z].c[t[z].c[1] == y] = x;
    t[x].p = z, t[y].c[d] = b, t[b].p = y;
    t[x].c[!d] = y, t[y].p = x;
    pull(y), pull(x);
  }
  void splay(int x) {
    st = {x};
    for (int y = x; !top(y); y = t[y].p) st.push_back(t[y].p);
    for (int i = sz(st); i--;) push(st[i]);
    while (!top(x)) {
      int y = t[x].p, z = t[y].p;
      if (!top(y))
        rot((t[y].c[0] == x) == (t[z].c[0] == y) ? y : x);
      rot(x);
    }
  }
  // Makes root..x one path. Returns the vertex where it met the
  // root's old path: the lca with the vertex accessed before.
  int access(int x) {
    int last = n;
    for (int y = x; y != n; last = y, y = t[y].p)
      splay(y), t[y].c[1] = last, pull(y);
    splay(x);
    return last;
  }
  void make_root(int x) { access(x), t[x].flip ^= 1; }
  int find_root(int x) {
    access(x);
    while (push(x), t[x].c[0] != n) x = t[x].c[0];
    splay(x);
    return x;
  }
  bool connected(int u, int v) {
    return find_root(u) == find_root(v);
  }
  void link(int u, int v) { make_root(u), t[u].p = v; }
  void cut(int u) { // u is last on its path: no right child
    access(u);
    t[t[u].c[0]].p = n, t[u].c[0] = n, t[u].agg = t[u].val;
  }
  void cut(int u, int v) { make_root(v), cut(u); }
  T query(int u, int v) {
    make_root(u), access(v);
    return t[v].agg;
  }
  void set(int v, T x) { splay(v), t[v].val = x, pull(v); }
  int lca(int u, int v) {
    access(u);
    return access(v);
  }
};
