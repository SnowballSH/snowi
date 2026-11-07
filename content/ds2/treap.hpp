/**
 * Title: Implicit treap
 * Description: A sequence with insert and erase at any position,
 *   range reverse, range sum and access by index. `split(x, k)`
 *   returns the first `k` elements and the rest; `merge(a, b)`
 *   concatenates. Node 0 is the empty tree.
 * Precondition: 0-indexed positions, half-open ranges `[l, r)`,
 *   `0 <= i <= size()` for `insert`, `i < size()` otherwise.
 * Time: $O(log n)$ expected per operation. Recursion depth is the
 *   tree height, $O(log n)$ expected (about 60 for $n = 10^6$).
 * Usage:
 *   Treap tr;
 *   tr.insert(i, v); tr.erase(i); tr.get(i);
 *   tr.reverse(l, r); ll s = tr.query(l, r); // sum of [l, r)
 *   vector<ll> out; tr.collect(tr.root, out); // in order
 *   auto [a, b] = tr.split(tr.root, k); // move the first k
 *   tr.root = tr.merge(b, a); // elements to the end
 * Tweak: Other lazy tags (range add, range assign) go next to
 *   `flip`: apply them to `val` and `sum` (add `x * cnt` to `sum`)
 *   when set on a node, pass them to the children in `push`. An
 *   ordered treap (a set) splits by key instead: compare `val`
 *   with the key where `split` compares `k` with the left size.
 * Example:
 *   // Cut [l, r) and paste it at position p of the rest.
 *   auto [a, bc] = tr.split(tr.root, l);
 *   auto [b, c] = tr.split(bc, r - l);
 *   auto [x, y] = tr.split(tr.merge(a, c), p);
 *   tr.root = tr.merge(tr.merge(x, b), y);
 * Notes: Think of it when a sequence is edited in the middle:
 *   block moves, reversals, "insert after the k-th", rope-like
 *   strings. `on(l, r, f)` splits out `[l, r)`, lets `f` read or
 *   replace that subtree, and merges back; new operations are a
 *   few lines on top of it. Node fields hold only aggregates of
 *   the subtree, so a lazy tag must already be applied to the
 *   node that carries it (here `flip` swaps children lazily, which
 *   does not change `sum`). Erased nodes are not reused; the pool
 *   grows by one node per insert. The priorities are seeded from
 *   the clock, so the expected height is $O(log n)$ for every
 *   input; a fixed seed lets a hacker who knows the code build a
 *   path. For only "k-th element" or "insert in sorted order", an
 *   ordered set is simpler.
 * Source: Standard implicit treap (cp-algorithms); kactl Treap.
 * Status: stress-tested against a vector
 */
#pragma once

struct Treap {
  struct Node {
    int l = 0, r = 0, cnt = 0;
    unsigned pri = 0;
    ll val = 0, sum = 0;
    bool flip = 0;
  };
  vector<Node> t{Node()};
  mt19937 rng{(unsigned)chrono::steady_clock::now()
                  .time_since_epoch()
                  .count()};
  int root = 0;
  int make(ll v) {
    t.push_back({0, 0, 1, (unsigned)rng(), v, v, 0});
    return sz(t) - 1;
  }
  void pull(int x) {
    Node &a = t[x], &L = t[a.l], &R = t[a.r];
    a.cnt = L.cnt + 1 + R.cnt, a.sum = L.sum + a.val + R.sum;
  }
  void push(int x) {
    if (!t[x].flip) return;
    swap(t[x].l, t[x].r);
    t[t[x].l].flip ^= 1, t[t[x].r].flip ^= 1, t[x].flip = 0;
  }
  pair<int, int> split(int x, int k) { // first k, rest
    if (!x) return {0, 0};
    push(x);
    int left = t[t[x].l].cnt;
    if (k <= left) {
      auto [a, b] = split(t[x].l, k);
      t[x].l = b, pull(x);
      return {a, x};
    }
    auto [a, b] = split(t[x].r, k - left - 1);
    t[x].r = a, pull(x);
    return {x, b};
  }
  int merge(int a, int b) {
    if (!a || !b) return a ^ b;
    if (t[a].pri > t[b].pri) {
      push(a), t[a].r = merge(t[a].r, b), pull(a);
      return a;
    }
    push(b), t[b].l = merge(a, t[b].l), pull(b);
    return b;
  }
  void on(int l, int r, auto f) { // f(int& subtree of [l, r))
    auto [a, bc] = split(root, l);
    auto [b, c] = split(bc, r - l);
    f(b), root = merge(merge(a, b), c);
  }
  void insert(int i, ll v) {
    on(i, i, [&](int& b) { b = make(v); });
  }
  void erase(int i) {
    on(i, i + 1, [](int& b) { b = 0; });
  }
  void reverse(int l, int r) {
    on(l, r, [&](int& b) { t[b].flip ^= 1; });
  }
  ll query(int l, int r) {
    ll s = 0;
    on(l, r, [&](int& b) { s = t[b].sum; });
    return s;
  }
  ll get(int i) {
    for (int x = root;;) {
      push(x);
      int left = t[t[x].l].cnt;
      if (i == left) return t[x].val;
      if (i < left) x = t[x].l;
      else i -= left + 1, x = t[x].r;
    }
  }
  void collect(int x, vector<ll>& out) {
    if (!x) return;
    push(x), collect(t[x].l, out);
    out.push_back(t[x].val), collect(t[x].r, out);
  }
  int size() { return t[root].cnt; }
};
