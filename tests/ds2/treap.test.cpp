#include "../testlib.hpp"
#include "../../content/ds2/treap.hpp"

vector<ll> contents(Treap& tr) {
  vector<ll> out;
  tr.collect(tr.root, out);
  return out;
}

void stress(int rounds, ll c) {
  rep(it, 0, iters(rounds)) {
    Treap tr;
    vector<ll> a;
    rep(op, 0, 60) {
      int n = sz(a), kind = (int)rnd(0, 6);
      int l = (int)rnd(0, n), r = (int)rnd(0, n);
      if (l > r) swap(l, r);
      if (kind <= 1 || n == 0) {
        ll v = rnd(-c, c);
        tr.insert(l, v), a.insert(a.begin() + l, v);
      } else if (kind == 2) {
        int i = (int)rnd(0, n - 1);
        tr.erase(i), a.erase(a.begin() + i);
      } else if (kind == 3) {
        tr.reverse(l, r), reverse(a.begin() + l, a.begin() + r);
      } else if (kind == 4) {
        CHECK_EQ(tr.query(l, r), accumulate(a.begin() + l, a.begin() + r, 0LL));
      } else if (kind == 5) {
        int i = (int)rnd(0, n - 1);
        CHECK_EQ(tr.get(i), a[i]);
      } else { // Usage: move the first k elements to the end
        auto [x, y] = tr.split(tr.root, l);
        tr.root = tr.merge(y, x);
        rotate(a.begin(), a.begin() + l, a.end());
      }
      CHECK_EQ(tr.size(), sz(a));
    }
    CHECK(contents(tr) == a);
  }
}

// The Example: cut [l, r) and paste it at position p of the rest.
void check_cut_paste() {
  rep(it, 0, iters(1000)) {
    int n = (int)rnd(0, 20);
    Treap tr;
    vector<ll> a(n);
    rep(i, 0, n) a[i] = rnd(0, 9), tr.insert(i, a[i]);
    if (rnd(0, 1)) {
      int l = (int)rnd(0, n), r = (int)rnd(0, n);
      if (l > r) swap(l, r);
      tr.reverse(l, r), reverse(a.begin() + l, a.begin() + r);
    }
    int l = (int)rnd(0, n), r = (int)rnd(l, n);
    int p = (int)rnd(0, n - (r - l));
    auto [x0, bc] = tr.split(tr.root, l);
    auto [b, c] = tr.split(bc, r - l);
    auto [x, y] = tr.split(tr.merge(x0, c), p);
    tr.root = tr.merge(tr.merge(x, b), y);
    vector<ll> mid(a.begin() + l, a.begin() + r), rest = a;
    rest.erase(rest.begin() + l, rest.begin() + r);
    rest.insert(rest.begin() + p, all(mid));
    CHECK(contents(tr) == rest);
  }
}

int height(Treap& tr) {
  int h = 0;
  vector<pii> stack{{tr.root, 0}};
  while (sz(stack)) {
    auto [x, d] = stack.back();
    stack.pop_back();
    if (!x) continue;
    h = max(h, d + 1);
    stack.push_back({tr.t[x].l, d + 1});
    stack.push_back({tr.t[x].r, d + 1});
  }
  return h;
}

// Regression: an adversary who knows the code inserts each node at
// the rank of the priority a default-seeded mt19937 gives it. With
// that seed the in-order priorities are sorted: a path of height n.
void check_unpredictable() {
  int n = 3000;
  mt19937 guess;
  vector<unsigned> seen;
  Treap tr;
  rep(i, 0, n) {
    unsigned p = guess();
    auto pos = lower_bound(all(seen), p);
    tr.insert(int(pos - seen.begin()), i);
    seen.insert(pos, p);
  }
  CHECK(height(tr) < 100);
}

int main() {
  check_unpredictable();
  stress(2500, 9);
  stress(1000, 1'000'000'000'000'000);
  check_cut_paste();
  { // empty ranges and a single element
    Treap tr;
    CHECK_EQ(tr.query(0, 0), 0LL);
    tr.reverse(0, 0);
    tr.insert(0, 5);
    tr.reverse(0, 1), tr.reverse(1, 1);
    CHECK_EQ(tr.get(0), 5LL);
    CHECK_EQ(tr.query(0, 1), 5LL);
    tr.erase(0);
    CHECK_EQ(tr.size(), 0);
    CHECK(contents(tr).empty());
  }
  // Large: 1e5 inserts at random positions, then 1e5 reversals,
  // each checked against the sums and ends it must preserve.
  elapsed();
  int n = 100'000;
  Treap tr;
  rep(i, 0, n) tr.insert((int)rnd(0, i), i);
  rep(it, 0, n) {
    int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
    ll s = tr.query(l, r), first = tr.get(l), last = tr.get(r - 1);
    tr.reverse(l, r);
    CHECK_EQ(tr.query(l, r), s);
    CHECK_EQ(tr.get(l), last);
    CHECK_EQ(tr.get(r - 1), first);
  }
  vector<ll> out = contents(tr);
  sort(all(out));
  rep(i, 0, n) CHECK_EQ(out[i], (ll)i);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/treap");
}
