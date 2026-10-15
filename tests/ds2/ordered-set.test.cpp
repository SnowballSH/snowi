#include "../testlib.hpp"
#include "../../content/ds2/ordered-set.hpp"

// Set of ints against a sorted vector of distinct values.
void stress_set() {
  rep(it, 0, iters(2000)) {
    Tree<int> t;
    vi ref;
    rep(op, 0, 40) {
      int x = (int)rnd(-10, 10);
      auto pos = lower_bound(all(ref), x);
      bool has = pos != ref.end() && *pos == x;
      if (rnd(0, 2)) {
        CHECK_EQ(t.insert(x).second, !has);
        if (!has) ref.insert(pos, x);
      } else {
        CHECK_EQ((int)t.erase(x), (int)has);
        if (has) ref.erase(pos);
      }
      CHECK_EQ(sz(t), sz(ref));
      int y = (int)rnd(-12, 12);
      CHECK_EQ((int)t.order_of_key(y),
               int(lower_bound(all(ref), y) - ref.begin()));
      int k = (int)rnd(0, sz(ref));
      if (k < sz(ref)) CHECK_EQ(*t.find_by_order(k), ref[k]);
      else CHECK(t.find_by_order(k) == t.end());
    }
  }
}

// The multiset trick from Usage, against a sorted vector.
void stress_multiset() {
  rep(it, 0, iters(2000)) {
    Tree<pii> ms;
    int id = 0;
    vi ref;
    rep(op, 0, 40) {
      int x = (int)rnd(-5, 5);
      if (rnd(0, 2) || ref.empty()) {
        ms.insert({x, id++});
        ref.insert(upper_bound(all(ref), x), x);
      } else {
        x = ref[rnd(0, sz(ref) - 1)];
        ms.erase(ms.lower_bound({x, INT_MIN}));
        ref.erase(lower_bound(all(ref), x));
      }
      int y = (int)rnd(-6, 6);
      int below = int(lower_bound(all(ref), y) - ref.begin());
      int upto = int(upper_bound(all(ref), y) - ref.begin());
      CHECK_EQ((int)ms.order_of_key({y, INT_MIN}), below);
      CHECK_EQ((int)ms.order_of_key({y, INT_MAX}), upto);
      CHECK_EQ(sz(ms), sz(ref));
      rep(k, 0, sz(ref)) CHECK_EQ(ms.find_by_order(k)->first, ref[k]);
    }
  }
}

// The Example: online inversion count.
void check_inversions() {
  rep(it, 0, iters(500)) {
    int n = (int)rnd(0, 30);
    vi a(n);
    for (int& x : a) x = (int)rnd(-4, 4);
    Tree<pii> seen;
    ll inv = 0, want = 0;
    rep(i, 0, n) {
      inv += i - seen.order_of_key({a[i], INT_MAX});
      seen.insert({a[i], i});
    }
    rep(i, 0, n)
      rep(j, 0, i) want += a[j] > a[i];
    CHECK_EQ(inv, want);
  }
}

void check_join_split() {
  Tree<int> t, t2;
  rep(i, 0, 10) t.insert(i), t2.insert(100 + i);
  t.join(t2); // larger keys
  CHECK(t2.empty() && sz(t) == 20);
  rep(i, 0, 5) t2.insert(-10 + i);
  t.join(t2); // smaller keys
  CHECK(t2.empty() && sz(t) == 25 && *t.begin() == -10);
  t.split(104, t2);
  CHECK_EQ(sz(t), 20);
  CHECK_EQ(*t2.begin(), 105);
  CHECK_EQ(*t.find_by_order(19), 104);
  t2.insert(50);
  bool threw = false;
  try {
    t.join(t2);
  } catch (const exception&) { threw = true; }
  CHECK(threw);
}

int main() {
  stress_set();
  stress_multiset();
  check_inversions();
  check_join_split();
  // Large: 1e5 inserts, then erases and order queries against a
  // sorted vector, with extreme keys.
  elapsed();
  int n = 100'000;
  Tree<ll> t;
  vector<ll> v = rnd_vec(n, LLONG_MIN, LLONG_MAX);
  v.push_back(LLONG_MIN), v.push_back(LLONG_MAX);
  for (ll x : v) t.insert(x);
  sort(all(v));
  v.erase(unique(all(v)), v.end());
  CHECK_EQ(sz(t), sz(v));
  rep(i, 0, sz(v)) {
    CHECK_EQ(*t.find_by_order(i), v[i]);
    CHECK_EQ((int)t.order_of_key(v[i]), i);
  }
  rep(i, 0, sz(v))
    if (i % 2) t.erase(v[i]);
  rep(i, 0, sz(t)) CHECK_EQ(*t.find_by_order(i), v[2 * i]);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/ordered-set");
}
