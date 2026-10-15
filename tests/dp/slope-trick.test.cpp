#include "../testlib.hpp"
#include "../../content/dp/slope-trick.hpp"

// f tabulated on [-G, G]; breakpoints stay well inside, so f is
// linear beyond the grid and value() extrapolates it exactly.
struct Grid {
  int G;
  vector<ll> v;
  Grid(int g) : G(g), v(2 * g + 1) {}
  ll& at(ll x) { return v[x + G]; }
  ll value(ll x) {
    if (x < -G) return at(-G) + (at(-G + 1) - at(-G)) * (x + G);
    if (x > G) return at(G) + (at(G) - at(G - 1)) * (x - G);
    return at(x);
  }
  template<class F> void add(F f) {
    rep(x, -G, G + 1) at(x) += f(x);
  }
  void prefix_min() {
    rep(x, -G + 1, G + 1) at(x) = min(at(x), at(x - 1));
  }
  void suffix_min() {
    for (int x = G - 1; x >= -G; x--) at(x) = min(at(x), at(x + 1));
  }
  void slide(ll a, ll b) {
    vector<ll> w(sz(v));
    rep(x, -G, G + 1) {
      w[x + G] = LLONG_MAX;
      for (ll y = x - b; y <= x - a; y++) w[x + G] = min(w[x + G], value(y));
    }
    v = w;
  }
};

void compare(SlopeTrick& st, Grid& g, int W) {
  CHECK(st.top_l() <= st.top_r());
  if (W <= 30)
    rep(x, -W, W + 1) CHECK_EQ(st.eval(x), g.at(x));
  else
    rep(t, 0, 20) {
      ll x = rnd(-W, W);
      CHECK_EQ(st.eval(x), g.at(x));
    }
  ll mn = *min_element(all(g.v));
  CHECK_EQ(st.min_f, mn);
  int lo = -g.G, hi = g.G;
  while (g.at(lo) != mn) lo++;
  while (g.at(hi) != mn) hi--;
  // An argmin reaching the grid edge with slope 0 is unbounded.
  if (lo == -g.G && g.at(-g.G + 1) == mn)
    CHECK_EQ(st.top_l(), -SlopeTrick::INF);
  else CHECK_EQ(st.top_l(), (ll)lo);
  if (hi == g.G && g.at(g.G - 1) == mn) CHECK_EQ(st.top_r(), SlopeTrick::INF);
  else CHECK_EQ(st.top_r(), (ll)hi);
}

void stress(int rounds, int ops, ll c, ll d) {
  rep(it, 0, iters(rounds)) {
    SlopeTrick st;
    Grid g(int(c + ops * 2 * d + 40));
    rep(op, 0, ops) {
      ll a = rnd(-c, c);
      switch (rnd(0, 7)) {
      case 0: {
        ll k = rnd(-5, 5);
        st.add_const(k), g.add([&](ll) { return k; });
        break;
      }
      case 1:
        st.add_x_minus_a(a);
        g.add([&](ll x) { return max(0LL, x - a); });
        break;
      case 2:
        st.add_a_minus_x(a);
        g.add([&](ll x) { return max(0LL, a - x); });
        break;
      case 3:
      case 4: st.add_abs(a), g.add([&](ll x) { return abs(x - a); }); break;
      case 5: st.prefix_min(), g.prefix_min(); break;
      case 6: st.suffix_min(), g.suffix_min(); break;
      default: {
        ll lo = rnd(-d, d), hi = rnd(lo, d);
        st.slide(lo, hi), g.slide(lo, hi);
      }
      }
      compare(st, g, int(c + (op + 1) * d + 3));
    }
  }
}

// Fewest +-1 steps making a non-decreasing, by DP over values.
ll brute_isotonic(const vector<ll>& a, ll V) {
  vector<ll> dp(V + 1, 0);
  for (ll x : a) {
    rep(v, 1, V + 1) dp[v] = min(dp[v], dp[v - 1]);
    rep(v, 0, V + 1) dp[v] += abs(x - v);
  }
  return *min_element(all(dp));
}

// Regression for the Example's variants: fewest +-1 steps so that
// |b[i] - b[i-1]| <= d (slide(-d, d) in place of prefix_min), or
// 0 <= b[i] - b[i-1] <= d (slide(0, d)), against every b in
// [min a, max a]^n (clamping b there keeps it feasible).
void bounded_steps(int rounds) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(1, 4), V = 6;
    vector<ll> a = rnd_vec(n, 0, V);
    ll d = rnd(0, 2);
    for (ll lo : {-d, 0LL}) {
      SlopeTrick st;
      for (ll x : a) st.slide(lo, d), st.add_abs(x);
      ll want = LLONG_MAX;
      vector<ll> b(n);
      int total = 1;
      rep(i, 0, n) total *= V + 1;
      rep(code, 0, total) {
        ll c = code, cost = 0;
        bool ok = true;
        rep(i, 0, n) {
          b[i] = c % (V + 1), c /= V + 1;
          cost += abs(b[i] - a[i]);
          if (i && (b[i] - b[i - 1] < lo || b[i] - b[i - 1] > d)) ok = false;
        }
        if (ok) want = min(want, cost);
      }
      CHECK_EQ(st.min_f, want);
    }
  }
}

int main() {
  stress(1500, 12, 3, 1);
  stress(1000, 25, 10, 3);
  bounded_steps(1000);
  // Extremes: only hinges, checked against the explicit sum at
  // every breakpoint (the minimum of a convex function is at one).
  rep(it, 0, iters(300)) {
    SlopeTrick st;
    ll e = 100'000'000'000'000'000;
    vector<pair<int, ll>> hs;
    rep(op, 0, 20) {
      int t = (int)rnd(0, 2);
      ll a = rnd(-e, e);
      hs.push_back({t, a});
      if (t == 0) st.add_x_minus_a(a);
      if (t == 1) st.add_a_minus_x(a);
      if (t == 2) st.add_abs(a);
    }
    auto f = [&](ll x) {
      ll s = 0;
      for (auto [t, a] : hs) {
        if (t != 1) s += max(0LL, x - a);
        if (t != 0) s += max(0LL, a - x);
      }
      return s;
    };
    ll mn = LLONG_MAX;
    for (auto [t, a] : hs) {
      mn = min(mn, f(a));
      CHECK_EQ(st.eval(a), f(a));
    }
    CHECK_EQ(st.min_f, mn);
    if (st.top_l() > -SlopeTrick::INF) CHECK_EQ(f(st.top_l()), mn);
    if (st.top_r() < SlopeTrick::INF) CHECK_EQ(f(st.top_r()), mn);
  }
  { // Empty function, and the Example.
    SlopeTrick st;
    CHECK_EQ(st.min_f, 0LL);
    CHECK_EQ(st.top_l(), -SlopeTrick::INF);
    CHECK_EQ(st.top_r(), SlopeTrick::INF);
    st.prefix_min(), st.suffix_min(), st.slide(-3, 5);
    CHECK_EQ(st.eval(12345), 0LL);
    for (ll x : {5, 2, 4, 1}) st.prefix_min(), st.add_abs(x);
    CHECK_EQ(st.min_f, 6LL);
  }
  // Large: isotonic regression against the DP over values, and a
  // long random mix of all operations for speed.
  elapsed();
  int n = 100'000;
  vector<ll> a = rnd_vec(n, 0, 30);
  SlopeTrick st;
  for (ll x : a) st.prefix_min(), st.add_abs(x);
  CHECK_EQ(st.min_f, brute_isotonic(a, 30));
  SlopeTrick mix;
  rep(i, 0, n) {
    ll v = rnd(-1'000'000'000, 1'000'000'000);
    switch (rnd(0, 9)) {
    case 0: mix.prefix_min(); break;
    case 1: mix.suffix_min(); break;
    case 2: mix.slide(-abs(v) / 1000, abs(v) / 1000); break;
    default: mix.add_abs(v);
    }
  }
  if (mix.top_l() > -SlopeTrick::INF)
    CHECK_EQ(mix.eval(mix.top_l()), mix.min_f);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK dp/slope-trick");
}
