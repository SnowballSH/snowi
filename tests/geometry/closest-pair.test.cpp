#include "../testlib.hpp"
#include "../../content/geometry/closest-pair.hpp"

// ASan's fake stacks make each by-value call of a small struct
// about 1000x slower on this toolchain; geometry code is all such
// calls and never returns pointers to locals.
extern "C" const char* __asan_default_options() {
  return "detect_stack_use_after_return=0";
}

typedef __int128 LL;
typedef P<ll> Pl;

Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
LL d2(Pl a, Pl b) {
  LL dx = a.x - b.x, dy = a.y - b.y;
  return dx * dx + dy * dy;
}
LL brute(const vector<Pl>& v) {
  LL best = -1;
  rep(i, 0, sz(v))
    rep(j, 0, i)
      if (best < 0 || d2(v[i], v[j]) < best) best = d2(v[i], v[j]);
  return best;
}
// The index recipe from the Tweak: find a, then b at another index.
pii indices(const vector<Pl>& v, Pl a, Pl b) {
  int i = int(find(all(v), a) - v.begin()), j = 0;
  while (j == i || v[j] != b) j++;
  return {i, j};
}
// The answer is two distinct input elements at the minimum distance.
void check_pair(const vector<Pl>& v, LL want) {
  auto [a, b] = closest_pair(v);
  CHECK(d2(a, b) == want);
  multiset<Pl> ms(all(v));
  CHECK(ms.count(a) >= 1u + (a == b) && ms.count(b) >= 1);
  auto [i, j] = indices(v, a, b);
  CHECK(i != j && v[i] == a && v[j] == b);
}

// The Manhattan variant from the Notes: L-infinity distance of the
// points (x + y, x - y), best itself as the window, from 4e9 + 1.
ll linf(Pl a, Pl b) { return max(abs(a.x - b.x), abs(a.y - b.y)); }
ll manhattan_pair(vector<Pl> v) {
  for (Pl& p : v) p = {p.x + p.y, p.x - p.y};
  sort(all(v), [](Pl a, Pl b) { return a.y < b.y; });
  pair<ll, array<Pl, 2>> best{4'000'000'001, {}};
  set<Pl> s;
  int j = 0;
  for (Pl p : v) {
    Pl d{best.first, 0};
    while (v[j].y < p.y - d.x) s.erase(v[j++]);
    auto lo = s.lower_bound(p - d), hi = s.upper_bound(p + d);
    for (; lo != hi; ++lo) best = min(best, {linf(*lo, p), {*lo, p}});
    s.insert(p);
  }
  auto [a, b] = best.second;
  CHECK(count(all(v), a) && count(all(v), b)); // input points
  return linf(a, b);
}
ll brute_manhattan(const vector<Pl>& v) {
  ll best = LLONG_MAX;
  rep(i, 0, sz(v))
    rep(j, 0, i) {
      best = min(best, abs(v[i].x - v[j].x) + abs(v[i].y - v[j].y));
    }
  return best;
}

int main() {
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(2, 12);
    ll C = it % 4 == 0 ? 1'000'000'000 : rnd(1, 6);
    vector<Pl> v(n);
    for (auto& p : v) p = rnd_pt(C);
    if (it % 5 == 1) // all on a vertical or horizontal line
      for (auto& p : v) (it % 2 ? p.x : p.y) = 0;
    if (it % 7 == 2) v.back() = v[0]; // a duplicate
    check_pair(v, brute(v));
  }
  rep(it, 0, iters(300)) { // larger inputs, tiny grid: many ties
    int n = (int)rnd(20, 200);
    vector<Pl> v(n);
    ll C = rnd(10, 60);
    for (auto& p : v) p = rnd_pt(C);
    check_pair(v, brute(v));
  }
  ll C = 1'000'000'000;
  // Extremes: the only pair spans the whole range.
  check_pair({{-C, -C}, {C, C}}, d2({-C, -C}, {C, C}));
  check_pair({{C, -C}, {-C, C}}, d2({C, -C}, {-C, C}));
  check_pair({{3, 4}, {3, 4}}, 0);
  check_pair({{-C, -C}, {C, C}, {C, -C}}, 4 * (LL)C * C);
  // The Tweak's index recipe on equal points.
  check_pair({{0, 0}, {5, 5}, {9, 1}, {5, 5}}, 0);
  // Manhattan variant: random, then distance 4e9 (the start - 1).
  rep(it, 0, iters(2000)) {
    vector<Pl> v(rnd(2, 12));
    ll c = it % 4 == 0 ? C : rnd(1, 6);
    for (auto& p : v) p = rnd_pt(c);
    CHECK_EQ(manhattan_pair(v), brute_manhattan(v));
  }
  for (ll s : {1, -1}) {
    CHECK_EQ(manhattan_pair({{C, s * C}, {-C, -s * C}}), 4 * C);
    CHECK_EQ(manhattan_pair({{-C, -s * C}, {C, s * C}}), 4 * C);
  }
  // Doubles: works unchanged, at every scale.
  vector<P<double>> vd{{0.5, 0.25}, {3.0, 1.0}, {0.75, 0.0}, {2.5, 2.5}};
  auto [x, y] = closest_pair(vd);
  CHECK(fabs((x - y).dist2() - 0.125) < 1e-12);
  rep(it, 0, iters(2000)) {
    double sc = pow(10.0, (double)rnd(-6, 3));
    vector<P<double>> w(rnd(2, 12));
    for (auto& p : w) p = {rnd_real(0, sc), rnd_real(0, sc)};
    double want = DBL_MAX;
    rep(i, 0, sz(w))
      rep(k, 0, i) want = min(want, (w[i] - w[k]).dist2());
    auto [a, b] = closest_pair(w);
    CHECK(fabs((a - b).dist2() - want) <= 1e-12 * want);
  }

  // Large: 1e5 random points, then 1e5 points on one vertical line
  // and on a narrow grid (worst cases for the window).
  int n = 100'000;
  vector<Pl> v(n);
  for (auto& p : v) p = rnd_pt(C);
  v[n - 1] = v[0] + Pl{1, 0};
  double start = elapsed();
  auto [a, b] = closest_pair(v);
  CHECK_EQ((a - b).dist2(), 1LL);
  rep(i, 0, n) v[i] = {7, 3 * (ll)i - C};
  v[n / 2].y++;
  auto [c, d] = closest_pair(v);
  CHECK_EQ((c - d).dist2(), 4LL);
  rep(i, 0, n) v[i] = {i % 3, i / 3};
  auto [e, f] = closest_pair(v);
  CHECK_EQ((e - f).dist2(), 1LL);
  // Dense doubles in [0, 3]^2: the window must shrink below 1.
  vector<P<double>> w(n);
  rep(i, 0, n) {
    v[i] = {rnd(0, 300'000), rnd(0, 300'000)};
    w[i] = {v[i].x * 1e-5, v[i].y * 1e-5};
  }
  auto [g, h] = closest_pair(w);
  auto [gl, hl] = closest_pair(v);
  double want = (double)(gl - hl).dist2() * 1e-10;
  CHECK(fabs((g - h).dist2() - want) <= 1e-9 * want);
  fprintf(stderr, "  4 x closest pair of 1e5 points: %.2fs\n",
          elapsed() - start);
  puts("OK geometry/closest-pair");
}
