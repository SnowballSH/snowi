#include "../testlib.hpp"
#include "../../content/offline/cdq.hpp"

typedef array<int, 3> P3;

vi brute(const vector<P3>& p) {
  int n = sz(p);
  vi cnt(n);
  rep(i, 0, n)
    rep(j, 0, n)
      if (j != i) {
        bool le = true;
        rep(d, 0, 3) le &= p[j][d] <= p[i][d];
        cnt[i] += le;
      }
  return cnt;
}

vector<P3> random_points(int n, int lo, int hi) {
  vector<P3> p(n);
  for (auto& x : p)
    for (int& v : x) v = (int)rnd(lo, hi);
  return p;
}

// The Example's CDQ on time: op t inserts (X[t], Y[t]) or, if
// ask[t], counts inserted points <= (X[t], Y[t]).
vector<ll> cdq_time(const vi& X, const vi& Y, const vi& ask) {
  int q = sz(X);
  vi idx(q), xs(X);
  iota(all(idx), 0);
  sort(all(xs));
  vi cx(q);
  rep(t, 0, q) cx[t] = int(lower_bound(all(xs), X[t]) - xs.begin());
  vector<ll> res(q);
  Fenwick fw(q);
  auto byy = [&](int i, int j) { return Y[i] < Y[j]; };
  auto rec = [&](auto&& self, int l, int r) -> void {
    if (r - l < 2) return;
    int m = (l + r) / 2, k = l;
    self(self, l, m), self(self, m, r);
    rep(i, m, r) {
      int t = idx[i];
      for (; k < m && Y[idx[k]] <= Y[t]; k++)
        if (!ask[idx[k]]) fw.add(cx[idx[k]], 1);
      if (ask[t]) res[t] += fw.query(cx[t] + 1);
    }
    rep(j, l, k)
      if (!ask[idx[j]]) fw.add(cx[idx[j]], -1);
    inplace_merge(idx.begin() + l, idx.begin() + m, idx.begin() + r, byy);
  };
  rec(rec, 0, q);
  return res;
}

void stress_time(int rounds, int maxq, int V) {
  rep(it, 0, iters(rounds)) {
    int q = (int)rnd(0, maxq);
    vi X(q), Y(q), ask(q);
    rep(t, 0, q)
      X[t] = (int)rnd(0, V), Y[t] = (int)rnd(0, V), ask[t] = rnd(0, 1);
    vector<ll> res = cdq_time(X, Y, ask);
    rep(t, 0, q)
      if (ask[t]) {
        ll want = 0;
        rep(s, 0, t) want += !ask[s] && X[s] <= X[t] && Y[s] <= Y[t];
        CHECK_EQ(res[t], want);
      }
  }
}

int main() {
  // Fixed cases: the packet's example, copies, n = 0 and 1.
  CHECK(dominance({{1, 3, 2}, {2, 1, 1}, {3, 2, 3}, {4, 4, 2}}) ==
        vi({0, 0, 1, 2}));
  CHECK(dominance({{1, 1, 1}, {2, 2, 2}, {1, 1, 1}, {0, 5, 0}}) ==
        vi({1, 2, 1, 0}));
  CHECK(dominance({}).empty());
  CHECK(dominance({{7, -3, 2}}) == vi({0}));
  CHECK(dominance(vector<P3>(5, {4, 4, 4})) == vi(5, 4));
  int M = INT_MAX, m = INT_MIN;
  CHECK(dominance({{M, M, M}, {m, m, m}, {M, m, M}, {m, M, m}}) ==
        vi({3, 0, 1, 1}));
  // Small ranges give many ties in every coordinate.
  for (auto [rounds, maxn, V] :
       {array<int, 3>{3000, 12, 2}, {3000, 25, 4}, {500, 60, 50}}) {
    rep(it, 0, iters(rounds)) {
      auto p = random_points((int)rnd(0, maxn), -V, V);
      CHECK(dominance(p) == brute(p));
    }
  }
  stress_time(2000, 30, 5);
  stress_time(300, 120, 1000);
  // Large: n = 1e5, checked on a sample of points.
  elapsed();
  int n = 100'000;
  auto p = random_points(n, 0, 1000);
  vi cnt = dominance(p);
  fprintf(stderr, "  n = 1e5: %.2fs\n", elapsed());
  rep(s, 0, 200) {
    int i = (int)rnd(0, n - 1), want = 0;
    rep(j, 0, n)
      want += j != i && p[j][0] <= p[i][0] && p[j][1] <= p[i][1] &&
              p[j][2] <= p[i][2];
    CHECK_EQ(cnt[i], want);
  }
  puts("OK offline/cdq");
}
