#include "../testlib.hpp"
#include "../../content/geometry/hull3d.hpp"

typedef P3<ll> Q;
typedef P3<double> Qd;
typedef vector<array<int, 3>> Faces;
typedef __int128 LL;

// fs is the hull surface of p: nonzero outward normals with every
// point behind or on each face's plane, each directed edge once and
// its reverse once (closed and consistently oriented), and Euler's
// formula V - E + F = 2. Returns 6 times the volume.
ll check_hull(const vector<Q>& p, const Faces& fs) {
  CHECK(sz(fs) >= 4);
  map<pii, int> edges;
  set<int> verts;
  ll vol6 = 0;
  for (auto [a, b, c] : fs) {
    CHECK(a != b && b != c && c != a);
    Q q = (p[b] - p[a]).cross(p[c] - p[a]);
    CHECK(q.x || q.y || q.z);
    for (const Q& r : p) CHECK(q.dot(r - p[a]) <= 0);
    edges[{a, b}]++, edges[{b, c}]++, edges[{c, a}]++;
    verts.insert(a), verts.insert(b), verts.insert(c);
    vol6 += p[a].dot(p[b].cross(p[c]));
  }
  for (auto [e, k] : edges) CHECK(k == 1 && edges.count({e.second, e.first}));
  CHECK(sz(verts) - sz(edges) / 2 + sz(fs) == 2);
  CHECK(vol6 > 0);
  return vol6;
}
// 6 * volume of the hull of points in general position: the facets
// are the triples with all other points strictly on one side.
double brute_vol6(const vector<Qd>& p) {
  int n = sz(p);
  double vol6 = 0;
  rep(i, 0, n)
    rep(j, i + 1, n)
      rep(k, j + 1, n) {
        Qd q = (p[j] - p[i]).cross(p[k] - p[i]);
        int pos = 0, neg = 0;
        rep(l, 0, n)
          if (l != i && l != j && l != k) {
            double s = q.dot(p[l] - p[i]);
            pos += s > 0, neg += s < 0;
          }
        double t = p[i].dot(p[j].cross(p[k]));
        if (pos == 0) vol6 += t;
        if (neg == 0) vol6 -= t;
      }
  return vol6;
}
bool coplanar(const vector<Q>& p) {
  for (auto& a : p)
    for (auto& b : p)
      for (auto& c : p) {
        Q q = (b - a).cross(c - a);
        for (auto& d : p)
          if (q.dot(d - a)) return false;
      }
  return true;
}

int main() {
  // Unit cube with face centres, edge midpoints and the centre
  // (coordinates doubled): every face is a flat square.
  vector<Q> cube;
  rep(x, 0, 3)
    rep(y, 0, 3)
      rep(z, 0, 3) cube.push_back({x, y, z});
  shuffle(all(cube), test_rng);
  CHECK_EQ(check_hull(cube, hull3d(cube)), 48);
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(4, it % 2 ? 12 : 30);
    ll c = vector<ll>{1, 2, 5, 100}[it % 4];
    vector<Q> p(n);
    for (Q& q : p) q = {rnd(-c, c), rnd(-c, c), rnd(-c, c)};
    if (it % 5 == 0) p.push_back(p[0]), p.push_back(p[1]);
    if (n <= 12 && coplanar(p)) continue;
    ll vol6 = check_hull(p, hull3d(p));
    if (n > 12) continue;
    // Perturbed copies are in general position; the volume moves by
    // at most about the perturbation times the surface area.
    vector<Qd> pd;
    for (Q& q : p)
      pd.push_back({q.x + rnd_real(-1e-7, 1e-7), q.y + rnd_real(-1e-7, 1e-7),
                    q.z + rnd_real(-1e-7, 1e-7)});
    CHECK(fabs(brute_vol6(pd) - vol6) < 1e-3 * max<ll>(1, c * c));
  }
  // Largest coordinates in the precondition: exact, no overflow.
  rep(it, 0, iters(50)) {
    ll c = 500'000;
    vector<Q> p(100);
    for (Q& q : p)
      q = it % 2 ? Q{rnd(-c, c), rnd(-c, c), rnd(-c, c)}
                 : Q{rnd(0, 1) ? -c : c, rnd(-c, c), rnd(-c, c)};
    check_hull(p, hull3d(p));
  }
  // Delaunay by lifting: lower faces, empty circumcircles, and the
  // triangles tile the 2D hull.
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(3, 25);
    ll c = vector<ll>{2, 4, 20'000}[it % 3];
    vector<Q> p(n);
    for (Q& q : p) {
      q.x = rnd(-c, c), q.y = rnd(-c, c);
      q.z = q.x * q.x + q.y * q.y;
    }
    bool line = true;
    rep(i, 0, n)
      line &= (p[1].x - p[0].x) * (p[i].y - p[0].y) ==
              (p[1].y - p[0].y) * (p[i].x - p[0].x);
    if (line || coplanar(p)) continue; // all on a line or circle
    LL area2 = 0;
    for (auto [a, b, cc] : hull3d(p)) {
      Q q = (p[b] - p[a]).cross(p[cc] - p[a]);
      if (q.z >= 0) continue;
      area2 -= q.z;    // clockwise in the plane
      for (Q& r : p) { // r not inside the circumcircle
        LL m[3][3];
        Q v[3] = {p[a] - r, p[b] - r, p[cc] - r};
        rep(i, 0, 3)
          m[i][0] = v[i].x, m[i][1] = v[i].y,
          m[i][2] = v[i].x * v[i].x + v[i].y * v[i].y;
        LL det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        CHECK(det >= 0); // clockwise triangle: inside iff det < 0
      }
    }
    vector<pair<ll, ll>> h; // 2D hull area by monotone chain
    for (Q& q : p) h.push_back({q.x, q.y});
    sort(all(h)), h.erase(unique(all(h)), h.end());
    auto cr = [](pair<ll, ll> o, pair<ll, ll> a, pair<ll, ll> b) {
      return (LL)(a.first - o.first) * (b.second - o.second) -
             (LL)(a.second - o.second) * (b.first - o.first);
    };
    vector<pair<ll, ll>> st;
    rep(pass, 0, 2) {
      size_t base = st.size();
      for (auto q : h) {
        while (st.size() >= base + 2 &&
               cr(st[st.size() - 2], st.back(), q) <= 0)
          st.pop_back();
        st.push_back(q);
      }
      st.pop_back();
      reverse(all(h));
    }
    LL want = 0;
    rep(i, 0, sz(st)) want += cr({0, 0}, st[i], st[(i + 1) % sz(st)]);
    CHECK(area2 == want);
  }
  int n = 2000;
  vector<Q> p(n);
  for (Q& q : p) { // near a sphere: most points are on the hull
    double t = rnd_real(0, 2 * M_PI), z = rnd_real(-1, 1);
    double r = sqrt(1 - z * z);
    q = {ll(1e5 * r * cos(t)), ll(1e5 * r * sin(t)), ll(1e5 * z)};
  }
  double start = elapsed();
  Faces fs = hull3d(p);
  fprintf(stderr, "  hull3d of 2000 points near a sphere: %.2fs (%d)\n",
          elapsed() - start, sz(fs));
  check_hull(p, fs);
  puts("OK geometry/hull3d");
}
