/**
 * Title: Gaussian elimination mod p
 * Description: `gauss(a, m)`: reduced row echelon form in place,
 *   pivots from the first `m` columns (later columns ride along);
 *   returns the pivot columns and, if `a` is square and `m = n`, the
 *   determinant.
 * Precondition: `MOD` prime, entries in $[0, "MOD")$. `a` is
 *   `n x m` with any `n`, `m`; size `x` to `m` before
 *   `solve_linear`.
 * Time: $O(n m min(n, m))$. For $n = m = 500$ `solve_linear` takes
 *   about 0.2 s and `mat_inv` 0.6 s.
 * Usage:
 *   vector<ll> x(m);
 *   int r = solve_linear(a, b, x); // -1 if none, else rank
 *   ll d = det(a);
 *   auto inv = mat_inv(a); // {} if singular
 *   // Last: gauss works in place, a becomes its RREF.
 *   vi piv = gauss(a, m).first; // rank sz(piv)
 * Tweak: Over $FF_2$ with many columns use rows of
 *   `bitset<2048>`: drop the scaling and eliminate with
 *   `a[i] ^= a[r]`, so a row update costs $m\/64$ word operations.
 * Example:
 *   // Kernel basis of an n x m matrix: one vector per free
 *   // column f, with x_f = 1 and x_c = -a[i][f] for pivot c.
 *   vector<vector<ll>> ker;
 *   vi piv = gauss(a, m).first, is_piv(m);
 *   for (int c : piv) is_piv[c] = 1;
 *   rep(f, 0, m) if (!is_piv[f]) {
 *     vector<ll> v(m);
 *     v[f] = 1;
 *     rep(i, 0, sz(piv)) v[piv[i]] = (MOD - a[i][f]) % MOD;
 *     ker.push_back(v);
 *   }
 * Notes: When consistent, $A x = b$ has $"MOD"^(m - r)$ solutions:
 *   the one returned (free variables 0) plus any kernel vector. Rank
 *   and solvability depend on the field: a matrix can be singular
 *   mod $p$ but not over the rationals. Cues: "count assignments
 *   satisfying linear conditions mod p" ($p^(m - r)$), "determinant
 *   / number of spanning trees mod p" (matrix-tree), random-weight
 *   rank tests (Schwartz-Zippel). With a composite modulus, pivots
 *   may have no inverse: for a determinant, clear each column with
 *   Euclid-style row subtractions in $O(n^3 log "MOD")$.
 * Source: Standard Gauss-Jordan; kactl MatrixInverse-mod (CC0).
 * Status: stress-tested against permutation determinants and minor
 *   ranks
 */
#pragma once
#include "../numtheory/modpow.hpp"

pair<vi, ll> gauss(vector<vector<ll>>& a, int m) {
  int n = sz(a), r = 0;
  vi piv;
  ll det = 1;
  for (int c = 0; c < m && r < n; c++) {
    int p = r;
    while (p < n && !a[p][c]) p++;
    if (p == n) continue;
    if (p != r) swap(a[p], a[r]), det = MOD - det;
    det = det * a[r][c] % MOD;
    ll inv = modinv(a[r][c]);
    for (ll& v : a[r]) v = v * inv % MOD;
    rep(i, 0, n)
      if (i != r && a[i][c]) {
        ll f = MOD - a[i][c];
        rep(j, c, sz(a[i]))
          a[i][j] = (a[i][j] + f * a[r][j]) % MOD;
      }
    piv.push_back(c), r++;
  }
  return {piv, r == n ? det : 0};
}

ll det(vector<vector<ll>> a) { return gauss(a, sz(a)).second; }

int solve_linear(vector<vector<ll>> a, const vector<ll>& b,
                 vector<ll>& x) {
  int n = sz(a), m = sz(x);
  rep(i, 0, n) a[i].push_back(b[i]);
  vi piv = gauss(a, m).first;
  int r = sz(piv);
  rep(i, r, n)
    if (a[i][m]) return -1;
  x.assign(m, 0);
  rep(i, 0, r) x[piv[i]] = a[i][m];
  return r;
}

vector<vector<ll>> mat_inv(vector<vector<ll>> a) {
  int n = sz(a);
  rep(i, 0, n) a[i].resize(2 * n), a[i][n + i] = 1;
  if (sz(gauss(a, n).first) < n) return {};
  for (auto& row : a) row.erase(row.begin(), row.begin() + n);
  return a;
}
