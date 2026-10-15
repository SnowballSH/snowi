/**
 * Title: Contest template
 * Description: Base template (from kactl). Every C++ snippet assumes
 *   these macros and typedefs. `cin.exceptions` turns a failed read
 *   into a crash instead of silently wrong answers. Compile and
 *   run commands: Contest setup.
 * Usage:
 *   int n; cin >> n; vi a(n);
 *   rep(i, 0, n) cin >> a[i];
 *   sort(all(a)); ll s = accumulate(all(a), 0LL);
 * Source: kactl (CC0).
 * Status: compiled and run by tests/contest/template.test.cpp
 */
#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for (int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

int main() {
  cin.tie(0)->sync_with_stdio(0);
  cin.exceptions(cin.failbit);
}
