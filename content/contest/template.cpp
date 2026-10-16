/**
 * Title: Contest template
 * Description: `cin.exceptions`: a failed read crashes instead of
 *   silently giving wrong answers.
 * Usage:
 *   int n; cin >> n; vi a(n); for (int& x : a) cin >> x;
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
