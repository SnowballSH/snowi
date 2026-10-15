// The template's prelude compiles and its macros behave.
#include "../testlib.hpp"

int main() {
  vi v = {3, 1, 2};
  sort(all(v));
  ll sum = 0;
  rep(i, 0, sz(v)) sum += v[i] * (i + 1);
  CHECK_EQ(sum, 14LL);
  CHECK_EQ(sz(vector<pii>(5)), 5);
  puts("OK contest/template");
}
