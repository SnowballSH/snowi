#include "../testlib.hpp"
#include "../../content/numtheory/modpow.hpp"

int main() {
  rep(it, 0, iters(2000)) {
    ll b = rnd(0, MOD - 1), e = rnd(0, 60), want = 1;
    rep(k, 0, e) want = want * b % MOD;
    CHECK_EQ(modpow(b, e), want);
    ll e1 = rnd(0, (ll)1e18), e2 = rnd(0, (ll)1e18);
    CHECK_EQ(modpow(b, e1) * modpow(b, e2) % MOD, modpow(b, e1 + e2));
    if (b) CHECK_EQ(b * modinv(b) % MOD, 1LL);
    ll big = rnd(0, LLONG_MAX); // unreduced base
    CHECK_EQ(modpow(big, e), modpow(big % MOD, e));
  }
  CHECK_EQ(modpow(0, 0), 1LL);
  CHECK_EQ(modpow(MOD + 2, 1), 2LL);
  puts("OK numtheory/modpow");
}
