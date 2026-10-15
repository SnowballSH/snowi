// Shared helpers for C++ stress tests. Include before the snippet.
#pragma once
// The template's main() is compiled under another name and never run.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main snowi_template_main
#include "../content/contest/template.cpp"
#undef main
#pragma GCC diagnostic pop

inline mt19937_64 test_rng = [] {
  const char* seed = getenv("SNOWI_SEED");
  return mt19937_64(seed ? strtoull(seed, nullptr, 10) : 20261015);
}();

// Uniform integer in [lo, hi].
inline ll rnd(ll lo, ll hi) {
  return uniform_int_distribution<ll>(lo, hi)(test_rng);
}
inline double rnd_real(double lo, double hi) {
  return uniform_real_distribution<double>(lo, hi)(test_rng);
}
inline vector<ll> rnd_vec(int n, ll lo, ll hi) {
  vector<ll> v(n);
  for (auto& x : v) x = rnd(lo, hi);
  return v;
}
// Scale iteration counts: SNOWI_SCALE=0.1 for a quick run.
inline int iters(int n) {
  const char* s = getenv("SNOWI_SCALE");
  return max(1, (int)(n * (s ? atof(s) : 1.0)));
}

#define CHECK(cond)                                              \
  do {                                                           \
    if (!(cond)) {                                               \
      fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__,     \
              __LINE__, #cond);                                  \
      exit(1);                                                   \
    }                                                            \
  } while (0)

#define CHECK_EQ(got, want)                                      \
  do {                                                           \
    auto&& got_ = (got);                                         \
    auto&& want_ = (want);                                       \
    if (!(got_ == want_)) {                                      \
      cerr << __FILE__ << ":" << __LINE__ << ": CHECK_EQ failed: " \
           << #got << " != " << #want << "\n";                   \
      if constexpr (requires { cerr << got_ << want_; })         \
        cerr << "  got " << got_ << "\n  want " << want_ << "\n"; \
      exit(1);                                                   \
    }                                                            \
  } while (0)

// Wall-clock seconds since the first call; for speed sanity checks.
inline double elapsed() {
  static auto start = chrono::steady_clock::now();
  return chrono::duration<double>(chrono::steady_clock::now() - start)
      .count();
}
