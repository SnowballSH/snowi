from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/suffix-array")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def lcp_of(s, a, b):
  k = 0
  while max(a, b) + k < len(s) and s[a + k] == s[b + k]:
    k += 1
  return k


def first_rank(S, s, pat, le):
  """As in the Example: first rank whose prefix is not < (<=) pat."""
  m, lo, hi = len(pat), 0, len(S.sa)
  while lo < hi:
    mid = (lo + hi) // 2
    key = s[S.sa[mid] : S.sa[mid] + m]
    if key < pat or (le and key == pat):
      lo = mid + 1
    else:
      hi = mid
  return lo


def check_brute(s, S):
  n = len(s)
  sa = sorted(range(n), key=lambda i: s[i:])
  check_eq(S.sa, sa, s)
  check_eq(S.lcp, [lcp_of(s, sa[i], sa[i + 1]) for i in range(n - 1)])
  check_eq([S.rank[i] for i in sa], list(range(n)))


for _ in range(iters(2000)):
  s = rnd_str(rnd(0, 30), rng.choice(["a", "ab", "abc", "az"]))
  n = len(s)
  S = lib.SuffixArray(s)
  check_brute(s, S)
  check_eq(lib.SuffixArray(s.encode()).sa, S.sa)
  subs = {s[l:r] for l in range(n) for r in range(l + 1, n + 1)}
  check_eq(n * (n + 1) // 2 - sum(S.lcp), len(subs))
  rot = [s[i:] + s[:i] for i in lib.SuffixArray(s + s).sa if i < n]
  check_eq(rot, sorted(s[i:] + s[:i] for i in range(n)))  # Tweak
  for _ in range(5):  # Example: LCP between ranks
    a, b = rnd(0, n), rnd(0, n)
    if a < n and b < n and a != b:
      i, j = sorted((S.rank[a], S.rank[b]))
      check_eq(min(S.lcp[i:j]), lcp_of(s, a, b))
  pat = s[rnd(0, n) :][: rnd(1, 3)] or "b"  # Example: occurrences
  lo, hi = first_rank(S, s, pat, False), first_rank(S, s, pat, True)
  check_eq(sorted(S.sa[lo:hi]), [i for i in range(n) if s.startswith(pat, i)])

# Lists: huge, negative and repeated values; compression keeps order.
for _ in range(iters(500)):
  a = [rng.choice([-(10**30), -1, 0, 7, 10**18]) for _ in range(rnd(0, 25))]
  check_brute(a, lib.SuffixArray(a))
  b = [tuple(rnd_str(rnd(0, 2), "ab")) for _ in range(rnd(0, 12))]
  check_brute(b, lib.SuffixArray(b))  # hashable, not only ints
check_brute([3, 1, 2, 1, 3], lib.SuffixArray([3, 1, 2, 1, 3]))
check_eq(lib.SuffixArray("").sa, [])
check_eq((lib.SuffixArray("z").sa, lib.SuffixArray("z").lcp), ([0], []))

n = 200_000
with Timer("suffix array of a^2e5"):
  S = lib.SuffixArray("a" * n)
check_eq(S.sa, list(range(n - 1, -1, -1)))
check_eq(S.lcp, list(range(1, n)))
for alpha in ("ab", "abcdefghijklmnopqrstuvwxyz"):
  s = rnd_str(n, alpha)
  with Timer(f"suffix array n = 2e5, alphabet {len(alpha)}"):
    S = lib.SuffixArray(s)
  for i in range(n - 1):
    a, b, h = S.sa[i], S.sa[i + 1], S.lcp[i]
    check(b + h < n and (a + h == n or s[a + h] < s[b + h]))
    check(h > 40 or s[a : a + h] == s[b : b + h])
check_eq(sorted(S.sa), list(range(n)))
print("OK strings/suffix-array")
