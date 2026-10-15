from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/prefix-function")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def brute_pi(s):
  return [
    max(k for k in range(i + 1) if s[:k] == s[i + 1 - k : i + 1])
    for i in range(len(s))
  ]


def brute_match(s, pat):
  m = len(pat)
  return [i for i in range(len(s) - m + 1) if s[i : i + m] == pat]


def brute_period(s):
  n = len(s)
  return next(p for p in range(1, n + 1) if s[p:] == s[: n - p])


for _ in range(iters(2000)):
  alpha = rng.choice(["a", "ab", "abc"])
  s = rnd_str(rnd(0, 24), alpha)
  pi = lib.prefix_function(s)
  check_eq(pi, brute_pi(s), s)
  if s:
    check_eq(len(s) - pi[-1], brute_period(s), s)
  if s:  # Notes: occurrences of every prefix of s in a text t
    t, cnt, k = rnd_str(rnd(0, 30), alpha), [0] * (len(s) + 1), 0
    for c in t:
      while k and s[k] != c:
        k = pi[k - 1]
      if s[k] == c:
        k += 1
      cnt[k] += 1
      if k == len(s):
        k = pi[k - 1]
    for L in range(len(s), 0, -1):
      cnt[pi[L - 1]] += cnt[L]
    for L in range(1, len(s) + 1):
      check_eq(cnt[L], len(brute_match(t, s[:L])))
  pat = rnd_str(rnd(1, 5), alpha)
  if s and rnd(0, 3) == 0:
    l = rnd(0, len(s) - 1)
    pat = s[l : rnd(l + 1, len(s))]
  check_eq(lib.kmp_match(s, pat), brute_match(s, pat), (s, pat))
  # bytes and lists of big ints behave like the str
  b, bp = s.encode(), pat.encode()
  check_eq(lib.kmp_match(b, bp), brute_match(s, pat))
  big = [10**30 * ord(c) for c in s]
  check_eq(lib.prefix_function(big), pi)
  check_eq(
    lib.kmp_match(big, [10**30 * ord(c) for c in pat]), brute_match(s, pat)
  )

check_eq(lib.kmp_match("ab", "abc"), [])
check_eq(lib.kmp_match("", "a"), [])
check_eq(lib.prefix_function(""), [])

n = 10**6
s = rnd_str(n, "ab")
pat = "abbabaab"
with Timer("prefix function n = 1e6 + match in 1e6"):
  pi = lib.prefix_function(s)
  occ = lib.kmp_match(s, pat)
want, p = [], s.find(pat)
while p >= 0:
  want.append(p)
  p = s.find(pat, p + 1)
check_eq(occ, want)
check(all(pi[i] <= pi[i - 1] + 1 for i in range(1, n)))
with Timer("prefix function of a^1e6, bytes"):
  pi = lib.prefix_function(b"a" * n)
check(pi == list(range(n)))
print("OK strings/prefix-function")
