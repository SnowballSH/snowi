from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/z-function")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def lcp(a, b):
  k = 0
  while k < len(a) and k < len(b) and a[k] == b[k]:
    k += 1
  return k


def brute_z(s):
  return [lcp(s, s[i:]) for i in range(len(s))]


def brute_pi(s):
  return [
    max(k for k in range(i + 1) if s[:k] == s[i + 1 - k : i + 1])
    for i in range(len(s))
  ]


for _ in range(iters(2000)):
  alpha = rng.choice(["a", "ab", "abc"])
  s = rnd_str(rnd(0, 24), alpha)
  n = len(s)
  z = lib.z_function(s)
  check_eq(z, brute_z(s), s)
  check_eq(lib.z_function(s.encode()), z)
  check_eq(lib.z_function([10**30 * ord(c) for c in s]), z)
  # Usage: matching through pat + "#" + s, and with lists and -1
  pat = rnd_str(rnd(1, 4), "ab")
  m = len(pat)
  zt = lib.z_function(pat + "#" + s)
  check_eq(
    [zt[m + 1 + i] for i in range(n)], [lcp(pat, s[i:]) for i in range(n)]
  )
  zl = lib.z_function(list(pat.encode()) + [-1] + list(s.encode()))
  want = [i for i in range(n) if s[i : i + m] == pat]
  check_eq([i for i in range(n) if zl[m + 1 + i] == m], want)
  # Notes and Example: periods (p = n has no z[p]) and pi from z
  for p in range(1, n + 1):
    check_eq(p == n or z[p] == n - p, s[p:] == s[: n - p])
  pi = [0] * n
  for i in range(1, n):
    if z[i]:
      pi[i + z[i] - 1] = max(pi[i + z[i] - 1], z[i])
  for i in range(n - 2, -1, -1):
    pi[i] = max(pi[i], pi[i + 1] - 1)
  check_eq(pi, brute_pi(s))

check_eq(lib.z_function(""), [])
check_eq(lib.z_function("a"), [1])
z = lib.z_function("abaab")  # Example: all periods, p = n included
check_eq([p for p in range(1, 6) if p == 5 or z[p] == 5 - p], [3, 5])

n = 10**6
s = rnd_str(n, "ab")
with Timer("z-function n = 1e6, str"):
  z = lib.z_function(s)
for i in range(1, n):  # z[i] is a common prefix that cannot grow
  k = z[i]
  check(s[:k] == s[i : i + k] and (i + k == n or s[k] != s[i + k]))
with Timer("z-function of a^1e6, bytes"):
  z = lib.z_function(b"a" * n)
check(z == list(range(n, 0, -1)))
print("OK strings/z-function")
