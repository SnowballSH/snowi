from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/manacher")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def is_pal(s, l, r):
  t = s[l:r]
  return t == t[::-1]


for _ in range(iters(2000)):
  s = rnd_str(rnd(0, 24), rng.choice(["a", "ab", "abc"]))
  n = len(s)
  plen = lib.manacher(s)
  check_eq(len(plen), max(2 * n - 1, 0))
  check_eq(lib.manacher(s.encode()), plen)
  check_eq(lib.manacher([10**30 * ord(c) for c in s]), plen)
  want = [0] * len(plen)
  for l in range(n):
    for r in range(l + 1, n + 1):
      pal = is_pal(s, l, r)
      check_eq(plen[l + r - 1] >= r - l, pal, (s, l, r))
      if pal:
        want[l + r - 1] = max(want[l + r - 1], r - l)
  check_eq(plen, want, s)
  pals = {s[l:r] for l in range(n) for r in range(l + 1, n + 1)}
  check(sum(t == t[::-1] for t in pals) <= n)  # Notes: at most n
  if n:  # Usage: longest palindromic substring
    c = max(range(len(plen)), key=plen.__getitem__)
    l = (c + 1 - plen[c]) // 2
    check(is_pal(s, l, l + plen[c]))
    check_eq(plen[c], max(want))
  # Notes: palindromes starting and ending at each position
  starts, ends = [0] * (n + 1), [0] * (n + 1)
  for c, p in enumerate(plen):
    starts[(c + 1 - p) // 2] += 1
    starts[c // 2 + 1] -= 1
    ends[(c + 1) // 2] += 1
    ends[(c + p - 1) // 2 + 1] -= 1
  for i in range(n):
    starts[i + 1] += starts[i]
    ends[i + 1] += ends[i]
    check_eq(starts[i], sum(is_pal(s, i, r) for r in range(i + 1, n + 1)))
    check_eq(ends[i], sum(is_pal(s, l, i + 1) for l in range(i + 1)))

check_eq(lib.manacher(""), [])
check_eq(lib.manacher("a"), [1])

n = 10**6
with Timer("manacher of a^1e6"):
  plen = lib.manacher("a" * n)
check(plen == [c + 1 if c < n else 2 * n - 1 - c for c in range(2 * n - 1)])
s = rnd_str(n, "ab")
with Timer("manacher n = 1e6, random str"):
  plen = lib.manacher(s)
for c in range(2 * n - 1):  # maximal: a palindrome that cannot grow
  p = plen[c]
  l = (c + 1 - p) // 2
  r = l + p
  check(l + r - 1 == c and (l == 0 or r == n or s[l - 1] != s[r]))
  check(p > 40 or is_pal(s, l, r))
print("OK strings/manacher")
