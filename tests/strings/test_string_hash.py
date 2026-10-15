from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/string-hash")
StrHash = lib.StrHash
M = (1 << 61) - 1


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def evaluate(vals, l, r):  # independent: big-int Horner
  h = 0
  for i in range(l, r):
    h = (h * StrHash.B + vals[i]) % M
  return h


check_eq(StrHash.M, M)
check(0 <= StrHash.B < M)
for _ in range(iters(20000)):  # mul against big-int arithmetic
  a = rng.choice([rnd(0, M - 1), rnd(M - 5, M - 1), rnd(0, 5)])
  b = rng.choice([rnd(0, M - 1), rnd(M - 5, M - 1), rnd(0, 5)])
  check_eq(StrHash.mul(a, b), a * b % M, (a, b))

for _ in range(iters(1000)):
  s = rnd_str(rnd(0, 14), rng.choice(["a", "ab", "abc"]))
  n = len(s)
  hs, hb = StrHash(s), StrHash(s.encode())
  vals = [ord(c) for c in s]
  for l in range(n + 1):
    for r in range(l, n + 1):
      h = hs.get(l, r)
      check_eq(h, evaluate(vals, l, r))
      check_eq(hb.get(l, r), h)
      for l2 in range(n - (r - l) + 1):
        same = s[l:r] == s[l2 : l2 + r - l]
        check_eq(hs.get(l2, l2 + r - l) == h, same)
  # Tweak: concatenation, with pw of the string holding y
  t = rnd_str(rnd(0, 8), "abc")
  ht, hst = StrHash(t), StrHash(s + t)
  for l in range(n + 1):
    for k in range(len(t) + 1):
      x, y = hs.get(l, n), ht.get(0, k)
      check_eq((ht.mul(x, ht.pw[k]) + y) % ht.M, hst.get(l, n + k))

for _ in range(iters(1000)):  # values near M, zeros, sums equal to M
  n = rnd(1, 12)
  a = [rnd(0, 2) if rnd(0, 3) else rnd(M - 3, M - 1) for _ in range(n)]
  h = 0
  for i in range(n):
    if h and rnd(0, 2) == 0:
      a[i] = M - h * StrHash.B % M
    h = (h * StrHash.B + a[i]) % M
  ha = StrHash(a)
  for l in range(n + 1):
    for r in range(l, n + 1):
      check_eq(ha.get(l, r), evaluate(a, l, r))

check_eq(StrHash("").get(0, 0), 0)
tm = "".join("ab"[bin(i).count("1") & 1] for i in range(2048))
co = tm.translate(str.maketrans("ab", "ba"))
check(StrHash(tm).get(0, 2048) != StrHash(co).get(0, 2048))

# Speed: n = 1e6 build, 1e6 gets; checked against slicing on a
# periodic string with a few changes, so many pairs are equal.
n = 10**6
s = list("abcab" * (n // 5))
for _ in range(20):
  s[rnd(0, n - 1)] = "c"
s = "".join(s)
queries = []
for _ in range(n // 2):
  ln = rnd(1, 50)
  a = rnd(0, n - ln)
  b = min(n - ln, rnd(0, (n - ln) // 5) * 5 + a % 5 * rnd(0, 1))
  queries.append((a, b, ln))
with Timer("build n = 1e6 + 1e6 get"):
  hs = StrHash(s)
  got = [hs.get(a, a + ln) == hs.get(b, b + ln) for a, b, ln in queries]
want = [s[a : a + ln] == s[b : b + ln] for a, b, ln in queries]
check(got == want)
check(sum(want) > n // 20)
print("OK strings/string-hash")
