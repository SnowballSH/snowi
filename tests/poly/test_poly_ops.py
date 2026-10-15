from testlib import Timer, check_eq, iters, load, rnd, rnd_list

lib = load("poly/poly-ops")
MOD = lib.MOD


# Independent O(n^2) series arithmetic modulo x^n.
def trunc(a, n):
  return a[:n] + [0] * (n - len(a))


def mul_trunc(a, b, n):
  c = [0] * n
  for i, x in enumerate(a[:n]):
    for j, y in enumerate(b[: n - i]):
      c[i + j] = (c[i + j] + x * y) % MOD
  return c


def brute_inv(a, n):
  a, b, i0 = trunc(a, n), [0] * n, pow(a[0], MOD - 2, MOD)
  for i in range(n):
    s = int(i == 0) - sum(a[j] * b[i - j] for j in range(1, i + 1))
    b[i] = s * i0 % MOD
  return b


def brute_log(a, n):  # i L_i = i a_i - sum_(0 < j < i) j L_j a_(i - j)
  a, L = trunc(a, n), [0] * n
  for i in range(1, n):
    s = a[i] * i - sum(j * L[j] * a[i - j] for j in range(1, i))
    L[i] = s * pow(i, MOD - 2, MOD) % MOD
  return L


def brute_exp(a, n):  # i e_i = sum_(0 < j <= i) j a_j e_(i - j)
  a, e = trunc(a, n), [1] * min(n, 1) + [0] * (n - 1)
  for i in range(1, n):
    s = sum(j * a[j] * e[i - j] for j in range(1, i + 1))
    e[i] = s * pow(i, MOD - 2, MOD) % MOD
  return e


def brute_pow(a, k, n):
  r = trunc([1], n)
  while k:
    if k & 1:
      r = mul_trunc(r, a, n)
    a, k = mul_trunc(a, a, n), k >> 1
  return r


def brute_divmod(a, b):  # long division
  a, n, m = a[:], len(a), len(b)
  q, ib = [0] * max(0, n - m + 1), pow(b[-1], MOD - 2, MOD)
  for i in range(n - m, -1, -1):
    q[i] = a[i + m - 1] * ib % MOD
    for j in range(m):
      a[i + j] = (a[i + j] - q[i] * b[j]) % MOD
  return q, trunc(a, m - 1)


def online(h, g, step):  # cdq as printed in the techniques note
  n, conv = len(h), lib.conv
  f, acc = [0] * n, [0] * n

  def finish(t, a):
    return step(h[t], a) % MOD

  def cdq(l, r):  # f, g, acc: lists of length n
    if r - l == 1:
      f[l] = finish(l, acc[l])
      return
    m = (l + r) // 2
    cdq(l, m)
    c = conv(f[l:m], g[: r - l])
    for t in range(m, r):
      acc[t] = (acc[t] + c[t - l]) % MOD
    cdq(m, r)

  if n:
    cdq(0, n)
  return f


def online_brute(h, g, step):  # f_t = step(h_t, sum_(0<i<=t) g_i f_(t-i))
  f = []
  for t in range(len(h)):
    a = sum(g[i] * f[t - i] for i in range(1, t + 1))
    f.append(step(h[t], a % MOD) % MOD)
  return f


def online_closed(h, g):  # linear step: f = h / (1 - g) with g_0 = 0
  n = len(h)
  d = [1] + [-x % MOD for x in g[1:n]]
  return lib.conv(h, lib.poly_inv(d, n))[:n]


def linear(x, a):
  return x + a


def squared(x, a):
  return x + a * a


def random_poly(length, it):
  if it % 9 == 0:
    return [MOD - 1] * length
  return rnd_list(length, 0, MOD - 1 if it % 4 else 2)


for it in range(iters(400)):
  n = rnd(0, 30 if it % 5 else 100)
  a = random_poly(rnd(1, n + 3), it)
  a[0] = rnd(1, MOD - 1)
  check_eq(lib.poly_inv(a, n), brute_inv(a, n))
  a[0] = 1
  L = lib.poly_log(a, n)
  check_eq(L, brute_log(a, n))
  check_eq(lib.poly_exp(L, n), trunc(a, n))
  s = lib.poly_pow(a, (MOD + 1) // 2, n)  # square root, a_0 = 1
  check_eq(mul_trunc(s, s, n), trunc(a, n))
  # Tweak for a_0 != 1: r^2 x^(2w) a has the root r x^w sqrt(a).
  w, r = rnd(0, 3), rnd(1, MOD - 1)
  sq = [0] * (2 * w) + [x * r * r % MOD for x in a]
  root = [0] * w + [x * r % MOD for x in s]
  check_eq(mul_trunc(root, root, n), trunc(sq, n))
  a[0] = 0
  check_eq(lib.poly_exp(a, n), brute_exp(a, n))
  m, zeros = min(n, 30), rnd(0, 4)
  p = random_poly(rnd(1, m + 3), it)
  p[: min(zeros, len(p))] = [0] * min(zeros, len(p))
  if it % 11 == 0:
    p = [0] * len(p)
  k = rnd(0, 6) if it % 3 else rnd(0, 10**18)
  if it % 7 == 0 and zeros and m:
    k = (m - 1) // zeros + rnd(0, 1)  # v k = n - 1 or >= n
  check_eq(lib.poly_pow(p, k, m), brute_pow(p, k, m), f"k={k}")
  x = random_poly(rnd(0, 50), it)
  y = random_poly(rnd(1, 25 if it % 4 else 60), it)
  if y[-1] == 0:
    y[-1] = 1
  check_eq(lib.poly_divmod(x, y), brute_divmod(x, y))

for n in range(3):
  check_eq(lib.poly_inv([5], n), brute_inv([5], n))
  check_eq(lib.poly_log([1], n), [0] * n)
  check_eq(lib.poly_exp([0], n), trunc([1], n))
  check_eq(lib.poly_pow([0], 0, n), trunc([1], n))
  check_eq(lib.poly_pow([0, 3], 5, n), [0] * n)

# Online convolution (techniques note), g_0 != 0 included.
h0, g0 = [1, 0, 0, 0, 0, 0], [5, 1, 2, 3, 4, 5]
check_eq(online(h0, g0, linear), [1, 1, 3, 8, 21, 55])
check_eq(online_closed(h0, g0), [1, 1, 3, 8, 21, 55])
for it in range(iters(150)):
  n = rnd(0, 30 if it % 10 else 200)
  h, g = random_poly(n, it), random_poly(n, it + 1)
  if n and it % 3 == 0:
    g[0] = 1
  want = online_brute(h, g, linear)
  check_eq(online(h, g, linear), want)
  check_eq(online_closed(h, g), want)
  check_eq(online(h, g, squared), online_brute(h, g, squared))

# The header's partition Example at size, against Euler's pentagonal
# recurrence. The timers show building the log is cheap next to exp.
n = 1 << 16
with Timer("partition log 2^16"):
  inv = [1] * n  # inv[j] = 1 / j
  for i in range(2, n):
    inv[i] = (MOD - MOD // i) * inv[MOD % i] % MOD
  L = [0] * n
  for k in range(1, n):
    for j in range(1, (n - 1) // k + 1):
      L[j * k] += inv[j]
with Timer("partition exp 2^16"):
  p = lib.poly_exp([x % MOD for x in L], n)
want = [1] + [0] * (n - 1)
for i in range(1, n):
  s, k = 0, 1
  while k * (3 * k - 1) // 2 <= i:
    g = k * (3 * k - 1) // 2
    t = want[i - g] + (want[i - g - k] if g + k <= i else 0)
    s += t if k % 2 else -t
    k += 1
  want[i] = s % MOD
check_eq(p, want)

n = 1 << 17
a = rnd_list(n, 0, MOD - 1)
a[0] = 1
with Timer("inverse 2^17"):
  b = lib.poly_inv(a, n)
check_eq(lib.conv(a, b)[:n], trunc([1], n))
with Timer("log 2^17"):
  L = lib.poly_log(a, n)
with Timer("exp 2^17"):
  check_eq(lib.poly_exp(L, n), a)
n = 1 << 14
a = a[:n]
sq = lib.conv(a, a)[:n]
check_eq(lib.poly_pow(sq, (MOD + 1) // 2, n), a)
b = a[: n // 2]
b[-1] = b[-1] or 1
with Timer("divmod 2^14 / 2^13"):
  q, r = lib.poly_divmod(a, b)
back = lib.conv(q, b)
check_eq([(x + y) % MOD for x, y in zip(back, r + [0] * n)], a)
print("OK poly/poly-ops")
