"""
Title: Binomial coefficients
Description: Factorials modulo `MOD` up to `n` give $binom(n, k)$
  (0 outside $0 <= k <= n$), $n!\/(n - k)!$, multinomials, Catalan
  numbers and $1\/i$ in $O(1)$. `Lucas` gives $binom(n, k)$ modulo
  a small prime `p` for huge `n` and `k`.
Precondition: `MOD` prime and $> n$. Arguments at most `n` (in
  `multinom` their sum): `catalan(k)` needs $2k <= n$, `inv(i)`
  needs $i >= 1$. `Lucas`: `p` prime up to about $10^6$;
  $n, k >= 0$.
Time: $O(n)$ to build, $O(1)$ per query; $n = 2 dot 10^6$ builds
  in about 0.15 s, $10^7$ in 0.6 s. `Lucas`: $O(p)$ to build,
  $O(log_p n)$ per query.
Usage:
  bc = Binom(2 * n)
  bc.C(5, 2), bc.perm(5, 2), bc.catalan(3)  # 10 20 5
  bc.multinom([2, 1, 1]), bc.inv(4)  # 12, 1/4
  Lucas(10007).C(10**18, 12345)
Tweak: In hot loops alias the lists, `f, fi = bc.fact, bc.ifact`,
  and inline `f[n] * fi[k] % MOD * fi[n - k] % MOD`.
Uses: numtheory/modpow
Example:
  >>> bc = Binom(10)
  >>> bc.C(5, 2), bc.C(5, 6), bc.C(-1, 0), bc.perm(5, 2)
  (10, 0, 0, 20)
  >>> [bc.catalan(k) for k in range(6)]
  [1, 1, 2, 5, 14, 42]
  >>> bc.multinom([2, 1, 1])  # arrangements of "aabc"
  12
  >>> bc.inv(3) * 3 % MOD
  1
  >>> n, k = 4, 3  # 4 equal balls into 3 boxes
  >>> bc.C(n + k - 1, k - 1)
  15
  >>> Lucas(3).C(14, 4)  # C(14, 4) = 1001 = 3 * 333 + 2
  2
  >>> Lucas(2).C(12, 4), Lucas(2).C(12, 3)  # 4 is a submask
  (1, 0)
Notes: Size the table for the largest argument: Catalan and
  stars and bars read up to $2n$ or $n + k$. For huge `n` and small
  `k` multiply $n (n - 1) dots (n - k + 1)$ by `ifact[k]` in $O(k)$.
  `C` returns 0 out of range, which handles unreachable cells and
  negative counts in inclusion-exclusion sums, but stars and bars
  with $n = k = 0$ needs the answer 1 by hand. A modulus $<= n$
  (like 10007) makes factorials vanish: use `Lucas`. Mod 2,
  $binom(n, k)$ is odd iff `(k & ~n) == 0`.
Source: Standard; Wilson's theorem seeds the inverse table of
  `Lucas`.
Status: stress-tested against Pascal's triangle and math.comb
"""


class Binom:
  def __init__(self, n):
    fact = [1] * (n + 1)
    for i in range(1, n + 1):
      fact[i] = fact[i - 1] * i % MOD
    ifact = [1] * (n + 1)
    ifact[n] = modinv(fact[n])
    for i in range(n, 0, -1):
      ifact[i - 1] = ifact[i] * i % MOD
    self.fact, self.ifact = fact, ifact

  def C(self, n, k):
    if k < 0 or k > n:
      return 0
    fi = self.ifact
    return self.fact[n] * fi[k] % MOD * fi[n - k] % MOD

  def perm(self, n, k):
    if k < 0 or k > n:
      return 0
    return self.fact[n] * self.ifact[n - k] % MOD

  def multinom(self, ks):  # sum(ks)! / prod(k!)
    r, fi = self.fact[sum(ks)], self.ifact
    for k in ks:
      r = r * fi[k] % MOD
    return r

  def catalan(self, k):
    return (self.C(2 * k, k) - self.C(2 * k, k + 1)) % MOD

  def inv(self, i):
    return self.ifact[i] * self.fact[i - 1] % MOD


class Lucas:
  def __init__(self, p):
    fact = [1] * p
    for i in range(1, p):
      fact[i] = fact[i - 1] * i % p
    ifact = [p - 1] * p  # Wilson: (p - 1)! = -1 mod p
    for i in range(p - 1, 0, -1):
      ifact[i - 1] = ifact[i] * i % p
    self.p, self.fact, self.ifact = p, fact, ifact

  def C(self, n, k):
    if k < 0 or k > n:
      return 0
    p, f, fi, r = self.p, self.fact, self.ifact, 1
    while k:
      a, b = n % p, k % p
      if b > a:
        return 0
      r = r * f[a] % p * fi[b] % p * fi[a - b] % p
      n //= p
      k //= p
    return r
