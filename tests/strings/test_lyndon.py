from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/lyndon")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def is_lyndon(w):
  return w != "" and all(w[i:] > w for i in range(1, len(w)))


def check_duval(s, st):  # with uniqueness these pin the answer down
  n = len(s)
  check(st == [] if n == 0 else st[0] == 0)
  ends = st[1:] + [n]
  prev = None
  for a, b in zip(st, ends):
    check(a < b)
    w = s[a:b]
    check(len(w) > 64 or is_lyndon(w), w)
    check(prev is None or w <= prev)
    prev = w


def brute_rotation(s):
  rots = [s[i:] + s[:i] for i in range(len(s))]
  return rots.index(min(rots)) if s else 0


for _ in range(iters(3000)):
  s = rnd_str(rnd(0, 20), rng.choice(["a", "ab", "abc"]))
  n = len(s)
  if n and rnd(0, 3) == 0:  # periodic: w^m w'
    p = rnd(1, 4)
    s = (s[:p] * n)[:n]
  st = lib.duval(s)
  check_duval(s, st)
  ends = st + [n]  # Usage: starts only, append n for the last end
  check_eq("".join(s[ends[j] : ends[j + 1]] for j in range(len(st))), s)
  if n:
    check_eq(s[st[-1] :], min(s[i:] for i in range(n)))
  check_eq(len(st) == 1, is_lyndon(s))
  r = lib.min_rotation(s)
  check_eq(r, brute_rotation(s), s)
  check_eq(lib.duval(s.encode()), st)
  big = [(ord(c) - 98) * 10**30 for c in s]
  check_eq(lib.duval(big), st)
  check_eq(lib.min_rotation(big), r)
  # Tweak: largest rotation from the flipped order
  f = "".join(chr(219 - ord(c)) for c in s)
  rf = lib.min_rotation(f)
  check(n == 0 or s[rf:] + s[:rf] == max(s[i:] + s[:i] for i in range(n)))

n = 10**6
s = rnd_str(n, "ab")
with Timer("duval n = 1e6"):
  st = lib.duval(s)
check_duval(s, st)
s = list(rnd_str(n, "bc"))
pos = rnd(0, n - 1)
for i in range(20):  # a unique run of 20 a's starts the minimum
  s[(pos + i) % n] = "a"
s = "".join(s)
with Timer("min_rotation n = 1e6"):
  r = lib.min_rotation(s)
check_eq(r, pos)
per = "aab" * (n // 3)
check_eq(lib.min_rotation(per[1:] + per[0]), 2)
with Timer("duval of a^1e6, bytes"):
  check_eq(len(lib.duval(b"a" * n)), n)
print("OK strings/lyndon")
