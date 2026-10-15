from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/eertree")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def occurrences(s, t):
  return sum(s.startswith(t, i) for i in range(len(s)))


def palindromes(s):
  n = len(s)
  return [
    s[l:r]
    for l in range(n)
    for r in range(l + 1, n + 1)
    if s[l:r] == s[l:r][::-1]
  ]


def check_brute(s, et):
  """Node strings from the edges (c + v + c) against the definitions."""
  n, m = len(s), len(et.len)
  pals = palindromes(s)
  check_eq(m - 2, len(set(pals)))
  check_eq((et.len[0], et.len[1]), (-1, 0))
  name = ["", ""] + [None] * (m - 2)
  for v in range(m):
    for c in range(26):
      w = et.nxt[26 * v + c]
      if w:
        check(w > v and name[w] is None)
        ch = chr(97 + c)
        name[w] = ch if v == 0 else ch + name[v] + ch
  check_eq(set(name[2:]), set(pals))
  node = {x: v for v, x in enumerate(name) if v >= 2}
  for v in range(2, m):
    x = name[v]
    check_eq(et.len[v], len(x))
    want = next((node[x[k:]] for k in range(1, len(x)) if x[k:] in node), 1)
    check_eq(et.link[v], want)
  check_eq(len(et.suf), n)
  for i in range(n):
    best = next(
      i + 1 - l for l in range(i + 1) if s[l : i + 1] == s[l : i + 1][::-1]
    )
    check_eq(name[et.suf[i]], s[i + 1 - best : i + 1])
  et.count()
  for v in range(2, m):
    check_eq(et.cnt[v], occurrences(s, name[v]))
  best = max((et.cnt[v] * et.len[v] for v in range(2, m)), default=0)
  check_eq(best, max((occurrences(s, x) * len(x) for x in pals), default=0))
  dep = [0] * m  # Example: palindromes with repeats
  for v in range(2, m):
    dep[v] = dep[et.link[v]] + 1
  check_eq(sum(dep[v] for v in et.suf), len(pals))


for _ in range(iters(1500)):
  s = rnd_str(rnd(0, 20), rng.choice(["a", "ab", "abc", "az"]))
  check_brute(s, lib.Eertree(s))
  on = lib.Eertree()
  for i in range(len(s)):
    on.add(ord(s[i]) - 97)
    check_eq(len(on.len) - 2, len(set(palindromes(s[: i + 1]))))
  on.count()  # count() twice, after more letters
  on.add(0)
  check_brute(s + "a", on)


def expand(s):
  """All palindromes by center expansion: total and longest ends."""
  m, total = len(s), 0
  longest = [0] * m
  for c in range(2 * m - 1):
    l, r = c // 2, (c + 1) // 2
    while l >= 0 and r < m and s[l] == s[r]:
      total += 1
      longest[r] = max(longest[r], r - l + 1)
      l -= 1
      r += 1
  return total, longest


n = 200_000
with Timer("eertree of a^2e5"):
  et = lib.Eertree("a" * n)
  et.count()
check_eq(len(et.len) - 2, n)
check(all(et.cnt[v] == n - et.len[v] + 1 for v in range(2, n + 2)))
check_eq(max(et.cnt[v] * et.len[v] for v in range(2, n + 2)), 10000100000)
f, g = "a", "ab"
while len(g) < n:
  f, g = g, g + f
for name, s in [
  ("Fibonacci", g),
  ("random 2", rnd_str(n, "ab")),
  ("random 26", rnd_str(n, "abcdefghijklmnopqrstuvwxyz")),
]:
  with Timer(f"eertree n = {len(s)}, {name}"):
    et = lib.Eertree(s)
    et.count()
  total, longest = expand(s)
  check(len(et.len) - 2 <= len(s))
  check_eq(sum(et.cnt[2:]), total)
  check_eq([et.len[v] for v in et.suf], longest)
print("OK strings/eertree")
