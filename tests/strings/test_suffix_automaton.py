from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/suffix-automaton", "strings/suffix-array")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def occurrences(s, t):
  return sum(s.startswith(t, i) for i in range(len(s)))


def check_brute(s, sam):
  """States are exactly the end-position classes; counts match."""
  n = len(s)
  endpos = {}
  for l in range(n):
    for r in range(l + 1, n + 1):
      endpos.setdefault(s[l:r], []).append(r)
  check_eq(sam.distinct(), len(endpos))
  state_of, class_of = {}, {}
  for t, pos in endpos.items():
    v = sam.find(t)
    check(v > 0)
    check(sam.len[sam.link[v]] < len(t) <= sam.len[v])
    check_eq(state_of.setdefault(tuple(pos), v), v)
    check_eq(class_of.setdefault(v, pos), pos)
  check_eq(len(sam.len), len(state_of) + 1)
  check_eq(sam.find(""), 0)
  sam.count()
  for t, pos in endpos.items():
    check_eq(sam.cnt[sam.find(t)], len(pos))
  for _ in range(10):
    t = rnd_str(rnd(1, 4), "abc")
    v = sam.find(t)
    check_eq(v >= 0, t in endpos)
    check_eq(sam.cnt[v] if v >= 0 else 0, occurrences(s, t))


def longest_common(s, t):
  return max(
    (
      j - i
      for i in range(len(t))
      for j in range(i + 1, len(t) + 1)
      if t[i:j] in s
    ),
    default=0,
  )


for _ in range(iters(1000)):
  alpha = rng.choice(["a", "ab", "abc", "az"])
  s = rnd_str(rnd(0, 14), alpha)
  sam = lib.SuffixAutomaton(s)
  check_brute(s, sam)
  on, subs, d = lib.SuffixAutomaton(), set(), 0
  for r in range(len(s)):
    for _ in range(3):  # finds between extends
      t = rnd_str(rnd(1, 3), "abcd")
      check_eq(on.find(t) >= 0, t in subs)
    on.extend(ord(s[r]) - 97)
    d += on.len[on.last] - on.len[on.link[on.last]]  # online count
    subs.update(s[l : r + 1] for l in range(r + 1))
    check_eq(on.distinct(), len(subs))
    check_eq(d, len(subs))
  check_brute(s, on)
  t = rnd_str(rnd(0, 14), alpha)  # Example: longest common substring
  v = l = best = 0
  for ch in t:
    c = ord(ch) - 97
    while v and not sam.nxt[26 * v + c]:
      v = sam.link[v]
      l = sam.len[v]
    if sam.nxt[26 * v + c]:
      v = sam.nxt[26 * v + c]
      l += 1
    best = max(best, l)
  check_eq(best, longest_common(s, t))

n = 200_000
with Timer("suffix automaton of a^2e5"):
  sam = lib.SuffixAutomaton("a" * n)
  sam.count()
check_eq(sam.distinct(), n)
check_eq(sam.cnt[sam.find("a" * 7)], n - 6)
for alpha in ("ab", "abcdefghijklmnopqrstuvwxyz"):
  s = rnd_str(n, alpha)
  with Timer(f"suffix automaton n = 2e5, alphabet {len(alpha)}"):
    sam, d = lib.SuffixAutomaton(), 0
    ln, link = sam.len, sam.link
    for ch in s:
      sam.extend(ord(ch) - 97)
      d += ln[sam.last] - ln[link[sam.last]]
    sam.count()
  check(len(sam.len) <= 2 * n - 1)
  check_eq(sam.distinct(), n * (n + 1) // 2 - sum(lib.SuffixArray(s).lcp))
  check_eq(d, sam.distinct())
  for it in range(30):
    m = rnd(1, 12)
    at = rnd(0, n - m)
    t = s[at : at + m] if it % 2 else rnd_str(m, "abc")
    v = sam.find(t)
    check_eq(sam.cnt[v] if v >= 0 else 0, occurrences(s, t))
print("OK strings/suffix-automaton")
