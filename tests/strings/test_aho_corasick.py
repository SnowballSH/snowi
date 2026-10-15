from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("strings/aho-corasick")


def rnd_str(n, alpha):
  return "".join(rng.choice(alpha) for _ in range(n))


def occurrences(s, t):
  return sum(s.startswith(t, i) for i in range(len(s)))


def check_automaton(pats, ac):
  """Moves, links and words against the definitions."""
  state = {}
  for p in pats:
    for k in range(len(p) + 1):
      v = 0
      for ch in p[:k]:
        v = ac.nxt[26 * v + ord(ch) - 97]
      check_eq(state.setdefault(p[:k], v), v)
  check_eq(len(state), len(ac.link))
  check_eq(sorted(ac.q), list(range(len(ac.link))))

  def longest_suffix(x, start):
    return next(
      state[x[k:]] for k in range(start, len(x) + 1) if x[k:] in state
    )

  for x, v in state.items():
    for c in range(26):
      want = longest_suffix(x + chr(97 + c), 0)
      check_eq(ac.nxt[26 * v + c], want)
    if v:
      check_eq(ac.link[v], longest_suffix(x, 1))
    check_eq(ac.words[v], sum(x.endswith(p) for p in pats))
  check_eq(ac.node, [state[p] for p in pats])


for _ in range(iters(1000)):
  alpha = rng.choice(["a", "ab", "abc", "az"])
  pats = [rnd_str(rnd(1, 5), alpha) for _ in range(rnd(1, 8))]
  if rnd(0, 3) == 0:
    pats.append(pats[0])
  ac = lib.AhoCorasick(pats)
  check_automaton(pats, ac)
  t = rnd_str(rnd(0, 40), alpha)
  occ = ac.count(t)
  check_eq(occ, [occurrences(t, p) for p in pats])
  v = total = 0  # Usage: matches of all patterns
  for ch in t:
    v = ac.nxt[26 * v + ord(ch) - 97]
    total += ac.words[v]
  check_eq(total, sum(occ))
  m = len(ac.link)  # Example: list every match
  pid, out = [-1] * m, [-1] * m
  for i, u in enumerate(ac.node):
    pid[u] = i
  for u in ac.q:
    out[u] = u if pid[u] >= 0 else out[ac.link[u]]
  v, got = 0, set()
  for i, ch in enumerate(t):
    v = ac.nxt[26 * v + ord(ch) - 97]
    u = out[v]
    while u >= 0:
      check((u, i) not in got)
      got.add((u, i))
      u = out[ac.link[u]]
  want = {
    (ac.node[j], e)
    for j, p in enumerate(pats)
    for e in range(len(t))
    if t[: e + 1].endswith(p)
  }
  check_eq(got, want)

check_eq(lib.AhoCorasick(["a"]).count(""), [0])
check_eq(lib.AhoCorasick(["a", "a"]).count("aaa"), [3, 3])
letters = [chr(97 + c) for c in range(26)]
pangram = "thequickbrownfoxjumpsoverthelazydog"
check_eq(
  lib.AhoCorasick(letters).count(pangram), [pangram.count(c) for c in letters]
)

n = 200_000
for alpha in ("ab", letters):  # the 26-letter trie is the slow one
  pats, total_len = [], 0
  while total_len < n:
    pats.append(rnd_str(rnd(1, 20), alpha))
    total_len += len(pats[-1])
  t = rnd_str(n, alpha)
  with Timer(f"Aho-Corasick L = |t| = 2e5, alphabet {len(alpha)}"):
    occ = lib.AhoCorasick(pats).count(t)
  for i in range(100):
    check_eq(occ[i], occurrences(t, pats[i]))
with Timer("Aho-Corasick a, aa, ..., a^600 on a^2e5"):
  occ = lib.AhoCorasick(["a" * k for k in range(1, 601)]).count("a" * n)
check_eq(occ, [n - k + 1 for k in range(1, 601)])
print("OK strings/aho-corasick")
