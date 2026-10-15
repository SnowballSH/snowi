#import "/typst/notes.typ": *

=== Compile and run
#cpp[#code("sh", "# debug build: sanitizers and a checked STL; o for timing
d() { g++ -std=gnu++20 -g -Wall -Wextra -Wshadow \\
  -D_GLIBCXX_DEBUG -fsanitize=address,undefined \\
  -fno-sanitize-recover $1.cpp -o $1; }
o() { g++ -std=gnu++20 -O2 $1.cpp -o $1; }
d a && ./a < in          # o a && time ./a < big > /dev/null")
Deep recursion needs a bigger local stack: `ulimit -s unlimited`.]
#py[#code("sh", "pypy3 a.py < in
time pypy3 a.py < big > /dev/null")
A crash in deep recursion can be a silent segmentation fault
(exit code 139) rather than a `RecursionError`.]

=== Stress testing
A generator takes a seed, so a failing case can be reproduced:
#lang-code(
  cpp: "int main(int, char** argv) {
  mt19937 rng(atoi(argv[1]));
  auto rnd = [&](int lo, int hi) {
    return uniform_int_distribution(lo, hi)(rng);
  };
  int n = rnd(1, 8);
  printf(\"%d\\n\", n);
  rep(i, 0, n) printf(\"%d \", rnd(1, 10));
  puts(\"\");
}",
  py: "import random, sys
random.seed(int(sys.argv[1]))
n = random.randint(1, 8)
print(n)
print(*[random.randint(1, 10) for _ in range(n)])",
)
Run until the first difference or crash, then print the input:
#cpp[#code("sh", "for i in $(seq 1000); do
  ./gen $i > in && ./brute < in > want && ./a < in > got &&
    diff -b want got > /dev/null ||
    { echo seed $i; cat in; break; }
done")]
#py[#code("sh", "for i in $(seq 1000); do
  pypy3 gen.py $i > in && pypy3 brute.py < in > want &&
    pypy3 a.py < in > got && diff -b want got > /dev/null ||
    { echo seed $i; cat in; break; }
done")]
Interactive problems: connect a local judge with a named pipe,
#cpp[#box[`mkfifo io; ./judge < io | ./a > io`]]#py[#box[`mkfifo io; ./judge < io | pypy3 a.py > io`]]
(`misc/interactive`).

=== Typing hash
Each item shows the first 6 hex digits of its typing hash, which
ignores whitespace and comments#py[ but keeps the indentation
structure]: retyped code gets the same hash. Check a file with
#cpp[#box[`sh hash.sh < a.cpp`]]#py[#box[`python3 hash.py < a.py`]], or in vim
select the lines and type
#cpp[#box[`:w !sh hash.sh`]]#py[#box[`:w !python3 hash.py`]].
#cpp[#code("sh", "# hash.sh
cpp -dD -P -fpreprocessed | tr -d '[:space:]' \\
  | md5sum | cut -c-6")]
#py[#code("py", "# hash.py
import hashlib, sys, tokenize as tk
skip = tk.COMMENT, tk.NL, tk.ENCODING, tk.ENDMARKER
mark = {tk.INDENT: \"{\", tk.DEDENT: \"}\", tk.NEWLINE: \";\"}
toks = tk.tokenize(sys.stdin.buffer.readline)
s = \"\".join(mark.get(t.type, t.string) for t in toks
            if t.type not in skip)
print(hashlib.md5(s.encode()).hexdigest()[:6])")]

=== .vimrc
#code("", "set nu ai et sw=2 ts=2 mouse=a | syntax on")
