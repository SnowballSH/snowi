#import "/typst/notes.typ": *

=== Compile and run
#cpp[#cmd("sh", "# debug build: sanitizers and a checked STL; o for timing
d() { {cxx} -std={std} -g -Wall -Wextra -Wshadow \\
  -D_GLIBCXX_DEBUG -fsanitize=address,undefined \\
  -fno-sanitize-recover $1.cpp -o $1; }
o() { {cxx} -std={std} -O2 $1.cpp -o $1; }
d {source} && ./{source} < in   # o {source} && time ./{source} < big > /dev/null")
Deep recursion needs a bigger local stack: `ulimit -s unlimited`.]
#py[#cmd("sh", "{pypy} {source}.py < in
time {pypy} {source}.py < big > /dev/null")
A crash in deep recursion can be a silent segmentation fault
(exit code 139) rather than a `RecursionError`.]

=== Stress testing
A generator takes a seed, so a failing case can be reproduced: see
Random test generators. Run until the first difference or crash,
then print the input:
#cpp[#cmd("sh", "for i in $(seq 1000); do
  ./{gen} $i > in && ./{brute} < in > want && ./{source} < in > got &&
    diff -b want got > /dev/null ||
    { echo seed $i; cat in; break; }
done")]
#py[#cmd("sh", "for i in $(seq 1000); do
  {pypy} {gen}.py $i > in && {pypy} {brute}.py < in > want &&
    {pypy} {source}.py < in > got && diff -b want got > /dev/null ||
    { echo seed $i; cat in; break; }
done")]
Interactive problems: connect a local judge with a named pipe,
#cpp[#cmd-inline("mkfifo io; ./{judge} < io | ./{source} > io")]#py[#cmd-inline("mkfifo io; ./{judge} < io | {pypy} {source}.py > io")]
(`misc/interactive`).

=== Typing hash
Each item shows the first 6 hex digits of its typing hash, which
ignores whitespace and comments#py[ but keeps the indentation
structure]: retyped code gets the same hash. Check a file with
#cpp[#cmd-inline("sh hash.sh < {source}.cpp")]#py[#cmd-inline("{python} hash.py < {source}.py")], or in vim
select the lines and type
#cpp[#cmd-inline(":w !sh hash.sh")]#py[#cmd-inline(":w !{python} hash.py")].
#cpp[#cmd("sh", "# hash.sh
{cpp} -dD -P -fpreprocessed | tr -d '[:space:]' \\
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
