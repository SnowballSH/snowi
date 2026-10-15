import subprocess
import sys

from testlib import CONTENT, check_eq

result = subprocess.run(
  [sys.executable, str(CONTENT / "contest" / "template.py")],
  input=b"3\n1 2 3\n",
  capture_output=True,
  check=True,
)
check_eq(result.stdout, b"6\n")
print("OK contest/template")
