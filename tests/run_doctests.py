"""Run the `>>>` lines of a snippet's Example field as doctests.

Usage: PYTHONPATH=tests pypy3 tests/run_doctests.py ds/segtree
"""

import doctest
import re
import sys

from testlib import CONTENT, load


def example_text(source):
  header = source.split('"""', 2)[1]
  match = re.search(r"^Example:(.*?)(?=^[A-Z][a-z]+:|\Z)", header, re.M | re.S)
  return match.group(1) if match else ""


def main():
  item_id = sys.argv[1]
  source = (CONTENT / (item_id + ".py")).read_text()
  globs = dict(vars(load(item_id)))
  parser = doctest.DocTestParser()
  test = parser.get_doctest(example_text(source), globs, item_id, None, 0)
  runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
  runner.run(test)
  failed, attempted = runner.summarize(verbose=False)
  if attempted == 0:
    sys.exit(f"{item_id}: Example has no doctest lines")
  print(f"OK doctest {item_id} ({attempted} examples)")
  sys.exit(1 if failed else 0)


main()
