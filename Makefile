# snowi: build the packets, run the stress tests, lint the content.
#   make pdf                 all four PDFs into build/pdf/
#   make test [ITEM=ds/segtree] [ONLY=cpp|py]
#   make lint [ITEM=...]     headers, 64 columns, clang-format, ruff
#   make format [ITEM=...]   apply clang-format and ruff format
#   make check               lint + test + pdf

PYTHON ?= python3
export PATH := $(HOME)/.local/bin:$(PATH)

TEST_ARGS := $(ITEM) $(if $(ONLY),--lang $(ONLY))

.PHONY: pdf pdf-separate manifest test lint format check clean

pdf:
	$(PYTHON) tools/build.py

pdf-separate:
	$(PYTHON) tools/build.py --no-bundle

manifest:
	$(PYTHON) tools/build.py --manifest-only

test:
	$(PYTHON) tools/test.py $(TEST_ARGS)

lint:
	$(PYTHON) tools/lint.py $(ITEM)

format:
	$(PYTHON) tools/lint.py --fix $(ITEM)

check: lint test pdf

clean:
	rm -rf build
