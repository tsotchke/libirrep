# Contributing to libirrep

Thanks for your interest.

## Ground rules

- Pure C11. No C++, no Python, no external dependencies beyond libc / libm and
 optional CPU SIMD intrinsics.
- Every new public function must ship with at least one positive and one
 negative test, in the same commit.
- Never land code with failing tests.
- Run `make asan ubsan` before submitting; sanitizer-clean is required.
- Formatting: `clang-format -i` with the repository's `.clang-format`. CI runs
 `make lint`.
- Keep public headers self-contained. `make test` verifies this.

## ABI discipline

- `libirrep` commits to a stable ABI within a major version.
- Changing exported symbol signatures, struct layouts, or enum values requires
 a MAJOR version bump.
- `scripts/check_abi.sh` must pass before release.

## Adding a module

Prefer editing an existing module. If a new module is truly warranted:
1. Add the header `include/irrep/<name>.h` with all intended public API.
2. Add the source `src/<name>.c` with stubs.
3. Add `tests/test_<name>.c` with at least one passing test.
4. Add the header to `include/irrep/irrep.h` and the Makefile.
5. Update `docs/API.md` and `docs/REFERENCES.md` if the module introduces
 non-trivial new formulas.
6. Open a PR with a clear description and the test output attached.

## Commit messages

Follow the style of existing commits. Imperative mood, ≤72-char subject,
body explaining *why* when it is not obvious.
