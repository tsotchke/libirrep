# SPDX-License-Identifier: MIT
# Makefile for libirrep: static + shared library, versioned install, sanitizer
# builds, benchmarks, examples, and release artifacts.
#
# Targets: all, lib, lib-static, lib-shared, test, bench, examples,
#          asan, ubsan, docs, lint, install, release, clean, distclean.

# ---------------------------------------------------------------------------
# Toolchain
# ---------------------------------------------------------------------------
CC       ?= cc
AR       ?= ar
RANLIB   ?= ranlib
UNAME_S  := $(shell uname -s)
UNAME_M  := $(shell uname -m)

VERSION  := $(shell cat VERSION)
SO_MAJOR := 1

# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
CFLAGS_COMMON  = -std=c11 -Iinclude -fvisibility=hidden
CFLAGS_WARN    = -Wall -Wextra -Wpedantic -Wno-unused-parameter
# -ffp-contract=on: force consistent a*b+c -> fma contraction across
# compilers so the NEON / AVX2 kernels (which use explicit fma
# intrinsics) stay bit-exact to the scalar reference under both clang
# (contract-on by default) and gcc (contract-off by default at -O2).
CFLAGS_OPT     = -O2 -fno-math-errno -ffp-contract=on
CFLAGS         = $(CFLAGS_COMMON) $(CFLAGS_WARN) $(CFLAGS_OPT)
LDFLAGS        = -lm

# Position-independent code for the shared library.
PIC            = -fPIC

# Arch-specific CFLAGS (runtime dispatch happens in src/simd_runtime.c; the
# macro USE_NEON_IF_AVAILABLE mirrors the spin_based_neural_network pattern).
ifeq ($(UNAME_M),arm64)
  ARCH_CFLAGS  = -DUSE_NEON_IF_AVAILABLE
  ARCH_TRIPLE  = arm64
else ifeq ($(UNAME_M),aarch64)
  ARCH_CFLAGS  = -DUSE_NEON_IF_AVAILABLE
  ARCH_TRIPLE  = aarch64
else ifeq ($(UNAME_M),x86_64)
  ARCH_CFLAGS  = -DUSE_SSE42_IF_AVAILABLE -DUSE_AVX2_IF_AVAILABLE
  ARCH_TRIPLE  = x86_64
else
  ARCH_CFLAGS  =
  ARCH_TRIPLE  = $(UNAME_M)
endif

ifeq ($(UNAME_S),Darwin)
  OS_TRIPLE    = macos
  SHLIB_EXT    = dylib
  SHLIB_SONAME = liblibirrep.$(SO_MAJOR).$(SHLIB_EXT)
  SHLIB_LINK   = liblibirrep.$(SHLIB_EXT)
  SHLIB_FLAGS  = -dynamiclib -install_name @rpath/$(SHLIB_SONAME) \
                 -compatibility_version $(SO_MAJOR) \
                 -current_version $(VERSION)
else ifeq ($(UNAME_S),Linux)
  OS_TRIPLE    = linux
  SHLIB_EXT    = so
  SHLIB_SONAME = liblibirrep.$(SHLIB_EXT).$(SO_MAJOR)
  SHLIB_LINK   = liblibirrep.$(SHLIB_EXT)
  SHLIB_FLAGS  = -shared -Wl,-soname,$(SHLIB_SONAME)
else ifneq (,$(findstring MINGW,$(UNAME_S)))
  OS_TRIPLE    = windows
  SHLIB_EXT    = dll
  SHLIB_SONAME = liblibirrep.$(SHLIB_EXT)
  SHLIB_LINK   = liblibirrep.$(SHLIB_EXT)
  # MinGW-w64 export control: DLL path emits import library via
  # --out-implib so consumers can link; dllexport attributes on public
  # symbols come from IRREP_API when IRREP_BUILDING_DLL is defined.
  SHLIB_FLAGS  = -shared -Wl,--out-implib,$(LIB_DIR)/liblibirrep.dll.a
  PIC          =
  CFLAGS      += -DIRREP_BUILDING_DLL
else ifneq (,$(findstring CYGWIN,$(UNAME_S)))
  OS_TRIPLE    = windows
  SHLIB_EXT    = dll
  SHLIB_SONAME = liblibirrep.$(SHLIB_EXT)
  SHLIB_LINK   = liblibirrep.$(SHLIB_EXT)
  SHLIB_FLAGS  = -shared -Wl,--out-implib,$(LIB_DIR)/liblibirrep.dll.a
  PIC          =
  CFLAGS      += -DIRREP_BUILDING_DLL
else
  # Unknown OS — fall back to POSIX-ish shared lib.
  OS_TRIPLE    = $(UNAME_S)
  SHLIB_EXT    = so
  SHLIB_SONAME = liblibirrep.$(SHLIB_EXT)
  SHLIB_LINK   = liblibirrep.$(SHLIB_EXT)
  SHLIB_FLAGS  = -shared
endif

TRIPLE = $(OS_TRIPLE)-$(ARCH_TRIPLE)

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
BUILD_DIR  = build
OBJ_DIR    = $(BUILD_DIR)/obj
LIB_DIR    = $(BUILD_DIR)/lib
BIN_DIR    = $(BUILD_DIR)/bin

$(OBJ_DIR) $(LIB_DIR) $(BIN_DIR):
	@mkdir -p $@

# ---------------------------------------------------------------------------
# Source file discovery
# ---------------------------------------------------------------------------
LIB_SRCS := $(wildcard src/*.c)
LIB_OBJS := $(patsubst src/%.c,$(OBJ_DIR)/src/%.o,$(LIB_SRCS))

# SIMD source files compile with architecture-specific flags (runtime dispatch
# decides whether to call them). For M1 these files do not yet exist; guard
# with wildcard so the rule is a no-op until M10/M11 add them.
NEON_SRCS := $(wildcard src/*_neon.c)
AVX2_SRCS := $(wildcard src/*_avx2.c)

TEST_SRCS := $(wildcard tests/test_*.c)
TEST_BINS := $(patsubst tests/%.c,$(BIN_DIR)/%,$(TEST_SRCS))

BENCH_SRCS := $(wildcard benchmarks/bench_*.c)
BENCH_BINS := $(patsubst benchmarks/%.c,$(BIN_DIR)/%,$(BENCH_SRCS))

EXAMPLE_SRCS := $(wildcard examples/*.c)
EXAMPLE_BINS := $(patsubst examples/%.c,$(BIN_DIR)/%,$(EXAMPLE_SRCS))

STATIC_LIB := $(LIB_DIR)/liblibirrep.a
SHARED_LIB := $(LIB_DIR)/$(SHLIB_SONAME)
SHARED_LIB_LINK := $(LIB_DIR)/$(SHLIB_LINK)

# ---------------------------------------------------------------------------
# Default / phony targets
# ---------------------------------------------------------------------------
.PHONY: all lib lib-static lib-shared test bench examples asan ubsan \
        coverage fuzz fuzz-driver fuzz-run docs lint check-headers check-abi \
        install release release-artifacts clean distclean dirs print-config

all: lib test examples check-headers check-abi

print-config:
	@echo "CC         = $(CC)"
	@echo "UNAME_S    = $(UNAME_S)"
	@echo "UNAME_M    = $(UNAME_M)"
	@echo "TRIPLE     = $(TRIPLE)"
	@echo "VERSION    = $(VERSION)"
	@echo "CFLAGS     = $(CFLAGS) $(ARCH_CFLAGS)"
	@echo "SHLIB      = $(SHARED_LIB)"
	@echo "STATIC_LIB = $(STATIC_LIB)"

# ---------------------------------------------------------------------------
# Library build rules
# ---------------------------------------------------------------------------
lib: lib-static lib-shared

lib-static: $(STATIC_LIB)
lib-shared: $(SHARED_LIB_LINK)

$(OBJ_DIR)/src/%.o: src/%.c | $(OBJ_DIR)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) $(PIC) -c $< -o $@

# Arch-specific SIMD TUs need their own per-file flags so the intrinsics
# compile on any host, even when the default build arch lacks them.
$(OBJ_DIR)/src/%_avx2.o: src/%_avx2.c | $(OBJ_DIR)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) $(PIC) \
	    $(if $(filter x86_64,$(UNAME_M)),-mavx2 -mfma) \
	    -c $< -o $@

$(STATIC_LIB): $(LIB_OBJS) | $(LIB_DIR)
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(SHARED_LIB): $(LIB_OBJS) | $(LIB_DIR)
	$(CC) $(SHLIB_FLAGS) -o $@ $^ $(LDFLAGS)

# On MinGW the SONAME and the bare link name are the same file
# (liblibirrep.dll), so there's nothing to symlink — the shared lib itself
# satisfies the target. On ELF / Mach-O the link name (liblibirrep.so /
# liblibirrep.dylib) is a versioned-symlink companion to the SONAME.
ifeq ($(SHLIB_SONAME),$(SHLIB_LINK))
$(SHARED_LIB_LINK): $(SHARED_LIB)
	@:
else
$(SHARED_LIB_LINK): $(SHARED_LIB)
	@cd $(LIB_DIR) && ln -sf $(notdir $(SHARED_LIB)) $(SHLIB_LINK)
endif

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
$(BIN_DIR)/test_%: tests/test_%.c $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) -Itests $< $(STATIC_LIB) $(LDFLAGS) -o $@

# test_downstream_compat links the lattice loader too; explicit rule
# overrides the generic %.c pattern and takes priority in make resolution.
$(BIN_DIR)/test_downstream_compat: tests/test_downstream_compat.c \
                                    tests/test_downstream_compat/lattice_loader.c \
                                    tests/test_downstream_compat/lattice_loader.h \
                                    $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) -Itests -Itests/test_downstream_compat \
	    tests/test_downstream_compat.c \
	    tests/test_downstream_compat/lattice_loader.c \
	    $(STATIC_LIB) $(LDFLAGS) -o $@

test: $(TEST_BINS)
	@echo "# Running $(words $(TEST_BINS)) test suites"
	@pass=0; fail=0; \
	for t in $(TEST_BINS); do \
	    if $$t; then pass=$$((pass+1)); \
	    else fail=$$((fail+1)); echo "# FAIL: $$t"; fi; \
	done; \
	echo "# $$pass passed, $$fail failed, $(words $(TEST_BINS)) total"; \
	test $$fail -eq 0

# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------
$(BIN_DIR)/bench_%: benchmarks/bench_%.c $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) -Ibenchmarks $< $(STATIC_LIB) $(LDFLAGS) -o $@

# Benchmarks that pull in downstream-compat helpers (lattice loader, golden
# fixtures). Explicit rule so the extra TU participates in the link.
$(BIN_DIR)/bench_downstream_shapes: benchmarks/bench_downstream_shapes.c \
                                     tests/test_downstream_compat/lattice_loader.c \
                                     tests/test_downstream_compat/lattice_loader.h \
                                     $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) -Ibenchmarks -Itests/test_downstream_compat \
	    benchmarks/bench_downstream_shapes.c \
	    tests/test_downstream_compat/lattice_loader.c \
	    $(STATIC_LIB) $(LDFLAGS) -o $@

bench: $(BENCH_BINS)
	@bash benchmarks/run_benchmarks.sh

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------
$(BIN_DIR)/%: examples/%.c $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) $< $(STATIC_LIB) $(LDFLAGS) -o $@

examples: $(EXAMPLE_BINS)

# ---------------------------------------------------------------------------
# Sanitizer rebuilds — wipes build dir then re-runs with -fsanitize=...
# ---------------------------------------------------------------------------
asan:
	$(MAKE) distclean
	$(MAKE) CFLAGS_OPT="-O1 -fno-math-errno -fsanitize=address -fno-omit-frame-pointer -g" \
	        LDFLAGS="-lm -fsanitize=address" test

ubsan:
	$(MAKE) distclean
	$(MAKE) CFLAGS_OPT="-O1 -fno-math-errno -fsanitize=undefined -fno-omit-frame-pointer -g" \
	        LDFLAGS="-lm -fsanitize=undefined" test

# ---------------------------------------------------------------------------
# Coverage. Rebuilds with LLVM source-based coverage instrumentation, runs
# every test binary with a unique LLVM_PROFILE_FILE, merges the .profraw
# stream into a .profdata blob, and prints a per-file report plus a
# summary for the src/ tree. Requires `llvm-profdata` and `llvm-cov` in
# PATH — both ship with Apple Xcode (`xcrun --find ...`) and Homebrew
# LLVM.
# ---------------------------------------------------------------------------
COV_DIR := $(BUILD_DIR)/coverage
LLVM_PROFDATA ?= $(shell xcrun --find llvm-profdata 2>/dev/null || command -v llvm-profdata)
LLVM_COV      ?= $(shell xcrun --find llvm-cov      2>/dev/null || command -v llvm-cov)

coverage:
	@command -v "$(LLVM_PROFDATA)" >/dev/null 2>&1 || { \
	    echo "llvm-profdata not found; install LLVM or Xcode"; exit 1; }
	$(MAKE) distclean
	$(MAKE) CFLAGS_OPT="-O1 -fno-math-errno -fprofile-instr-generate -fcoverage-mapping -g" \
	        LDFLAGS="-lm -fprofile-instr-generate -fcoverage-mapping" lib test bench
	@mkdir -p $(COV_DIR)/raw
	@rm -f $(COV_DIR)/raw/*.profraw
	@# Run each test binary with its own profraw output
	@for t in $(TEST_BINS) $(BENCH_BINS); do \
	    name=$$(basename $$t); \
	    LLVM_PROFILE_FILE="$(COV_DIR)/raw/$$name.profraw" $$t > /dev/null 2>&1 || :; \
	done
	@"$(LLVM_PROFDATA)" merge -sparse $(COV_DIR)/raw/*.profraw \
	    -o $(COV_DIR)/merged.profdata
	@"$(LLVM_COV)" report \
	    -instr-profile=$(COV_DIR)/merged.profdata \
	    $(firstword $(TEST_BINS)) $(addprefix -object=,$(wordlist 2,$(words $(TEST_BINS)),$(TEST_BINS))) \
	    src/ > $(COV_DIR)/report.txt 2>&1 || :
	@awk '/^TOTAL/ { \
	    printf "# coverage: %s line, %s function, %s region, %s branch\n", \
	        $$10, $$7, $$4, $$13 \
	  }' $(COV_DIR)/report.txt | tee $(COV_DIR)/summary.txt
	@echo "# per-file report at $(COV_DIR)/report.txt"

# ---------------------------------------------------------------------------
# Header self-containedness: each public header must compile standalone
# under strict warnings, with no transitive include dependencies.
# ---------------------------------------------------------------------------
HEADER_SRCS   := $(wildcard include/irrep/*.h)
HEADER_CHECKS := $(patsubst include/irrep/%.h,$(OBJ_DIR)/hdr/%.ok,$(HEADER_SRCS))

$(OBJ_DIR)/hdr/%.ok: include/irrep/%.h | $(OBJ_DIR)
	@mkdir -p $(dir $@)
	@printf '/* auto-generated: header self-containedness check */\n#include <irrep/%s.h>\ntypedef int irrep_hdr_ok_%s_t;\n' '$*' '$*' > $(dir $@)/$*.c
	@$(CC) -std=c11 -Wall -Wextra -Wpedantic -Werror -Iinclude \
	    -fsyntax-only $(dir $@)/$*.c
	@touch $@

check-headers: $(HEADER_CHECKS)
	@echo "# $(words $(HEADER_CHECKS)) public headers compile self-contained"

# ---------------------------------------------------------------------------
# Fuzz targets.
#
# `make fuzz`             — builds libFuzzer-instrumented binaries (needs a
#                           clang that ships libclang_rt.fuzzer_osx.a; run
#                           e.g. ./build/bin/fuzz_cg -max_total_time=60).
# `make fuzz-driver`      — builds libFuzzer-less binaries that link the
#                           same harnesses against tests/fuzz/driver.c, a
#                           deterministic PCG-style random-byte driver.
#                           Portable across Apple clang / Homebrew clang /
#                           GCC; the preferred path on systems where
#                           libFuzzer isn't available.
# `make fuzz-run`         — builds the driver binaries and runs each for
#                           1_000_000 iterations with seed 42. Writes a
#                           report to build/fuzz-report.txt.
# ---------------------------------------------------------------------------
FUZZ_SRCS := $(wildcard tests/fuzz/fuzz_*.c)
FUZZ_BINS := $(patsubst tests/fuzz/%.c,$(BIN_DIR)/%,$(FUZZ_SRCS))
FUZZ_DRIVER_BINS := $(patsubst tests/fuzz/%.c,$(BIN_DIR)/%_driver,$(FUZZ_SRCS))
FUZZ_FLAGS = -fsanitize=fuzzer,address,undefined -O1 -g
FUZZ_DRIVER_FLAGS = -fsanitize=address,undefined -O1 -g

$(BIN_DIR)/fuzz_%: tests/fuzz/fuzz_%.c $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS_WARN) $(ARCH_CFLAGS) $(FUZZ_FLAGS) \
	    $< $(STATIC_LIB) -lm -o $@

$(BIN_DIR)/fuzz_%_driver: tests/fuzz/fuzz_%.c tests/fuzz/driver.c $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS_WARN) $(ARCH_CFLAGS) $(FUZZ_DRIVER_FLAGS) \
	    $< tests/fuzz/driver.c $(STATIC_LIB) -lm -o $@

fuzz: $(FUZZ_BINS)
fuzz-driver: $(FUZZ_DRIVER_BINS)

# ---------------------------------------------------------------------------
# ABI regression gate. The first `check-abi` run after a clone writes
# release/BASELINE_ABI_HASH; subsequent runs compare the current hash to it
# and fail the build on drift. Bumping MAJOR in VERSION and refreshing the
# baseline is the intentional-break workflow.
# ---------------------------------------------------------------------------
check-abi: lib
	@bash scripts/check_abi.sh

fuzz-run: fuzz-driver
	@mkdir -p $(BUILD_DIR)
	@: > $(BUILD_DIR)/fuzz-report.txt
	@for b in $(FUZZ_DRIVER_BINS); do \
	    echo "## $$b" >> $(BUILD_DIR)/fuzz-report.txt; \
	    $$b 42 1000000 256 >> $(BUILD_DIR)/fuzz-report.txt 2>&1 || \
	        { echo "FAILED: $$b"; exit 1; }; \
	done
	@echo "# fuzz run ok; report in $(BUILD_DIR)/fuzz-report.txt"

# ---------------------------------------------------------------------------
# Docs / lint
# ---------------------------------------------------------------------------
docs:
	@command -v doxygen >/dev/null 2>&1 || { \
	    echo "doxygen not installed; skipping"; exit 0; }
	@mkdir -p $(BUILD_DIR)/docs && doxygen Doxyfile

lint:
	@command -v clang-format >/dev/null 2>&1 || { \
	    echo "clang-format not installed; skipping"; exit 0; }
	@clang-format --dry-run --Werror $(LIB_SRCS) $(TEST_SRCS) \
	    $(wildcard include/irrep/*.h)

# ---------------------------------------------------------------------------
# Install / release
# ---------------------------------------------------------------------------
PREFIX ?= /usr/local
DESTDIR ?=

install: lib
	install -d $(DESTDIR)$(PREFIX)/include/irrep
	install -m 644 include/irrep/*.h $(DESTDIR)$(PREFIX)/include/irrep/
	install -d $(DESTDIR)$(PREFIX)/lib
	install -m 644 $(STATIC_LIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(SHARED_LIB) $(DESTDIR)$(PREFIX)/lib/
	cd $(DESTDIR)$(PREFIX)/lib && ln -sf $(notdir $(SHARED_LIB)) $(SHLIB_LINK)
	install -d $(DESTDIR)$(PREFIX)/lib/pkgconfig
	sed -e 's|@PREFIX@|$(PREFIX)|g' -e 's|@VERSION@|$(VERSION)|g' \
	    scripts/libirrep.pc.in > $(DESTDIR)$(PREFIX)/lib/pkgconfig/libirrep.pc

release: release-artifacts

release-artifacts:
	@bash scripts/build_release.sh $(VERSION)

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

distclean:
	@rm -rf $(BUILD_DIR)
