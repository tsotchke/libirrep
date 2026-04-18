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
CFLAGS_OPT     = -O2 -fno-math-errno
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
else
  OS_TRIPLE    = $(UNAME_S)
  SHLIB_EXT    = dll
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
        docs lint install release release-artifacts clean distclean \
        dirs print-config

all: lib test examples

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

$(STATIC_LIB): $(LIB_OBJS) | $(LIB_DIR)
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(SHARED_LIB): $(LIB_OBJS) | $(LIB_DIR)
	$(CC) $(SHLIB_FLAGS) -o $@ $^ $(LDFLAGS)

$(SHARED_LIB_LINK): $(SHARED_LIB)
	@cd $(LIB_DIR) && ln -sf $(notdir $(SHARED_LIB)) $(SHLIB_LINK)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
$(BIN_DIR)/test_%: tests/test_%.c $(STATIC_LIB) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) -Itests $< $(STATIC_LIB) $(LDFLAGS) -o $@

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
# Docs / lint
# ---------------------------------------------------------------------------
docs:
	@command -v doxygen >/dev/null 2>&1 || { \
	    echo "doxygen not installed; skipping"; exit 0; }
	@mkdir -p $(BUILD_DIR)/docs && cd $(BUILD_DIR)/docs && doxygen $(CURDIR)/Doxyfile

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
