#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Produce release/<version>/ artifacts for the current host triple.
# Cross-compilation is handled by CI fan-out (see .github/workflows/release.yml);
# this script builds only for the runner it executes on.
set -euo pipefail

VERSION="${1:-$(cat VERSION)}"
UNAME_S="$(uname -s | tr '[:upper:]' '[:lower:]')"
UNAME_M="$(uname -m)"
case "${UNAME_S}" in
    darwin)            OS="macos"   ;;
    linux)             OS="linux"   ;;
    mingw*|msys*|cygwin*) OS="windows" ;;
    *)                 OS="${UNAME_S}" ;;
esac

# Normalise CPU arch names so triples read consistently across runners:
#   macOS: arm64, x86_64  (uname -m reports these directly)
#   Linux: x86_64, aarch64 — rename aarch64 → arm64 for uniformity
#   Windows: x86_64 (under MSYS2/MINGW64), aarch64 / arm64 on win-arm runners
case "${UNAME_M}" in
    aarch64) ARCH="arm64"  ;;
    amd64)   ARCH="x86_64" ;;
    *)       ARCH="${UNAME_M}" ;;
esac
TRIPLE="${OS}-${ARCH}"

ROOT="release/${VERSION}"
LIB_OUT="${ROOT}/lib/${TRIPLE}"
INC_OUT="${ROOT}/include/irrep"
mkdir -p "${LIB_OUT}" "${INC_OUT}" "${ROOT}/pkgconfig"

echo "# building libirrep ${VERSION} for ${TRIPLE}"

# Two-pass build: pass 1 generates an ABI hash against the placeholder
# binary; pass 2 bakes that hash in via -DIRREP_BAKED_ABI_HASH so runtime
# irrep_abi_hash() returns the same value `pkg-config --variable=abi_hash`
# reports.
make distclean
make lib
PASS1_HASH="$(bash scripts/generate_abi_hash.sh build/lib)"
make distclean
make lib CFLAGS_OPT="-O2 -fno-math-errno -DIRREP_BAKED_ABI_HASH=\"\\\"${PASS1_HASH}\\\"\""

cp build/lib/*.a     "${LIB_OUT}/" 2>/dev/null || true
cp build/lib/*.so*   "${LIB_OUT}/" 2>/dev/null || true
cp build/lib/*.dylib "${LIB_OUT}/" 2>/dev/null || true
cp build/lib/*.dll   "${LIB_OUT}/" 2>/dev/null || true
cp build/lib/*.dll.a "${LIB_OUT}/" 2>/dev/null || true
cp include/irrep/*.h "${INC_OUT}/"

echo "${VERSION}" > "${ROOT}/VERSION"
echo "${PASS1_HASH}" > "${ROOT}/ABI_HASH"
ABI_HASH_VALUE="${PASS1_HASH}"

# pkg-config — substitute prefix, version, ABI hash, and the per-target
# Libs line. Windows (PE/COFF) has no rpath concept and MinGW's ld warns
# on `-Wl,-rpath`, so we drop the directive there; ELF and Mach-O targets
# keep it so consumers don't need `-Wl,-rpath,${libdir}` themselves.
case "${OS}" in
    windows) PKGCONFIG_LIBS='-L${libdir} -llibirrep -lm' ;;
    *)       PKGCONFIG_LIBS='-L${libdir} -Wl,-rpath,${libdir} -llibirrep -lm' ;;
esac

sed -e "s|@PREFIX@|/usr/local|g" \
    -e "s|@VERSION@|${VERSION}|g" \
    -e "s|@ABI_HASH@|${ABI_HASH_VALUE}|g" \
    -e "s|@PKGCONFIG_LIBS@|${PKGCONFIG_LIBS}|g" \
    scripts/libirrep.pc.in > "${ROOT}/pkgconfig/libirrep.pc"

# SPDX 2.3 software bill of materials — enumerates every source file shipped
# in the release with its SHA-256 and declared license.
bash scripts/make_sbom.sh "${VERSION}" > "${ROOT}/SBOM.spdx.json"

# Checksums — last, so the SBOM and pkg-config stub are both hashed.
(cd "${ROOT}" && find . -type f \! -name CHECKSUMS -exec shasum -a 256 {} + | sort > CHECKSUMS)

# Tarball.
TAR="liblibirrep-${VERSION}-${TRIPLE}.tar.gz"
(cd release && tar -czf "${TAR}" "${VERSION}")
echo "# produced release/${TAR}"
