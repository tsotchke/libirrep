#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Produce release/<version>/ artifacts for the current host triple.
# Full per-triple cross-compilation lands in M13.
set -euo pipefail

VERSION="${1:-$(cat VERSION)}"
UNAME_S="$(uname -s | tr '[:upper:]' '[:lower:]')"
UNAME_M="$(uname -m)"
case "${UNAME_S}" in
    darwin) OS="macos" ;;
    linux)  OS="linux" ;;
    *)      OS="${UNAME_S}" ;;
esac
TRIPLE="${OS}-${UNAME_M}"

ROOT="release/${VERSION}"
LIB_OUT="${ROOT}/lib/${TRIPLE}"
INC_OUT="${ROOT}/include/irrep"
mkdir -p "${LIB_OUT}" "${INC_OUT}" "${ROOT}/pkgconfig"

echo "# building libirrep ${VERSION} for ${TRIPLE}"
make lib

cp build/lib/*.a   "${LIB_OUT}/" 2>/dev/null || true
cp build/lib/*.so* "${LIB_OUT}/" 2>/dev/null || true
cp build/lib/*.dylib "${LIB_OUT}/" 2>/dev/null || true
cp include/irrep/*.h "${INC_OUT}/"

echo "${VERSION}" > "${ROOT}/VERSION"
bash scripts/generate_abi_hash.sh "${LIB_OUT}" > "${ROOT}/ABI_HASH"

# Checksums.
(cd "${ROOT}" && find . -type f \! -name CHECKSUMS -exec shasum -a 256 {} + | sort > CHECKSUMS)

# pkg-config stub.
sed -e "s|@PREFIX@|/usr/local|g" -e "s|@VERSION@|${VERSION}|g" \
    scripts/libirrep.pc.in > "${ROOT}/pkgconfig/libirrep.pc"

# Tarball.
TAR="liblibirrep-${VERSION}-${TRIPLE}.tar.gz"
(cd release && tar -czf "${TAR}" "${VERSION}")
echo "# produced release/${TAR}"
