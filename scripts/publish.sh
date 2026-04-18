#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Stage release tarballs + checksums + signatures for publishing.
# Never pushes automatically — user invokes manually after review.
set -euo pipefail

VERSION="${1:-$(cat VERSION)}"
ROOT="release/${VERSION}"
if [ ! -d "${ROOT}" ]; then
    echo "error: ${ROOT} not found — run scripts/build_release.sh first" >&2
    exit 1
fi

STAGE="${ROOT}/publish"
mkdir -p "${STAGE}"
cp release/liblibirrep-"${VERSION}"-*.tar.gz "${STAGE}/" 2>/dev/null || true
(cd "${STAGE}" && shasum -a 256 *.tar.gz > CHECKSUMS)

if command -v gpg >/dev/null 2>&1 && [ -n "${LIBIRREP_GPG_KEY:-}" ]; then
    for f in "${STAGE}"/*.tar.gz; do
        gpg --detach-sign --armor -u "${LIBIRREP_GPG_KEY}" "${f}"
    done
fi

echo "# staged to ${STAGE} — review, then upload manually"
