#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Orchestrator: bump VERSION, refresh CHANGELOG, build artifacts, tag.
# Does not push. Full implementation lands in M13.
set -euo pipefail

NEW_VERSION="${1:?usage: scripts/release.sh <new-version>}"

echo "${NEW_VERSION}" > VERSION
echo "# VERSION -> ${NEW_VERSION}"

bash scripts/build_release.sh "${NEW_VERSION}"
bash scripts/check_abi.sh

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git add VERSION CHANGELOG.md
    git commit -m "release: v${NEW_VERSION}" || true
    git tag -a "v${NEW_VERSION}" -m "libirrep v${NEW_VERSION}"
    echo "# tagged v${NEW_VERSION} (not pushed)"
fi
