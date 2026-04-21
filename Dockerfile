# SPDX-License-Identifier: MIT
# Reproducible Linux build environment for libirrep.
#
# This Dockerfile mirrors the ubuntu-22.04 CI runner so anyone can
# reproduce the Linux test + coverage + fuzz runs locally without
# matching toolchain versions to their host distro.
#
#   docker build -t libirrep-dev .
#   docker run --rm -v "$(pwd):/w" -w /w libirrep-dev make all
#   docker run --rm -v "$(pwd):/w" -w /w libirrep-dev make asan
#   docker run --rm -v "$(pwd):/w" -w /w libirrep-dev make fuzz
#
# For the aarch64 CI matrix entry, cross-arch the image:
#   docker build --platform linux/arm64 -t libirrep-dev-arm64 .
#   docker run --platform linux/arm64 --rm -v "$(pwd):/w" -w /w libirrep-dev-arm64 make all
#
# Apple Silicon users can exercise the x86_64 AVX2 path under Rosetta
# the same way:
#   docker build --platform linux/amd64 -t libirrep-dev-x64 .
#   docker run --platform linux/amd64 --rm -v "$(pwd):/w" -w /w libirrep-dev-x64 make all

FROM ubuntu:22.04

# Pin to the same Debian-apt packages CI uses. The `build-essential`
# meta-package gives us gcc / g++ / make / libc-dev. `clang` is the
# ASan / UBSan / libFuzzer toolchain. `cmake` / `ninja-build` cover the
# second build surface. `python3` is used by the SBOM generator and the
# fuzz driver. `valgrind` + `gdb` are optional debugging aids.
RUN apt-get update -qq \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq --no-install-recommends \
      build-essential \
      clang \
      cmake \
      ninja-build \
      doxygen \
      git \
      python3 \
      valgrind \
      gdb \
      ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Default shell is bash rather than /bin/sh so `set -euxo pipefail`
# works in interactive sessions too.
SHELL ["/bin/bash", "-c"]

# Non-root build user — matches what GH Actions runners provide and
# avoids file-ownership issues with bind-mounted working directories.
ARG USER_ID=1000
ARG GROUP_ID=1000
RUN groupadd --gid $GROUP_ID dev 2>/dev/null || true \
 && useradd --uid $USER_ID --gid $GROUP_ID --shell /bin/bash --create-home dev
USER dev
WORKDIR /w

CMD ["make", "all"]
