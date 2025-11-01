#!/usr/bin/env bash
# Build inquiSTR as a statically linked MUSL binary
# Target: x86_64-unknown-linux-musl
# Usage:
#   scripts/build_musl.sh
# Requires:
#   - rustup target add x86_64-unknown-linux-musl
#   - Either `cross` installed (recommended) OR `musl-tools` (musl-gcc) on the host
set -euo pipefail

TARGET=x86_64-unknown-linux-musl
BIN=target/${TARGET}/release/inquiSTR

rustup target add "${TARGET}" >/dev/null 2>&1 || true

# Prefer `cross` for portable builds (uses Docker images with musl)
if command -v cross >/dev/null 2>&1; then
  echo "[build-musl] Using cross"
  OPENSSL_STATIC=1 \
  LIBZ_SYS_STATIC=1 \
  BZIP2_STATIC=1 \
  ZSTD_STATIC=1 \
  LZMA_API_STATIC=1 \
  CURL_STATIC=1 \
  cross build --release --target "${TARGET}"
else
  echo "[build-musl] Using cargo (host toolchain)"
  if ! command -v musl-gcc >/dev/null 2>&1; then
    echo "[build-musl] musl-gcc not found. On Ubuntu/Debian: sudo apt-get install -y musl-tools" >&2
    exit 1
  fi
  # Encourage static linking of common C libraries used by -sys crates
  OPENSSL_STATIC=1 \
  LIBZ_SYS_STATIC=1 \
  BZIP2_STATIC=1 \
  ZSTD_STATIC=1 \
  LZMA_API_STATIC=1 \
  CURL_STATIC=1 \
  cargo build --release --target "${TARGET}"
fi

if [ -f "${BIN}" ]; then
  echo "[build-musl] Built: ${BIN}"
  file "${BIN}" || true
else
  echo "[build-musl] Build did not produce ${BIN}" >&2
  exit 1
fi
