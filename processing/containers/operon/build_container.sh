#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="operon-annotation"
IMAGE_TAG="1.0"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Building ${IMAGE_NAME}:${IMAGE_TAG}"
docker build --platform linux/amd64 -t "${IMAGE_NAME}:${IMAGE_TAG}" "${SCRIPT_DIR}"

echo "Build complete"
echo "Run test command:"
echo "  docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} bash -lc 'UniOP --help || true'"
