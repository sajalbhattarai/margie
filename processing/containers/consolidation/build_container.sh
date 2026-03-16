#!/bin/bash
#
# Build consolidation-annotation container
#

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

IMAGE_NAME="consolidation-annotation"
IMAGE_TAG="1.0"

echo "Building $IMAGE_NAME:$IMAGE_TAG container..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Build container
echo "Building Docker image..."
docker build \
    --platform linux/amd64,linux/arm64 \
    -t "${IMAGE_NAME}:${IMAGE_TAG}" \
    -t "${IMAGE_NAME}:latest" \
    .

echo ""
echo "✓ Container built successfully!"
echo "  Image: ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""
echo "To test:"
echo "  docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} --help"
