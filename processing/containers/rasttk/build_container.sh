#!/bin/bash
# Build BV-BRC RASTtk Container
# Independent container with BV-BRC tools for gene calling and functional annotation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE_NAME="rasttk-annotation"
IMAGE_TAG="1.0"

echo "=========================================="
echo "Building BV-BRC RASTtk Container"
echo "=========================================="
echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
echo "Platform: linux/amd64"
echo ""

cd "$SCRIPT_DIR"

# Build the container
docker build --platform linux/amd64 -t ${IMAGE_NAME}:${IMAGE_TAG} .

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Container built successfully!"
    echo ""
    echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
    echo "Size: $(docker images ${IMAGE_NAME}:${IMAGE_TAG} --format '{{.Size}}')"
    echo ""
    echo "To test the container:"
    echo "  docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} rast-create-genome --help"
else
    echo "✗ Build failed"
    exit 1
fi
