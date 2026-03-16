#!/bin/bash
# Build UniProt annotation container

set -euo pipefail

CONTAINER_NAME="uniprot-annotation"
VERSION="1.0"

echo "Building ${CONTAINER_NAME}:${VERSION} container..."
echo ""

# Build the container for Apple Silicon / Linux x86_64
docker build --platform linux/amd64 -t "${CONTAINER_NAME}:${VERSION}" .

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Container built successfully!"
    echo ""
    echo "Container: ${CONTAINER_NAME}:${VERSION}"
    echo ""
    
    # Show container size
    SIZE=$(docker images "${CONTAINER_NAME}:${VERSION}" --format "{{.Size}}")
    echo "Size: $SIZE"
    echo ""
    
    # Show installed tool versions
    echo "Installed tools:"
    docker run --rm "${CONTAINER_NAME}:${VERSION}" upimapi --version
    docker run --rm "${CONTAINER_NAME}:${VERSION}" diamond --version
    docker run --rm "${CONTAINER_NAME}:${VERSION}" python --version
else
    echo ""
    echo "✗ Container build failed!"
    exit 1
fi
