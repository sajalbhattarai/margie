#!/bin/bash
# Build COG annotation container

set -e

CONTAINER_NAME="cog-annotation"
VERSION="1.0"

echo "=========================================="
echo "Building COG Annotation Container"
echo "=========================================="
echo "Container: ${CONTAINER_NAME}:${VERSION}"
echo ""

# Build the container
docker build --platform linux/amd64 -t "${CONTAINER_NAME}:${VERSION}" .

echo ""
echo "=========================================="
echo "✓ Container built successfully!"
echo "=========================================="
echo ""
echo "Container: ${CONTAINER_NAME}:${VERSION}"
echo ""

# Show container size
docker images "${CONTAINER_NAME}:${VERSION}" --format "Size: {{.Size}}"
echo ""
echo "Installed tools:"
docker run --rm "${CONTAINER_NAME}:${VERSION}" COGclassifier --version
docker run --rm "${CONTAINER_NAME}:${VERSION}" rpsblast -version | head -n 1
echo ""
echo "You can now use this container for COG annotation."
echo "Required database directory: cddid.tbl.gz, Cog_LE.tar.gz, cog-24.def.tab, cog-24.fun.tab"
