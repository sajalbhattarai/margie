#!/bin/bash
# Build TCDB annotation container

set -e

CONTAINER_NAME="tcdb-annotation"
VERSION="1.0"

echo "=========================================="
echo "Building TCDB Annotation Container"
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
echo "You can now use this container for TCDB annotation."
echo "Required databases: tcdb_blast.* (official BLAST), tcdb_families.tsv, tcdb_substrates.tsv"
