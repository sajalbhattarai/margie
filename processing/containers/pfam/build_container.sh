#!/bin/bash
# Build script for Pfam annotation container

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "========================================"
echo "Building Pfam Annotation Container"
echo "========================================"
echo "Directory: $SCRIPT_DIR"
echo ""

# Check prerequisites
echo "Checking prerequisites..."
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is not installed or not in PATH"
    exit 1
fi
echo "✓ Docker found"

# Build container
echo ""
echo "Building container..."
docker build --platform linux/amd64 -t pfam-annotation:1.0 .

# Verify build
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✓ Container built successfully!"
    echo "========================================"
    echo ""
    echo "Container information:"
    docker images pfam-annotation:1.0 --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    echo ""
    echo "Test the container:"
    echo "  docker run --rm pfam-annotation:1.0 hmmscan -h"
    echo ""
else
    echo ""
    echo "========================================"
    echo "ERROR: Container build failed"
    echo "========================================"
    exit 1
fi
