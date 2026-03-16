#!/bin/bash
# Build script for KEGG annotation container

set -e

IMAGE_NAME="kegg-annotation:1.0"
PLATFORM="linux/amd64"

echo "========================================"
echo "  Building KEGG Annotation Container   "
echo "========================================"
echo "Image: $IMAGE_NAME"
echo "Platform: $PLATFORM"
echo "========================================"
echo ""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Build the container
echo "Building container..."
docker build --platform "$PLATFORM" -t "$IMAGE_NAME" .

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✓ Build successful!"
    echo "========================================"
    echo ""
    
    # Show image size
    echo "Image details:"
    docker images "$IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    echo ""
    
    # Verify HMMER and KofamScan
    echo "Verifying tools..."
    docker run --rm "$IMAGE_NAME" hmmsearch -h > /dev/null && echo "✓ HMMER installed"
    docker run --rm "$IMAGE_NAME" exec_annotation --help > /dev/null 2>&1 && echo "✓ KofamScan installed"
    
    echo ""
    echo "Container ready for use!"
    echo "Test with: ./run_container.sh"
else
    echo ""
    echo "========================================"
    echo "✗ Build failed"
    echo "========================================"
    exit 1
fi
