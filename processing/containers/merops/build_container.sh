#!/bin/bash
# Build MEROPS annotation container

set -e

IMAGE_NAME="merops-annotation"
IMAGE_TAG="1.0"

echo "========================================"
echo "Building MEROPS Annotation Container"
echo "========================================"
echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""

if [ ! -f "Dockerfile" ]; then
    echo "Error: Dockerfile not found"
    exit 1
fi

echo "Building container..."
docker build --platform linux/amd64 -t "${IMAGE_NAME}:${IMAGE_TAG}" .

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✓ Container built successfully!"
    echo "========================================"
    echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
    
    echo ""
    echo "Image details:"
    docker images "${IMAGE_NAME}:${IMAGE_TAG}" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    
    echo ""
    echo "Test the container:"
    echo "  docker run --rm ${IMAGE_NAME}:${IMAGE_TAG}"
else
    echo ""
    echo "Error: Container build failed"
    exit 1
fi
