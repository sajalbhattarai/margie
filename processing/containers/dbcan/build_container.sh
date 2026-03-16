#!/bin/bash
# Pull the official dbCAN runner image used by Margie setup.

set -euo pipefail

IMAGE_NAME="haidyi/run_dbcan:latest"

echo "========================================"
echo "Pulling official dbCAN image"
echo "========================================"
echo "Image: $IMAGE_NAME"
echo ""

if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: Docker is not installed or not in PATH"
    exit 1
fi

docker pull --platform linux/amd64 "$IMAGE_NAME"

echo ""
echo "========================================"
echo "✓ Official image ready"
echo "========================================"
docker images "$IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
echo ""
echo "Margie setup uses this official image for dbCAN builds and database preparation."
