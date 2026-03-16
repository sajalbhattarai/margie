#!/bin/bash
# Pull the official eggNOG-mapper biocontainer used by Margie setup.

set -euo pipefail

IMAGE_NAME="quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"

echo "=========================================="
echo "Pulling official EggNOG container"
echo "=========================================="
echo "Image: ${IMAGE_NAME}"
echo ""

if ! command -v docker >/dev/null 2>&1; then
	echo "ERROR: Docker is not installed or not in PATH"
	exit 1
fi

docker pull --platform linux/amd64 "$IMAGE_NAME"

echo ""
echo "=========================================="
echo "✓ Official image ready"
echo "=========================================="
docker images "$IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
echo ""
echo "Margie setup uses this official eggNOG-mapper image for runtime setup."
