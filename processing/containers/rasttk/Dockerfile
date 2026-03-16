FROM --platform=linux/amd64 ubuntu:22.04

# Container metadata
LABEL maintainer="BV-BRC RASTtk Container"
LABEL description="Independent BV-BRC RASTtk container for gene calling and functional annotation"
LABEL version="1.0"

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    python3 \
    python3-pip \
    git \
    less \
    vim \
    nano \
    ca-certificates \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create directory for installation
WORKDIR /tmp

# Download and install BV-BRC CLI (.deb package)
# Using version 1.040 for stability
RUN wget -q https://github.com/BV-BRC/BV-BRC-CLI/releases/download/1.040/bvbrc-cli-1.040.deb && \
    apt-get update && \
    apt-get install -y ./bvbrc-cli-1.040.deb && \
    rm bvbrc-cli-1.040.deb && \
    rm -rf /var/lib/apt/lists/*

# Verify BV-BRC installation
RUN which rast-create-genome || (echo "BV-BRC installation failed" && exit 1)

# Create container directories
WORKDIR /
RUN mkdir -p /container/scripts && \
    mkdir -p /container/db && \
    mkdir -p /container/input && \
    mkdir -p /container/output

# Copy pipeline scripts
COPY scripts/ /container/scripts/
RUN find /container/scripts -maxdepth 1 -type f -name "*.sh" -exec chmod +x {} \;

# Note: Database will be mounted from host at /container/db/
# This keeps the container size small and allows database updates without rebuild

# Set working directory for user data
WORKDIR /container

# Add scripts to PATH
ENV PATH="/container/scripts:${PATH}"

# Entrypoint
CMD ["/bin/bash"]
