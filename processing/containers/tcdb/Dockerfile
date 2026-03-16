# TCDB (Transporter Classification Database) Annotation Container
# Supports official BLAST+ workflow and optional DIAMOND acceleration

FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Create working directories
RUN mkdir -p /container/input \
    /container/output \
    /container/db \
    /container/scripts

# Copy scripts
COPY scripts/ /container/scripts/

# Make scripts executable
RUN find /container/scripts -maxdepth 1 -type f -name "*.sh" -exec chmod +x {} \;

# Set working directory
WORKDIR /container

# Default command
CMD ["/bin/bash"]
