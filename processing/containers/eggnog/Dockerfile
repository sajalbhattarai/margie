# EggNOG (evolutionary genealogy of genes) Annotation Container
# Uses Diamond for sequence search and eggNOG databases for annotation

FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz \
    && mv diamond /usr/local/bin/ \
    && rm diamond-linux64.tar.gz

# Create working directories
RUN mkdir -p /container/input \
    /container/output \
    /container/db \
    /container/scripts

# Copy scripts
COPY scripts/ /container/scripts/

# Make scripts executable
RUN chmod +x /container/scripts/*.sh

# Set working directory
WORKDIR /container

# Default command
CMD ["/bin/bash"]
