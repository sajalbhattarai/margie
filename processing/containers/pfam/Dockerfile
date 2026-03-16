FROM ubuntu:22.04

# Metadata
LABEL maintainer="Bioinformatics Pipeline"
LABEL description="Pfam protein family annotation with HMMER 3.4"
LABEL version="1.0"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    tar \
    gzip \
    make \
    gcc \
    g++ \
    perl \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Download and install HMMER 3.4
RUN cd /tmp && \
    wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz && \
    tar -xzf hmmer-3.4.tar.gz && \
    cd hmmer-3.4 && \
    ./configure && \
    make -j$(nproc) && \
    make install && \
    cd / && \
    rm -rf /tmp/hmmer-3.4*

# Create working directories
RUN mkdir -p /container/scripts \
             /container/db \
             /container/input \
             /container/output

# Copy container scripts
COPY scripts/run_hmmscan.sh /container/scripts/
COPY scripts/consolidate_pfam.py /container/scripts/

# Make scripts executable
RUN chmod +x /container/scripts/*.sh

# Set working directory
WORKDIR /container

# Entry point
CMD ["/bin/bash"]
