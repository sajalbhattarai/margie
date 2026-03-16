# TIGRfam Protein Family Annotation Container
# Uses HMMER for TIGRfam annotation
# Expected database mounts:
#   - /container/db/TIGRFAMs_15.0_HMM.LIB (main HMM database)
#   - /container/db/tigrfam_definitions.txt (family info)

FROM ubuntu:22.04

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
    bc \
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
RUN mkdir -p /container/input /container/output /container/db /container/scripts

# Copy scripts
COPY scripts/ /container/scripts/

# Make scripts executable
RUN chmod +x /container/scripts/*.sh

# Set working directory
WORKDIR /container

# Default command shows usage
CMD ["cat", "/container/scripts/README.txt"]
