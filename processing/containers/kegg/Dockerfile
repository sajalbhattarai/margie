# KEGG/KOfam Annotation Container
# Uses KofamScan tool (Ruby wrapper) with HMMER for KEGG Orthology assignment
# Lean architecture: databases mounted externally from host db/kegg/
# Platform specified as amd64 for HMMER SSE support
FROM --platform=linux/amd64 ubuntu:22.04

LABEL maintainer="Bioinformatics Pipeline"
LABEL description="KEGG annotation using KofamScan and HMMER"
LABEL version="1.0"

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies: Ruby, Python3, parallel, wget, basic build tools for HMMER
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    ruby \
    python3 \
    python3-pip \
    parallel \
    build-essential \
    gzip \
    tar \
    && rm -rf /var/lib/apt/lists/*

# Install HMMER 3.3.2 (required by KofamScan)
RUN cd /tmp && \
    wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz && \
    tar -xzf hmmer-3.3.2.tar.gz && \
    cd hmmer-3.3.2 && \
    ./configure --prefix=/usr/local && \
    make -j4 && \
    make install && \
    cd / && \
    rm -rf /tmp/hmmer-*

# Verify HMMER installation
RUN hmmsearch -h > /dev/null && echo "HMMER installed successfully"

# Download and install official KOfamScan from Genome.jp INSTALL path
RUN set -eux; \
    cd /tmp; \
    KOFAM_SCAN_TARBALL="kofam_scan-1.3.0.tar.gz"; \
    wget -O "$KOFAM_SCAN_TARBALL" "https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz"; \
    rm -rf /tmp/kofam_scan_extract /opt/kofam_scan; \
    mkdir -p /tmp/kofam_scan_extract /opt/kofam_scan; \
    tar -xzf "$KOFAM_SCAN_TARBALL" -C /tmp/kofam_scan_extract; \
    if [ -f /tmp/kofam_scan_extract/exec_annotation ]; then \
        cp -a /tmp/kofam_scan_extract/. /opt/kofam_scan/; \
    else \
        EXTRACT_TOP="$(find /tmp/kofam_scan_extract -mindepth 1 -maxdepth 1 -type d | head -n 1)"; \
        test -n "$EXTRACT_TOP"; \
        cp -a "$EXTRACT_TOP"/. /opt/kofam_scan/; \
    fi; \
    test -x /opt/kofam_scan/exec_annotation; \
    ln -sf /opt/kofam_scan/exec_annotation /usr/local/bin/exec_annotation; \
    rm -rf /tmp/kofam_scan_extract /tmp/kofam_scan* /tmp/kofamscan*

# Create container directory structure
# Database will be mounted from host at /container/db/
RUN mkdir -p /container/scripts \
             /container/db \
             /container/input \
             /container/output

# Copy processing scripts
COPY scripts/*.sh /container/scripts/
COPY scripts/*.py /container/scripts/

# Make scripts executable
RUN find /container/scripts -maxdepth 1 -type f -name "*.sh" -exec chmod +x {} \;

# Set working directory
WORKDIR /container

# Default command
CMD ["/bin/bash"]
