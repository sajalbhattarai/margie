# MEROPS Peptidase Annotation Container
# Uses Diamond for fast peptidase family identification

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install Diamond and Python
RUN apt-get update && \
    apt-get install -y \
        diamond-aligner \
        python3 \
        python3-pip \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create working directories
RUN mkdir -p /container/input /container/output /container/db /container/scripts

# Copy scripts
COPY scripts/ /container/scripts/

# Make scripts executable
RUN find /container/scripts -maxdepth 1 -type f -name "*.sh" -exec chmod +x {} \;

WORKDIR /container

CMD ["cat", "/container/scripts/README.txt"]
