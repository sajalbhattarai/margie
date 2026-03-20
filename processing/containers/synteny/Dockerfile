FROM --platform=linux/amd64 ubuntu:22.04

LABEL maintainer="Margie"
LABEL description="Synteny workflow container (MCScanX + DIAMOND + plotting)"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    ncbi-blast+ \
    python3 \
    python3-pip \
    tar \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install DIAMOND binary
RUN set -eux; \
    cd /tmp; \
    wget -O diamond.tar.gz "http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz"; \
    tar -xzf diamond.tar.gz; \
    mv diamond /usr/local/bin/diamond; \
    chmod +x /usr/local/bin/diamond; \
    rm -rf /tmp/diamond*

# Install MCScanX from official release
RUN set -eux; \
    cd /tmp; \
    wget -O mcscanx.tar.gz "https://github.com/wyp1125/MCScanX/archive/refs/tags/v1.0.0.tar.gz"; \
    tar -xzf mcscanx.tar.gz; \
    cd MCScanX-1.0.0; \
    g++ struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc msa.cc permutation.cc -o MCScanX; \
    cp MCScanX /usr/local/bin/MCScanX; \
    chmod +x /usr/local/bin/MCScanX; \
    rm -rf /tmp/mcscanx.tar.gz /tmp/MCScanX-1.0.0

RUN pip3 install --no-cache-dir pandas matplotlib

RUN mkdir -p /container/scripts /container/output

COPY scripts/ /container/scripts/
RUN chmod +x /container/scripts/*.sh /container/scripts/*.py

WORKDIR /container
CMD ["/bin/bash"]
