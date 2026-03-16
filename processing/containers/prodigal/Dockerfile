FROM ubuntu:22.04

LABEL maintainer="Context-Based Annotation Pipeline"
LABEL description="Prodigal gene-calling container"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV APT_LISTCHANGES_FRONTEND=none

RUN apt-get -o Acquire::ForceIPv4=true -o Acquire::Retries=5 update && apt-get install -y --no-install-recommends \
    prodigal \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /container/input /container/output
WORKDIR /container

CMD ["/bin/bash"]
