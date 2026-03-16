FROM python:3.11-slim

LABEL maintainer="Margie Annotation Pipeline"
LABEL description="Operon prediction container based on UniOP source build"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV UNIOP_HOME=/opt/UniOP

RUN apt-get update && apt-get install -y \
    bash \
    ca-certificates \
    git \
    libgomp1 \
    prodigal \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/hongsua/UniOP.git "$UNIOP_HOME" && \
    python -m pip install --no-cache-dir numpy pandas scikit-learn && \
    printf '#!/usr/bin/env bash\nexec python3 /opt/UniOP/src/UniOP "$@"\n' > /usr/local/bin/UniOP && \
    chmod +x /usr/local/bin/UniOP

RUN mkdir -p /container/scripts /container/input /container/output

COPY scripts/ /container/scripts/

RUN chmod +x /container/scripts/*.sh || true

WORKDIR /container

CMD ["/bin/bash"]
