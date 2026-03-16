FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

RUN apt-get update && \
    apt-get install -y \
    bash \
    ca-certificates \
    coreutils \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y -c conda-forge -c bioconda upimapi=1.13.3 && \
    conda clean -afy

WORKDIR /container

COPY scripts/ /container/scripts/

RUN find /container/scripts -maxdepth 1 -type f -name "*.sh" -exec chmod +x {} \;

CMD ["/bin/bash"]
