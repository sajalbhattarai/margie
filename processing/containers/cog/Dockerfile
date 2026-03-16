# COG (Clusters of Orthologous Groups) Annotation Container
# Uses official COGclassifier workflow (RPS-BLAST + COG functional classification)
FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

RUN apt-get update && apt-get install -y \
    bash \
    ca-certificates \
    coreutils \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y -c conda-forge -c bioconda cogclassifier blast && \
    conda clean -afy

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
