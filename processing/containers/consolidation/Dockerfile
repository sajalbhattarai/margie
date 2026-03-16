FROM python:3.11-slim

LABEL maintainer="annotation-pipeline"
LABEL description="Consolidation container for merging all genome annotations"

# Install pandas
RUN pip install --no-cache-dir pandas

# Copy consolidation script
COPY consolidate_all_annotations.py /usr/local/bin/consolidate_annotations
RUN chmod +x /usr/local/bin/consolidate_annotations

# Set working directory
WORKDIR /workspace

# Entry point
ENTRYPOINT ["python3", "/usr/local/bin/consolidate_annotations"]
