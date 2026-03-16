# Margie GUI Container
# Provides a self-contained Streamlit web interface for the Margie annotation pipeline.
# The repo is mounted at runtime using Docker volumes; no pipeline data is baked in.
#
# Build:
#   docker build --platform linux/amd64 -t margie-gui:1.0 .
#
# Run:
#   docker run --rm -it -p 8501:8501 \
#       -v /path/to/margie:/margie:ro \
#       -v /path/to/margie/input:/margie/input:rw \
#       -v /path/to/margie/output:/margie/output:rw \
#       margie-gui:1.0
#   Then open http://localhost:8501

FROM python:3.11-slim

LABEL maintainer="Margie Pipeline"
LABEL description="Streamlit GUI for the Margie genome annotation pipeline"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_ENABLE_CORS=false \
    STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION=false

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies for the GUI
RUN pip install --no-cache-dir \
    "streamlit>=1.33.0" \
    "pandas>=2.0.0" \
    "altair>=5.0.0" \
    "watchdog"

# The Margie repo will be bind-mounted at /margie at runtime.
# We create the expected directory structure so Streamlit can start even if
# the volume has not been mounted yet (shows a helpful message instead of crashing).
RUN mkdir -p /margie/app/gui /margie/input /margie/output /margie/db

# Lightweight placeholder so the app can be imported without the full repo
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

WORKDIR /margie

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

ENTRYPOINT ["/entrypoint.sh"]
