FROM python:3.11-slim-bookworm

# Dependencias del sistema (R + libs para compilar paquetes)
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    gfortran \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Instalar SoupX en R
RUN R -e "install.packages('SoupX', repos='https://cloud.r-project.org/')"

WORKDIR /app

# Dependencias Python (copiadas antes del código para cache de layers)
COPY requirements.txt requirements-scrna.txt ./
RUN pip install --no-cache-dir meson-python meson ninja cython numpy && \
    pip install --no-cache-dir -r requirements.txt -r requirements-scrna.txt

# Código de la app
COPY . .

# Variable para desactivar auto-shutdown en Docker
ENV BEGINSEQ_DOCKER=1
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')" || exit 1

CMD ["streamlit", "run", "app.py", \
     "--server.address=0.0.0.0", \
     "--server.port=8501", \
     "--server.headless=true"]
