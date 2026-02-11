FROM rocker/r2u:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Rome

# Reduce memory spikes from BLAS/OMP
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1

WORKDIR /workspace

# -----------------------------
# System deps + R binary packages (r2u)
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-venv python3-pip \
    git curl wget ca-certificates locales \
    \
    # HDF5 + compression + ATAC utils
    libhdf5-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    tabix samtools gzip \
    \
    # Rmd HTML
    pandoc \
    \
    # devtools deps
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libcairo2-dev \
    \
    # -------------------------
    # R packages as Ubuntu binaries (no compile)
    # -------------------------
    r-cran-seurat \
    r-cran-signac \
    r-cran-matrix \
    r-cran-data.table \
    r-cran-dplyr \
    r-cran-ggplot2 \
    r-cran-patchwork \
    r-cran-future \
    r-cran-hdf5r \
    r-cran-rmarkdown \
    r-cran-knitr \
    r-cran-irkernel \
    r-cran-roxygen2 \
    r-cran-usethis \
    r-cran-devtools \
    r-cran-pkgdown \
    r-cran-biocmanager \
    r-cran-mclust \
    \
    # Bioconductor binaries
    r-bioc-dropletutils \
    r-bioc-singlecellexperiment \
    r-bioc-summarizedexperiment \
    r-bioc-biocparallel \
    r-bioc-delayedarray \
    r-bioc-s4vectors \
    r-bioc-biocgenerics \
    \
    # Genomics deps
    r-bioc-genomicranges \
    r-bioc-genomeinfodb \
    r-bioc-rtracklayer \
    r-bioc-rsamtools \
    r-bioc-genomicalignments \
    && rm -rf /var/lib/apt/lists/*

# Locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# -----------------------------
# JupyterLab in venv
# -----------------------------
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
 && /opt/venv/bin/pip install --no-cache-dir jupyterlab
ENV PATH="/opt/venv/bin:${PATH}"

# Register R kernel for Jupyter
RUN R -q -e "IRkernel::installspec(user=FALSE)"

EXPOSE 8888

CMD ["bash","-lc", "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password=''"]
