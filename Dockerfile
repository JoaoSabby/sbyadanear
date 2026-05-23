FROM rocker/r-ver:4.4.3

ENV DEBIAN_FRONTEND=noninteractive \
    R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    qpdf \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("RcppHNSW", "cli", "FNN", "generics", "recipes", "rlang", "Rfast", "tibble"), repos = "https://cloud.r-project.org")'

WORKDIR /workspace/sbyadanear

CMD ["bash"]
