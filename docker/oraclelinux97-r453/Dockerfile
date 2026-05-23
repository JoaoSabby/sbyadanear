FROM container-registry.oracle.com/os/oraclelinux:9

ARG R_VERSION=4.5.3
ARG CRAN=https://cloud.r-project.org
ARG SBY_USE_ONEAPI_MKL=0

ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=America/Sao_Paulo \
    CRAN=${CRAN} \
    R_VERSION=${R_VERSION} \
    R_REMOTES_NO_ERRORS_FROM_WARNINGS=true \
    R_KEEP_PKG_SOURCE=yes \
    MAKEFLAGS=-j2 \
    SBY_USE_ONEAPI_MKL=${SBY_USE_ONEAPI_MKL} \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    BLIS_NUM_THREADS=1

RUN cat /etc/oracle-release

RUN dnf -y update && \
    dnf -y install \
      which \
      tar \
      gzip \
      xz \
      bzip2 \
      unzip \
      zip \
      git \
      make \
      gcc \
      gcc-c++ \
      gcc-gfortran \
      patch \
      diffutils \
      findutils \
      procps-ng \
      hostname \
      file \
      curl \
      wget \
      rsync \
      openssl-devel \
      libcurl-devel \
      libxml2-devel \
      zlib-devel \
      bzip2-devel \
      xz-devel \
      pcre2-devel \
      readline-devel \
      libjpeg-turbo-devel \
      libpng-devel \
      libtiff-devel \
      cairo-devel \
      freetype-devel \
      fontconfig-devel \
      harfbuzz-devel \
      fribidi-devel \
      libicu-devel \
      glibc-langpack-en \
      glibc-langpack-pt && \
    dnf clean all && \
    rm -rf /var/cache/dnf

RUN curl -fsSL "https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz" \
      -o "/tmp/R-${R_VERSION}.tar.gz" && \
    tar -xzf "/tmp/R-${R_VERSION}.tar.gz" -C /tmp && \
    cd "/tmp/R-${R_VERSION}" && \
    ./configure \
      --prefix="/opt/R/${R_VERSION}" \
      --enable-R-shlib \
      --with-blas \
      --with-lapack \
      --with-readline=yes \
      --with-x=no && \
    make -j"$(nproc)" && \
    make install && \
    ln -s "/opt/R/${R_VERSION}/bin/R" /usr/local/bin/R && \
    ln -s "/opt/R/${R_VERSION}/bin/Rscript" /usr/local/bin/Rscript && \
    rm -rf "/tmp/R-${R_VERSION}" "/tmp/R-${R_VERSION}.tar.gz"

RUN Rscript -e 'install.packages(c("cli", "FNN", "generics", "recipes", "rlang", "RcppHNSW", "Rfast", "tibble", "glue", "stringr", "testthat", "knitr", "rmarkdown", "withr", "callr", "bench", "microbenchmark", "BiocManager"), repos = Sys.getenv("CRAN"), Ncpus = max(1L, parallel::detectCores() - 1L))'

RUN Rscript -e 'BiocManager::install(c("BiocNeighbors", "BiocParallel"), ask = FALSE, update = FALSE)'

RUN R --version

WORKDIR /workspace/sbyadanear

COPY . /workspace/sbyadanear

RUN SBY_USE_ONEAPI_MKL=${SBY_USE_ONEAPI_MKL} R CMD INSTALL --preclean .

CMD ["bash"]
