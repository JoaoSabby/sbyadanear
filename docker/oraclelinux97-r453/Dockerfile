FROM oraclelinux:9
ARG R_VERSION=4.5.3
SHELL ["/bin/bash", "-c"]
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 R_VERSION=${R_VERSION}
RUN dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled ol9_codeready_builder && \
    dnf -y install gcc gcc-c++ gcc-gfortran make which tar gzip bzip2 xz \
      curl ca-certificates diffutils findutils perl python3 git \
      readline-devel zlib-devel bzip2-devel xz-devel pcre2-devel \
      libcurl-devel openssl-devel libxml2-devel libicu-devel \
      libpng-devel libjpeg-turbo-devel libtiff-devel cairo-devel \
      fontconfig-devel freetype-devel harfbuzz-devel fribidi-devel && \
    dnf clean all
RUN printf '%s\n' '[oneAPI]' 'name=Intel oneAPI' \
      'baseurl=https://yum.repos.intel.com/oneapi' 'enabled=1' \
      'gpgcheck=1' 'repo_gpgcheck=1' \
      'gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB' \
      > /etc/yum.repos.d/oneAPI.repo && \
    dnf -y install intel-oneapi-compiler-dpcpp-cpp intel-oneapi-mkl-devel && \
    dnf clean all
# R uses the same Intel BLAS. Disable its optional OpenMP to avoid introducing
# libgomp through the core or dependency defaults. The package enables Intel
# OpenMP explicitly in its Makevars.
RUN source /opt/intel/oneapi/setvars.sh && \
    mkl_lib=$(find -L "$MKLROOT" -name libmkl_core.so -printf '%h\n' -quit) && \
    omp_lib=$(find -L /opt/intel/oneapi/compiler -name libiomp5.so -printf '%h\n' -quit) && \
    test -n "$mkl_lib" && test -n "$omp_lib" && \
    printf '%s\n' "$mkl_lib" "$omp_lib" > /etc/ld.so.conf.d/intel-oneapi.conf && ldconfig && \
    curl -fsSL "https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz" -o /tmp/R.tar.gz && \
    tar -xzf /tmp/R.tar.gz -C /tmp && cd "/tmp/R-${R_VERSION}" && \
    CC=icx CXX=icpx FC=gfortran F77=gfortran ./configure \
      --prefix=/usr/local --enable-R-shlib --without-x \
      --disable-openmp \
      --with-blas="-L${mkl_lib} -Wl,-rpath,${mkl_lib} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl" \
      --with-lapack && make -j2 && make install && \
    cd / && rm -rf "/tmp/R-${R_VERSION}" /tmp/R.tar.gz
RUN source /opt/intel/oneapi/setvars.sh && \
    Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages(c("Rcpp","tibble","recipes","rlang","generics","testthat","roxygen2"), Ncpus=2L)'
COPY docker/oneapi-entrypoint.sh /usr/local/bin/oneapi-entrypoint
RUN chmod +x /usr/local/bin/oneapi-entrypoint
WORKDIR /workspace/sbyadanear
ENTRYPOINT ["/usr/local/bin/oneapi-entrypoint"]
CMD ["bash"]
