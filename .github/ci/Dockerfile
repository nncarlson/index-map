FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV MPICH_VERSION=4.3.1
ENV OPENCOARRAYS_VERSION=2.10.3
ENV MPICH_PREFIX=/opt/mpich/${MPICH_VERSION}
ENV OPENCOARRAYS_PREFIX=/opt/opencoarrays/${OPENCOARRAYS_VERSION}
ENV PATH=${OPENCOARRAYS_PREFIX}/bin:${MPICH_PREFIX}/bin:${PATH}
ENV LD_LIBRARY_PATH=${OPENCOARRAYS_PREFIX}/lib64:${MPICH_PREFIX}/lib:${LD_LIBRARY_PATH:-}

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    gfortran \
    git \
    ninja-build \
    python3 \
    python3-pip \
    wget \
  && wget -q https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz -O /tmp/mpich.tar.gz \
  && tar -C /tmp -xf /tmp/mpich.tar.gz \
  && cd /tmp/mpich-${MPICH_VERSION} \
  && ./configure \
       FCFLAGS="-fallow-argument-mismatch" \
       FFLAGS="-fallow-argument-mismatch" \
       --prefix="${MPICH_PREFIX}" \
  && make -j"$(nproc)" FCFLAGS="-fallow-argument-mismatch" FFLAGS="-fallow-argument-mismatch" \
  && make install \
  && git clone --branch "${OPENCOARRAYS_VERSION}" --depth 1 https://github.com/sourceryinstitute/OpenCoarrays.git /tmp/OpenCoarrays \
  && cmake -S /tmp/OpenCoarrays -B /tmp/OpenCoarrays/build \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_C_COMPILER=gcc \
       -DCMAKE_Fortran_COMPILER=gfortran \
       -DCMAKE_INSTALL_PREFIX="${OPENCOARRAYS_PREFIX}" \
       -DMPI_C_COMPILER="${MPICH_PREFIX}/bin/mpicc" \
       -DMPI_Fortran_COMPILER="${MPICH_PREFIX}/bin/mpifort" \
  && cmake --build /tmp/OpenCoarrays/build --parallel \
  && cmake --install /tmp/OpenCoarrays/build \
  && cd / \
  && rm -rf /tmp/OpenCoarrays /tmp/mpich-${MPICH_VERSION} /tmp/mpich.tar.gz \
  && pip3 install --break-system-packages --no-cache-dir fypp \
  && rm -rf /var/lib/apt/lists/*
