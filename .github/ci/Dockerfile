FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --yes --no-install-recommends \
    ca-certificates \
    cmake \
    gfortran \
    libopenmpi-dev \
    ninja-build \
    openmpi-bin \
    python3 \
    python3-pip \
  && pip3 install --break-system-packages --no-cache-dir fypp \
  && rm -rf /var/lib/apt/lists/*
