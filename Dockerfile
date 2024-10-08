# Image for building EXP with MPI and CUDA based on the hpcbase recipe
#    and our standard Ubuntu 22.04 development environment for Docker
# 
# Contents:
# 
#   * Ubuntu 22.04 
#   * GNU compilers (upstream)
#   * Python 2 and 3 (upstream)
#   * OpenMPI, FFTW3, HDF5, Eigen3, PNG (upstream)
# 
# Build notes:
# 
#   * The container build requires the following commands:
#     $ hpccm --recipe exp_all_deb.py --format docker > Dockerfile
#     $ docker build -t exp-test -f Dockerfile .
#   or
#     $ docker buildx build --platform=linux/amd64,linux/arm64 -t the9cat/exp:latest --push -f Dockerfile .
# 
#   * You will need to put the EXP.tar.gz file in the build directory. I
#     like to make a fresh clone and run `git submodule update --init
#     --recursive` before tarring it up.
# 
# 

FROM ubuntu:22.04 AS devel

# Python
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python2 \
        python3 && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        cmake \
        git \
        libeigen3-dev \
        libfftw3-dev \
        libhdf5-103 \
        libhdf5-dev \
        libopenmpi-dev \
        libpng-dev \
        make \
        openmpi-bin \
        python3-dev \
        tar \
        wget && \
    rm -rf /var/lib/apt/lists/*

# git@github.com:EXP-code/EXP.git
RUN mkdir -p /var/tmp && cd /var/tmp && git clone --depth=1 git://github.com/EXP-code/EXP.git EXP && cd - && \
    cd /var/tmp/EXP && \
    git config --global --add safe.directory /var/tmp/EXP && \
    git config --global --add safe.directory /var/tmp/EXP/extern/HighFive && \
    git config --global --add safe.directory /var/tmp/EXP/extern/pybind11 && \
    git config --global --add safe.directory /var/tmp/EXP/extern/yaml-cpp && \
    git config --global --add safe.directory /var/tmp/EXP/extern/HighFive/deps/catch2 && \
    mkdir -p /usr/local/EXP/doc && \
    cp -a /var/tmp/EXP/sphinx/* /usr/local/EXP/doc && \
    mkdir -p /var/tmp/EXP/build && cd /var/tmp/EXP/build && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/EXP -D CMAKE_BUILD_TYPE=Release -D ENABLE_CUDA=NO -D ENABLE_USER=YES -D ENABLE_SLURM=NO -D ENABLE_PNG=NO -D ENABLE_VTK=NO -D FFTW_INCLUDE_DIRS=/usr/include/fftw3 -D Eigen3_DIR=/usr/share/eigen3/cmake -D CMAKE_INSTALL_PREFIX=/usr/local/EXP /var/tmp/EXP && \
    cmake --build /var/tmp/EXP/build --target all -- -j$(nproc) && \
    cmake --build /var/tmp/EXP/build --target install -- -j$(nproc) && \
    rm -rf /var/tmp/EXP

FROM ubuntu:22.04

# Python
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python2 \
        python3 && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler runtime
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libgfortran5 \
        libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# git@github.com:EXP-code/EXP.git
COPY --from=devel /usr/local/EXP /usr/local/EXP
ENV LD_LIBRARY_PATH=/usr/local/EXP/lib \
    LIBRARY_PATH=/usr/local/EXP/lib \
    PATH=/usr/local/EXP/bin:${PATH}

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        busybox \
        dvipng \
        ffmpeg \
        git \
        less \
        libeigen3-dev \
        libfftw3-3 \
        libgsl-dev \
        libhdf5-103 \
        libhdf5-cpp-103 \
        libhdf5-dev \
        libopenmpi-dev \
        libpython3.10-dev \
        make \
        nano \
        openmpi-bin \
        python3.10-dev \
        unzip && \
    rm -rf /var/lib/apt/lists/*

COPY --from=devel /usr/local/EXP/bin /usr/local/EXP/bin

COPY --from=devel /usr/local/EXP/lib /usr/local/EXP/lib

COPY --from=devel /usr/local/EXP/doc /var/www/html

ENV PATH=/usr/local/EXP/bin:$PATH

ENV LIBRARY_PATH=/usr/local/EXP/lib

ENV LD_LIBRARY_PATH=/usr/local/EXP/lib

ENV PYTHONPATH=/usr/local/EXP/lib/python3.10/site-packages

# pip
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3-pip \
        python3-pip-whl \
        python3-setuptools \
        python3-wheel && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 --no-cache-dir install --upgrade pip && \
    pip3 --no-cache-dir install numpy scipy matplotlib jupyter h5py mpi4py PyYAML k3d pandas astropy gala galpy pynbody jupyterlab ipyparallel

# Add a user with a home directory

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
    
# Make sure the contents of our repo are in ${HOME}
# COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
