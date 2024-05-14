FROM ghcr.io/rocker-org/devcontainer/r-ver:4.4

RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libxtst6 \
    libxt6 \
    libmagick++-dev \
    python3-pip \
    python-is-python3 \
    libcurl4-openssl-dev \
    cmake \
    mpich \
    git

# Install RangerBasediRF
RUN git clone https://github.com/jailGroup/RangerBasediRF.git
RUN cd RangerBasediRF/cpp_version && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_CXX_COMPILER=mpiCC -DCMAKE_C_COMPILER=mpicc .. && \
    make
RUN ln -s /RangerBasediRF/cpp_version/build/ranger /usr/local/bin/ranger

COPY requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt

RUN install2.r --error --deps TRUE renv tidyverse here httpgd

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'TCGAbiolinks', 'sva'))"
