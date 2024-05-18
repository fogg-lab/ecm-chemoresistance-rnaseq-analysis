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
    git \
    libglpk40

# Install RangerBasediRF
RUN git clone https://github.com/jailGroup/RangerBasediRF.git
RUN cd RangerBasediRF/cpp_version && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_CXX_COMPILER=mpiCC -DCMAKE_C_COMPILER=mpicc .. && \
    make
RUN ln -s /RangerBasediRF/cpp_version/build/ranger /usr/local/bin/ranger

RUN install2.r --error --deps TRUE renv tidyverse here httpgd ggpubr parallelly styler plotly htmlwidgets

RUN R -e "BiocManager::install(c('DESeq2', 'apeglm'))" || exit 1
RUN R -e "BiocManager::install(c('edgeR', 'limma', 'sva'))" || exit 1
RUN R -e "BiocManager::install(c('TCGAbiolinks', 'clusterProfiler', 'org.Hs.eg.db'))" || exit 1
RUN R -e "BiocManager::install(c('fgsea', 'graph', 'msigdbr'))" || exit 1
RUN R -e "BiocManager::install(c('ReactomePA', 'rrvgo', 'enrichplot'))" || exit 1

RUN install2.r --error --deps TRUE pheatmap treemap

# Set R_LIBS for the current command to ensure dependencies are found
RUN git clone https://github.com/fogg-lab/ecm-chemoresistance-rnaseq-analysis.git && \
    cd ecm-chemoresistance-rnaseq-analysis/src/R && \
    R CMD build src && \
    R CMD INSTALL src_0.1.0.tar.gz && \
    pip install -e ../python

RUN pip install --upgrade pip wheel setuptools
COPY requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt
