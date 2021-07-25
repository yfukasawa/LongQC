# build minimap2-coverage
FROM continuumio/miniconda3

### MAINTAINER ###
MAINTAINER Yoshinori Fukasawa <yoshinori.fukasawa@kaust.edu.sa>
### LABELS ###
LABEL base_image="miniconda3"
LABEL software="LongQC docker"
LABEL software.version="1.2"

### Basic dependency ###
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
    git \
    build-essential \
    libc6-dev \
    zlib1g-dev && \
    apt-get clean && \
    apt-get purge

### LongQC installation ###
ADD https://api.github.com/repos/yfukasawa/longqc/git/refs/heads/minimap2_update version.json
RUN git clone https://github.com/yfukasawa/LongQC.git $HOME/LongQC
RUN cd $HOME/LongQC/minimap2-coverage && make
RUN cd $HOME/LongQC && \
    sed -i \
    -e '1{s;^;#!/opt/conda/bin/python\n;}' \
    longQC.py && \
    chmod +x longQC.py

### install LongQC's dependency ###
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -y \
    python=3.9 \
    numpy \
    pandas'>=0.24.0' \
    scipy \
    jinja2 \
    h5py \
    matplotlib'>=2.1.2' \
    scikit-learn && \
    conda install -y -c bioconda \
    edlib \
    pysam \
    python-edlib

### Define PATH ###
ENV PATH="/root/LongQC:$PATH"

### Entry point
ENTRYPOINT ["longQC.py"]