# build minimap2-coverage
FROM continuumio/miniconda3
### MAINTAINER ###
MAINTAINER Yoshinori Fukasawa <yoshinori.fukasawa@kaust.edu.sa>

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

### ENV VARS ###
ENV USER user
ENV HOME /home/${USER}

### LABELS ###
LABEL base_image="miniconda3"
LABEL version="1.1"
LABEL software="LongQC docker"
LABEL software.version="1.1"

# add a general user account
RUN useradd -m ${USER}
# define a password for user
RUN echo "${USER}:test_pass" | chpasswd

ADD https://api.github.com/repos/yfukasawa/longqc/git/refs/heads/master version.json
RUN git clone https://github.com/yfukasawa/LongQC.git $HOME/LongQC
RUN cd $HOME/LongQC/minimap2_mod && make extra

# install dependency
RUN conda update -y conda

RUN conda install -y numpy
RUN conda install -y pandas'>=0.24.0'
RUN conda install -y scipy
RUN conda install -y jinja2
RUN conda install -y h5py
RUN conda install -y matplotlib'>=2.1.2'
RUN conda install -y scikit-learn

RUN conda install -y -c bioconda pysam
RUN conda install -y -c bioconda edlib
RUN conda install -y -c bioconda python-edlib

# change user to "user" defined above
USER ${USER}

# define a working dir
WORKDIR $HOME
