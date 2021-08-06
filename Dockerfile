FROM docker.io/library/ubuntu:18.04

LABEL eu.ifom.bioinformatics.packages.docker.author="Raoul Jean Pierre Bonnal"
LABEL eu.ifom.bioinformatics.packages.docker.maintainer="raoul.bonnal@ifom.eu"
LABEL authors="Akiko Oguchi <aoguchi@kuhp.kyoto-u.ac.jp>,Raoul J.P. Bonnal <raoul.bonnal@ifom.eu>,Yasuhiro Murakawa <yasuhiro.murakawa@riken.jp>" \
      description="Docker image containing all software requirements for the nf-core/reaptec pipeline"

LABEL version="0.0.3"


ENV TINI_VERSION="v0.18.0"
ENV MINICONDA_VERSION="3-py39_4.9.2"
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y \
            apt-utils \
            bzip2 \
            build-essential \
	    ca-certificates \
	    curl \
	    git \
	    libfontconfig1 \
	    wget \
	    procps \
	    tzdata \
            unzip \
	    uuid-runtime 

WORKDIR /opt

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda${MINICONDA_VERSION}-Linux-x86_64.sh -O miniconda.sh && \
    chmod ugo+x miniconda.sh

RUN ./miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN wget --quiet https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /usr/bin/tini
RUN chmod +x /usr/bin/tini


SHELL ["/bin/bash", "-c","-l"]
RUN conda update conda
RUN conda config --add channels defaults &&\
    conda config --add channels bioconda &&\
    conda config --add channels conda-forge

SHELL ["/bin/sh", "-c"]
RUN mkdir -p /run/shm
RUN apt-get install -y bzip2 \
		       dh-autoreconf \
		       gcc \
		       libbz2-dev \
		       liblzma-dev \
		       libncurses5-dev \
		       libssl-dev \
		       make \
		       zlib1g-dev \
                       libbz2-dev \
                       libcurl4-openssl-dev

RUN conda install -y cython \
                     future \
                     matplotlib \
                     pandas \
                     scipy \
                     setuptools &&\
    conda install -y -c bioconda regex pysam

RUN cd /opt &&\
    git clone https://github.com/samtools/htslib.git &&\
    cd htslib &&\
    git checkout 1.10 &&\
    autoreconf -i  &&\ 
    ./configure --prefix=/usr/local   &&\ 
    make &&\
    make install


RUN git clone https://github.com/CGATOxford/UMI-tools.git &&\
    cd UMI-tools &&\
    git checkout 1.0.1 &&\
    python setup.py install


RUN cd /opt &&\
    git clone https://github.com/samtools/samtools.git &&\
    cd samtools &&\
    git checkout 1.10 &&\
    autoheader &&\
    autoconf -Wno-syntax &&\
    ./configure --prefix=/usr/local &&\
    make &&\
    make install


RUN cd /opt &&\
    git clone https://github.com/marcelm/cutadapt &&\
    cd cutadapt &&\
    git checkout v2.8 &&\
    python setup.py install

RUN cd /opt &&\
    git clone https://github.com/alexdobin/STAR &&\
    cd STAR &&\
    git checkout 2.6.0c &&\
    cd source &&\
    make && \
    cd ../bin/ &&\
    cp Linux_x86_64_static/* /usr/local/bin/ &&\
    cd /opt/ &&\
    rm -rf /opt/STAR

RUN cd /opt &&\
    git clone https://github.com/arq5x/bedtools2 &&\
    cd bedtools2 &&\
    git checkout v2.29.2 &&\
    make &&\
    cp bin/* /usr/local/bin/

ADD http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig /usr/local/bin/
ADD http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed /usr/local/bin/
ADD http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph /usr/local/bin/
ADD http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge /usr/local/bin/

RUN chmod +x /usr/local/bin/bedGraphToBigWig \
             /usr/local/bin/bigWigAverageOverBed \
	     /usr/local/bin/bigWigToBedGraph \
	     /usr/local/bin/bigWigMerge

RUN apt-get install -y pigz

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
