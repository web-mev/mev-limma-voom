From debian:stretch

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libopenblas-dev

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/software/miniconda.sh && \
  bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# Install limma and dependencies:
ADD spec_file.txt /opt/software/
RUN conda create --name limma_env --file /opt/software/spec_file.txt
ADD limma.R /opt/software/

# add to the PATH. Otherwise, conda can swallow return codes
ENV PATH="/opt/software/miniconda/envs/limma_env/bin:${PATH}"

ENTRYPOINT ["/bin/bash"]

