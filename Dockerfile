FROM mambaorg/micromamba:2.3.0-ubuntu22.04 AS app

ARG MYCOTOOLS_VER="1.0.0"
USER root

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="mambaorg/micromamba:2.3.0-ubuntu22.04"
LABEL dockerfile.version="1"
LABEL software="Mycotools"
LABEL software.version="${MYCOTOOLS_VER}"
LABEL description="Mycotools is a compilation of computational biology tools and database (MycotoolsDB/MTDB) software that facilitate large-scale comparative genomics."
LABEL website="https://github.com/xonq/mycotools"
LABEL license="https://github.com/xonq/mycotools/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="konkelzach@protonmail.com"

# this is unfortunately necessary to install ete4 
RUN micromamba install --name base -c conda-forge -c bioconda -c defaults legacy-cgi pip mycotools=${MYCOTOOLS_VER} && \
  eval "$(micromamba shell hook --shell bash)" && \
  micromamba activate base && \
  python3 -m pip install dna_features_viewer && \
  micromamba clean -a -f -y && \
  mkdir /data

ENV PATH="/opt/conda/bin/:${PATH}" \
    LC_ALL=C.UTF-8

# 'CMD' instructions set a default command when the container is run. This is typically 'tool --help.'
CMD [ "mtdb", "--help" ]

# 'WORKDIR' sets working directory
WORKDIR /data