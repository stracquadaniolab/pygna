FROM python:3.8-slim

LABEL org.opencontainers.image.title="pygna"
LABEL org.opencontainers.image.description="PyGNA: a Python framework for Geneset Network Analysis."
LABEL org.opencontainers.image.url="https://github.com/stracquadaniolab/pygna"
LABEL org.opencontainers.image.documentation="https://github.com/stracquadaniolab/pygna"
LABEL org.opencontainers.image.source="https://github.com/stracquadaniolab/pygna"
LABEL org.opencontainers.image.revision="3.4.0"
LABEL org.opencontainers.image.vendor="stracquadaniolab"
LABEL org.opencontainers.image.authors="Viola Fanfani, Fabio Cassano, Giovanni Stracquadanio"

# update sources and install tini
RUN apt-get update --fix-missing && apt-get autoremove \
    && apt-get clean \
    && apt-get install tini \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# copy the package and install pygna
WORKDIR /opt
COPY . .
RUN python setup.py install

ENTRYPOINT ["/usr/bin/tini", "--", "pygna"]

