FROM r-base:3.5.1

RUN apt-get update &&\
    apt-get install -y --no-install-recommends --allow-downgrades libssl-dev\
    libssh2-1-dev pandoc libcurl4=7.61.0-1  libcurl4-openssl-dev libxml2-dev\
    libudunits2-dev

ADD . /CEMiTool

WORKDIR /CEMiTool

RUN /usr/bin/Rscript docker/install-deps.R
