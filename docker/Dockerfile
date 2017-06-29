FROM r-base:3.4.0

RUN apt-get update &&\
    apt-get install -y --no-install-recommends libssl-dev\
    libssh2-1-dev libcurl4-openssl-dev pandoc

RUN R --slave -e "install.packages('devtools');\
                  source('https://bioconductor.org/biocLite.R');\
                  biocLite(c('WGCNA', 'gRbase'));\
                  devtools::install_github('csbl-usp/CEMiTool')"

ENTRYPOINT ["/usr/bin/Rscript", "/usr/local/lib/R/site-library/CEMiTool/exec/CEMiTool.R"]
