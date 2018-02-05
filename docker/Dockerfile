FROM csblusp/cemitool:builder

RUN R CMD install .

ENTRYPOINT ["/usr/bin/Rscript", "/usr/local/lib/R/site-library/CEMiTool/exec/CEMiTool.R"]
