## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      cache=TRUE)

## ------------------------------------------------------------------------
library("CEMiTool")
# load expression data
data(expr)
head(expr[,1:4])

## ---- results='hide'-----------------------------------------------------
cem <- cemitool(expr)

## ------------------------------------------------------------------------
print(cem)
# or, more conveniently, just call `cem`
cem

## ------------------------------------------------------------------------
# inspect modules
nmodules(cem)
head(module_genes(cem))

## ---- eval=FALSE---------------------------------------------------------
#  generate_report(cem, output_format=c("pdf_document", "html_document"))

## ------------------------------------------------------------------------
write_files(cem, directory="./reports", force=TRUE)

## ------------------------------------------------------------------------
# load your sample annotation data
data(sample_annot)
head(sample_annot)

## ---- results='hide'-----------------------------------------------------
# run cemitool with sample annotation
cem <- cemitool(expr,sample_annot)

## ------------------------------------------------------------------------
# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")

## ------------------------------------------------------------------------
# perform GSEA of modules across your experimental classes
cem <- mod_gsea(cem)

## ------------------------------------------------------------------------
# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]

## ------------------------------------------------------------------------
# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

## ------------------------------------------------------------------------
# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

## ------------------------------------------------------------------------
# plot ora results
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[1]

## ------------------------------------------------------------------------
# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

## ------------------------------------------------------------------------
# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module

## ---- eval=FALSE---------------------------------------------------------
#  # run cemitool
#  library(ggplot2)
#  cem <- cemitool(expr, sample_annot, gmt_in, interactions=int_df,
#                  filter=TRUE, plot=TRUE)
#  # create report
#  generate_report(cem, output_format=c("pdf_document", "html_document"))
#  
#  # write analysis results into files
#  write_files(cem, directory=".", force=TRUE)

