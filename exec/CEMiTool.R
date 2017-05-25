#!/usr/bin/env Rscript

"CEMiTool - Co-Expression Modules Identification Tool - Co-Expression Module analysis made easy

Usage: cemitool.R EXPRSFILE  --output=<DIR> [--gene-column=<GENECOL> --sample-annot=<ANNOT> --samples-column=<SAMPLECOL> --class-column=<CLASSCOL> --dontfilter (--filter-pval=<p>|--genenum=<GENNUM>) --network-type=<NETTYPE> --tom-type=<NETTYPE> --gmtfile=<GMT> --interact=<INT> --correlation=<COR> --dontmerge --ora-pvalue=<p> --min-module=<MIN> --diss-thresh=<THRESH> --directed --verbose]

Input:
  EXPRSFILE                         a normalized expression file .tsv format

Options:
  -h --help                         show this help message
  --version                         show program version
  --gene-column=<GENECOL>           the column name containing gene symbols in expression file [default: Symbol]
  -s <ANNOT> --sample-annot=<ANNOT> sample annotation, must have a column with sample names and class
  --samples-column=<SAMPLECOL>      the column name containing sample names in template file [default: SampleName]
  --class-column=<CLASSCOL>         the column name containing classes in template file [default: Class]
  --dontfilter                      filter the expression data.frame
  --filter-pval=<p>                 p-value to be used in the filtering step [default: 0.1]
  --genenum=<GENNUM>                number of genes remaining after filtering
  --network-type=<NETTYPE>          network type, 'signed' or 'unsigned' [default: 'unsigned']
  --tom-type=<NETTYPE>              TOM type, 'signed' or 'unsigned' [default: 'signed']
  --gmtfile=<GMT>                   GMT file name (Gene Matrix Transposed format)
  -i <INT> --interact=<INT>         gene interaction file, must have two columns
  -c <COR> --correlation=<COR>      correlation method (spearman or pearson) [default: pearson]
  --dontmerge                       merge related modules based on eigengene similarity
  --ora-pvalue=<p>                  p-value cutoff to be used on over representation analysis [default: 0.05]
  --min-module=<MIN>                minimum module size [default: 30]
  --diss-thresh=<THRESH>             module merging correlation threshold for eigengene similarity [default: 0.8]
  --directed                        the interactions are directed
  -o <DIR> --output=<DIR>           output directory
  -v --verbose                      display progress messages [default: TRUE]

Authors:
  Gustavo R Ferreira - gustavo.rodrigues.ferreira at usp.br
  Pedro S T Russo - pedro.russo at usp.br
  Matheus C Burger - burger at usp.br
  Thiago D C Hirata - thiagodch at gmail.com
  Lucas E Cardozo - lucasecardozo at usp.br
  Helder I Nakaya - hnakaya at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

if (!interactive()) {
    # Get and check arguments.
    suppressMessages(library("docopt"))
    arg <- docopt(doc, version="0.0.1\n", strict=TRUE)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))
    parameters <- arg

    print(parameters)

    ## RUN
    library("CEMiTool")

    # parameters
    p <- list()

    # verbosity
    p$verbose <- parameters[["verbose"]]

    # gene column
    gene_column <- parameters[["gene_column"]]

    # expression file
    if(p$verbose){
        message("Reading expression file ...")
    }
    p$expr <- read.delim(parameters[["exprsfile"]], row.names=NULL)
    if(!gene_column %in% colnames(p$expr)) {
        stop("Please give a valid column containing gene symbols.")
    }

    # remove the column containing gene symbols
    if(p$verbose) {
        message("Setting row names ...")
    }
    rownames(p$expr) <- p$expr[[gene_column]]
    col_to_remove <- which(gene_column %in% colnames(p$expr))
    p$expr <- p$expr[, -col_to_remove]

    # verify if all columns are numeric
    if(!all(sapply(p$expr, is.numeric))){
        stop("Please make sure that your expression file have only numeric values.")
    }

    # sample annotation file
    if("sample_annot" %in% names(parameters)){
        if(p$verbose){
            message("Reading sample annotation file...")
        }
        p$annot <- read.delim(parameters[["sample_annot"]], row.names=NULL)

        # sample name column in sample annotation file
        p$sample_name_column <- parameters[["samples_column"]]
        if(!p$sample_name_column %in% colnames(p$annot)) {
            stop("Please give a valid column containing sample names in sample annotation file.")
        }

        # class column in sample annotation file
        p$class_column <- parameters[["class_column"]] 
        if(!p$class_column %in% colnames(p$annot)){
            stop("Please give a valid column containing classes in sample annotation file.")
        }
    }


    # Should filter ?
    p$filter <- !parameters[["dontfilter"]]

    # filter p-value or number of genes after filtering
    if("genenum" %in% names(parameters)) {
        p$n_genes <- as.numeric(parameters[["genenum"]])
    } else {
        p$filter_pval <- as.numeric(parameters[["filter_pval"]])
    }


    # gmt list
    if("gmtfile" %in% names(parameters)){
        if(p$verbose){
            message("Reading GMT file ...")
        }
        p$gmt <- read_gmt(parameters[["gmtfile"]])
    }

    # interactions file
    if("interact" %in% names(parameters)){
        if(p$verbose){
            message("Reading interactions file ...")
        }
        p$interactions <- read.delim(parameters[["interact"]], row.names=NULL)
    }

    # correlation method
    p$cor_method <- parameters[["correlation"]] 
      
    # network type
    p$network_type <- parameters[["network_type"]]
      
    # TOM type
    p$tom_type <- parameters[["tom_type"]]
      
    # should similar modules be merged ?
    p$merge_similar <- !parameters[["dontmerge"]]

    # p-value cutoff to be used on over representation analysis
    p$ora_pval <- as.numeric(parameters[["ora_pvalue"]])

    # minimum size of a module
    p$min_ngen <- as.numeric(parameters[["min_module"]])

    # dissimilarity threshold
    p$diss_thresh <- as.numeric(parameters[["diss_thresh"]])

    # should create figures ?
    p$plot <- TRUE

    # Is the interactions file directed ?
    p$directed <- parameters[["directed"]]

    # CEMiTool
    if(p$verbose){
        message("Running CEMiTool ...")
    }
    save(file="CEMiTool.RData", list=ls())
    cem <- do.call(cemitool, p)

    # Save files
    if(p$verbose){
        message("Writing CEMiTool results...")
    }
    write_files(cem, directory=parameters[["output"]], force=TRUE)

    # Generate reports
    generate_report(cem, directory=parameters[["output"]])

}
