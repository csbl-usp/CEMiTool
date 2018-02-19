#!/usr/bin/env Rscript

"CEMiTool - Co-Expression Modules identification Tool

Usage: cemitool.R EXPRSFILE  --output=<dir> [--sample-annot=<annot> --samples-column=<samplecol> --class-column=<classcol> --no-filter (--filter-pval=<p>|--ngenes=<ngenes>) --vst --eps --network-type=<nettype> --tom-type=<tomtype> --interactions=<inter> --pathways=<gmt> --ora-pvalue=<p> --cor-method=<cor> --no-merge --rank-method --min-module-size=<min> --diss-thresh=<thresh> --center-func=<fun> --directed --verbose]

Input:
  EXPRSFILE                         a normalized expression file .tsv format

Options:
  -h --help                          show this help message
  --version                          show program version
  -s <annot> --sample-annot=<annot>  sample annotation, must have a column with sample names and class
  --samples-column=<samplecol>       the column name containing sample names in template file [default: SampleName]
  --class-column=<classcol>          the column name containing classes in template file [default: Class]
  -i <int> --interactions=<int>      gene interaction file, must have two columns
  -p <gmt> --pathways=<gmt>          GMT file name (Gene Matrix Transposed format)
  --ora-pvalue=<p>                   p-value cutoff to be used on over representation analysis [default: 0.05]
  --no-filter                        does not filter the expression data.frame
  --filter-pval=<p>                  p-value to be used in the filtering step [default: 0.1]
  --vst                              apply Variance Stabilizing Transformation
  --ngenes=<ngenes>                  number of genes remaining after filtering
  --eps=<eps>                        epsilon [default: 0.1]
  -c <cor> --cor-method=<cor>        correlation method (spearman or pearson) [default: pearson]
  --network-type=<nettype>           network type, 'signed' or 'unsigned' [default: unsigned]
  --tom-type=<nettype>               TOM type, 'signed' or 'unsigned' [default: signed]
  --no-merge                         does not merge related modules based on eigengene similarity
  --rank-method                      rank method [default: mean]
  --min-module-size=<min>            minimum module size [default: 30]
  --diss-thresh=<thresh>             module merging correlation threshold for eigengene similarity [default: 0.8]
  --center-func=<fun>                metric used for centering [default: mean]
  --directed                         the interactions are directed
  -o <dir> --output=<dir>            output directory
  -v --verbose                       display progress messages

Authors:
  Pedro S T Russo - pedro.russo at usp.br
  Gustavo R Ferreira - gustavo.rodrigues.ferreira at usp.br
  Lucas E Cardozo - lucasecardozo at usp.br
  Matheus C Burger - burger at usp.br
  Thiago D C Hirata - thiagodch at gmail.com
  Diogenes S Lima - diogenes.lima at usp.br
  Fernando M Passos - fmarcon at usp.br
  Raul A Carrasco 
  Melissa Lever - melissalever at gmail.com 
  Vinicius Maracaja-Coutinho
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

    # expression file
    if(p$verbose){
        message("Reading expression file ...")
    }
    p$expr <- data.table::fread(parameters[["exprsfile"]], data.table=FALSE)
    
    # remove the column containing gene symbols
    if(p$verbose) {
        message("Setting row names ...")
    }
    rownames(p$expr) <- p$expr[,1]
    p$expr[,1] <- NULL

    # verify if all columns are numeric
    if(!all(sapply(p$expr, is.numeric))){
        stop("Please make sure that your expression file have only numeric values.")
    }

    # sample annotation file
    if("sample_annot" %in% names(parameters)){
        if(p$verbose){
            message("Reading sample annotation file...")
        }
        p$annot <- data.table::fread(parameters[["sample_annot"]], data.table=FALSE)

        # sample name column in sample annotation file
        p$sample_name_column <- parameters[["samples_column"]]

        # class column in sample annotation file
        p$class_column <- parameters[["class_column"]] 
    }


    # Should filter ?
    p$filter <- !parameters[["no_filter"]]

    # filter p-value or number of genes after filtering
    if("ngenes" %in% names(parameters)) {
        p$n_genes <- as.numeric(parameters[["ngenes"]])
    } else {
        p$filter_pval <- as.numeric(parameters[["filter_pval"]])
    }

    p$apply_vst <- parameters[["vst"]]

    # gmt list
    if("pathways" %in% names(parameters)){
        if(p$verbose){
            message("Reading GMT file ...")
        }
        p$gmt <- read_gmt(parameters[["pathways"]])
    }

    # interactions file
    if("interactions" %in% names(parameters)){
        if(p$verbose){
            message("Reading interactions file ...")
        }
        p$interactions <- data.table::fread(parameters[["interactions"]], data.table=FALSE)
    }

    # correlation method
    p$cor_method <- parameters[["cor_method"]]

    # epsilon
    p$eps <- parameters[["eps"]]
      
    # network type
    p$network_type <- parameters[["network_type"]]
      
    # TOM type
    p$tom_type <- parameters[["tom_type"]]
      
    # should similar modules be merged ?
    p$merge_similar <- !parameters[["no_merge"]]

    # p-value cutoff to be used on over representation analysis
    p$ora_pval <- as.numeric(parameters[["ora_pvalue"]])

    # minimum size of a module
    p$min_ngen <- as.numeric(parameters[["min_module_size"]])

    # dissimilarity threshold
    p$diss_thresh <- as.numeric(parameters[["diss_thresh"]])

    p$rank_method <- parameters[["rank_method"]]

    p$center_func <- parameters[["center_func"]]
    # should create figures ?
    p$plot <- TRUE

    # Is the interactions file directed ?
    p$directed <- parameters[["directed"]]

    # CEMiTool
    if(p$verbose){
        message("Running CEMiTool ...")
    }

    cem <- do.call(cemitool, p)

    # Save files
    if(p$verbose){
        message("Writing CEMiTool results...")
    }
    write_files(cem, directory=parameters[["output"]], force=TRUE)

    # Generate reports
    generate_report(cem, directory=parameters[["output"]])

}
