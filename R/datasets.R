#' Yellow Fever gene expression data from GEO study GSE13485
#'
#' Modified data from a yellow fever vaccination study by Querec et al, 2009.
#' In order to reduce package size, only the 4000 genes with the highest 
#' variance were selected for this dataset. 
#'
#' @name expr 
#' @docType data
#' @usage data(expr)
#' @format An object of class \code{data.frame}
#' @keywords datasets
#' @references Querec TD, Akondy RS, Lee EK, Cao W et al. Systems biology 
#' approach predicts immunogenicity of the yellow fever vaccine in humans. 
#' Nat Immunol 2009 Jan;10(1):116-25. PMID: 19029902
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/19029902}{PubMed}
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse13485}{GEO}
#' @examples
#' data(expr)
#' # Run CEMiTool analysis
#' \dontrun{cemitool(expr)}
"expr"

#' Yellow Fever Sample Annotation data 
#'
#' Modified data from a yellow fever vaccination study by Querec et al, 
#' 2009. This dataset, together with \code{expr} can be used as input for 
#' CEMiTool functions
#'
#' @name sample_annot
#' @docType data
#' @usage data(sample_annot)
#' @format An object of class \code{data.frame} 
#' @keywords datasets
#' @references Querec TD, Akondy RS, Lee EK, Cao W et al. Systems biology 
#' approach predicts immunogenicity of the yellow fever vaccine in humans. 
#' Nat Immunol 2009 Jan;10(1):116-25. PMID: 19029902
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/19029902}{PubMed}
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse13485}{GEO}
#' @examples
#' data(expr)
#' data(sample_annot)
#' # Run CEMiTool analysis
#' \dontrun{cemitool(expr, sample_annot)}
"sample_annot"

#' CEMiTool Object 
#'
#' This object can be used as input for CEMiTool functions. Data used are from 
#' \code{expr} and \code{sample_annot}.
#'
#' @name cem
#' @docType data
#' @usage data(cem)
#' @format An object of class \code{CEMiTool}
#' @keywords data
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Read example gmt file
#' gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
#' gmt_in <- clusterProfiler::read.gmt(gmt_fname)
#' # Read example interactions file
#' int_df <- read.delim(system.file("extdata", "interactions.tsv", 
#'         package = "CEMiTool"))
#' # Insert interactions data
#' interactions_data(cem) <- int_df
#' # Run analyses
#' cem <- mod_gsea(cem)
#' cem <- plot_gsea(cem)
#' cem <- mod_ora(cem, gmt_in)
#' cem <- plot_ora(cem)
#' \dontrun{generate_report(cem)}
"cem"



















