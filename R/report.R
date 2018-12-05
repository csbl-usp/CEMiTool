#' @importFrom rmarkdown render
#' @import knitr
#' @importFrom DT datatable
#' @import htmltools
#'
NULL

#' CEMiTool report
#'
#' Creates report for CEMiTool results
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param directory Directory name for results.
#' @param max_rows_ora maximum number of rows in Over Representation Analysis table results
#' @param title Character string with the title of the report.
#' @param force If the directory exists, execution will not stop. 
#' @param ... parameters to rmarkdown::render
#' 
#' @return An HTML file with an interactive report of CEMiTool analyses.
#'
#' @examples
#' \dontrun{  
#' # Get example CEMiTool object
#' data(cem)
#' generate_report(cem, output_format=c("pdf_document", "html_document"))
#' }
#'
#' @rdname generate_report
#' @export
setGeneric('generate_report', function(cem, ...) {
    standardGeneric('generate_report')
})

#' @rdname generate_report
setMethod('generate_report', signature('CEMiTool'),
    function(cem, max_rows_ora=50, title="Report", directory="./Reports/Report", force=FALSE, ...) {
        if(is.null(unique(cem@module$modules))){
            stop("No modules in CEMiTool object! Did you run find_modules()?")
        }
        if(dir.exists(directory)){
            if(!force){
                stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
              }
        }else{
            dir.create(directory, recursive=TRUE)
        }
        rmd <- system.file("report", "report.Rmd", package = "CEMiTool")
        rmarkdown::render(rmd, output_dir=directory, intermediates_dir=directory, quiet=TRUE, ...)
    })
