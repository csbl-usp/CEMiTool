#' @import rmarkdown
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
#' @param ... parameters to rmarkdown::render
#'
#' @examples
#' # Get example expression data
#' data(expr)
#' # Run CEMiTool with analysis plots
#' cem <- cemitool(expr, plot=TRUE)
#' \dontrun{
#' generate_report(cem, output <- format=c("pdf_document", "html_document"))
#' }
#' @rdname generate_report
#' @export
setGeneric('generate_report', function(cem, ...) {
               standardGeneric('generate_report')
})

#' @rdname generate_report
setMethod('generate_report', signature('CEMiTool'),
          function(cem, max_rows_ora=50, title="Report", directory="./reports", ...) {
              rmd <- system.file("report", "report.Rmd", package = "CEMiTool")
              rmarkdown::render(rmd, output_dir=directory, intermediates_dir=directory, ...)
          })
