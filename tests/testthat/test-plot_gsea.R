context('plot_gsea')

data(cem)
cem0 <- new_cem()

test_that('plot_gsea throws an error when there are no modules', {
    expect_error(cem0 <- plot_gsea(cem0))
})

test_that('plot_gsea throws an error when mod_gsea has not been run', {
    expect_error(cem <- plot_gsea(cem))
})

test_that('plot_gsea throws a warning when dfs in cem@enrichment are empty', {
    bad_cem <- mod_gsea(cem)
    bad_cem@enrichment <- lapply(bad_cem@enrichment, function(x) x <- x[c(), ])
    expect_warning(bad_cem <- plot_gsea(bad_cem))
})
