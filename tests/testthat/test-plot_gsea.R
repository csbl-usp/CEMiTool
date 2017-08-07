context('plot_gsea')

data(cem)
cem0 <- new_cem()

test_that('plot_gsea throws an error when there are no modules', {
    expect_error(cem0 <- plot_gsea(cem0))
})

test_that('plot_gsea throws an error when mod_gsea has not been run', {
    expect_error(cem <- plot_gsea(cem))
})
