context('filter')

data(cem)
data(expr0)
cem0 <- new_cem()

test_that('filter_genes throws an error when there is no expression file', {
	expect_error(filter_genes(data.frame()))
})

test_that('filter_genes throws an error when pct is zero', {
    expect_error(filter_genes(expr0, pct=0))
})

test_that('select_genes throws a warning when no genes are left after filtering', {
	expr_f <- filter_genes(expr0)
	expect_warning(selected <- select_genes(expr_f, filter_pval=-100))
})








