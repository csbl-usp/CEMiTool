context('filter_expr')

data(cem)
cem0 <- new_cem()

test_that('expression file is getting the correct number of genes post-filtering', {
	cem <- filter_expr(cem)
	expect_equal(nrow(expr_data(cem)), length(cem@selected_genes))
})

test_that('filter_expr throws an error when both pval and n_genes are given', {
	expect_error(cem <- filter_expr(cem, pval=0.1, n_genes=200))
})

test_that('filter_expr throws an error when there is no expression file', {
	expect_error(filter_expr(cem0))
})

test_that('filter_expr throws a warning when no genes are left after filtering', {
	expect_warning(cem <- filter_expr(cem, pval=0))
})
