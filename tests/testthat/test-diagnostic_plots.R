context("Diagnostic plots")

data(expr0)
cem0 <- new_cem()

test_that('diagnostic plot functions throw errors when there is no expression file', {
	expect_error(plot_sample_tree(cem0))
	expect_error(plot_hist(cem0))
	expect_error(plot_qq(cem0))
	expect_error(plot_mean_var(cem0))
})
