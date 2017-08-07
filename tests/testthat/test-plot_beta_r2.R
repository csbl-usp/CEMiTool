context('plot_beta_r2')

cem0 <- new_cem()

test_that('plot_beta_r2 throws an error if there is no fit_indices slot', {
	expect_error(plot_beta_r2(cem0))
})


