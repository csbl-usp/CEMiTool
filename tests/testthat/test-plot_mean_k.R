context('plot_mean_k')

cem0 <- new_cem()

test_that('plot_mean_k throws an error if there is no fit_indices slot', {
	expect_error(plot_mean_k(cem0))
}) 
