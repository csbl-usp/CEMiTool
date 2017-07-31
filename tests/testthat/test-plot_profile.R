context('plot_profile')

data(cem)
cem0 <- new_cem()

test_that('plot_profile throws an error when there are no modules', {
	expect_error(plot_profile(cem0))
})

test_that('plot_profile returns as many plots as modules', {
	cem <- plot_profile(cem)
    expect_equal(length(module_names(cem)), length(cem@profile_plot))
})
