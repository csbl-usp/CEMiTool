context('plot_ora')

data(cem)
cem0 <- new_cem()

test_that('plot_ora throws an error when there are no modules', {
	expect_error(cem0 <- plot_ora(cem0))
})

test_that('plot_ora throws an error when mod_ora has not been run', {
	expect_error(cem <- plot_ora(cem))
})
