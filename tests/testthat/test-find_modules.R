context('find_modules')

data(expr0)
cem0 <- new_cem()

test_that('find_modules functions throw errors when there is no expression file', {
	expect_error(find_modules(cem0))
	expect_error(get_beta_data(cem0))
	expect_error(get_mods(cem0))
	expect_error(get_merged_mods(cem0))
})

cem <- cem0
expr_data(cem) <- expr0
#cem <- filter_expr(cem)

test_that('find_modules functions throw errors when there is no fit_indices slot', {
	expect_error(get_cemitool_r2_beta(cem))
	expect_error(get_phi(cem))
	expect_error(get_connectivity(cem))

})

test_that('find_modules throws errors with innapropriate set_beta values', {
	expect_error(invisible(capture.output(find_modules(cem, set_beta=NA))))
	expect_error(invisible(capture.output(find_modules(cem, set_beta=13))))
	expect_error(invisible(capture.output(find_modules(cem, set_beta="foo"))))
})

test_that('find_modules throws error when set_beta and force_beta are given together', {
	expect_error(invisible(capture.output(find_modules(cem, set_beta=7, force_beta=TRUE))))
})

