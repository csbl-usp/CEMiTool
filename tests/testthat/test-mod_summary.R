context('mod_summary')

cem0 <- new_cem()

test_that('mod_summary returns an error when there are no modules', {
	expect_error(mod_summary(cem0))
})
