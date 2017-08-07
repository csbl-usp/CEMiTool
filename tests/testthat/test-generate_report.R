context('generate_report')

cem0 <- new_cem()

test_that('generate_report throws an error when there are no modules', {
	expect_error(generate_report(cem0))
})
