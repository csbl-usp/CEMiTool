context('expr_data')

data(expr0)

test_that('expr_data gives exactly the object received', {
	cem <- new_cem()
	expr_data(cem) <- expr0
	expect_identical(expr_data(cem), expr0)
})
