context('expr_data')

data(expr)

test_that('expr_data gives exactly the object received', {
	cem <- new_cem()
	expr_data(cem) <- expr
	expect_identical(expr_data(cem), expr)
})
