context('find_modules')

cem0 <- new_cem()

test_that('find_modules throws an error when there is no expression file', {
	expect_error(find_modules(cem0))
})
