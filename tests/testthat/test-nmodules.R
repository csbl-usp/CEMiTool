context('nmodules')

test_that('a warning is thrown if there are no modules', {
	cem <- new_cem()
	expect_warning(nmodules(cem))
})
