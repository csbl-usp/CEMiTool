context('mod_names')

test_that('a warning is thrown if there are no modules', {
	cem <- new_cem()
	expect_warning(mod_names(cem))
})    
