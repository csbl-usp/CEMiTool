context('mod_names')

test_that('a warning is thrown if there are no modules', {
    cem <- new_cem()
    expect_warning(mod_names(cem))
})

test_that('modules are sorted by size', {
    data(cem)
    expect_true(identical(mod_names(cem), names(sort(table(cem@module$modules), decreasing=T))))
})
