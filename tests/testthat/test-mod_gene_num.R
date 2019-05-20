context('mod_gene_num')

test_that('an error is thrown if there are no modules', {
    cem <- new_cem()
    expect_error(mod_gene_num(cem))
})

test_that('an error is thrown if given a module not present in the object', {
    data(cem)
    expect_error(mod_gene_num(cem, "M99"))
})
