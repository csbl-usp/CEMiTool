context('get_hubs')

data(cem)
cem0 <- new_cem()

test_that('get_hubs throws an error when there are no elements in adjacency matrix', {
	expect_error(get_hubs(cem0))
})

test_that('get_hubs returns a named list', {
    expect_is(get_hubs(cem), "list")
    expect_equal(sort(names(get_hubs(cem))), sort(mod_names(cem)))
})

test_that('get_hubs returns all genes in a module when parameter "all == TRUE"', {
    expect_equal(sort(names(get_hubs(cem)[["M1"]])), sort(module_genes(cem, module="M1")))
})
