context('get_hubs')

cem0 <- new_cem()

test_that('get_hubs throws an error when there are no elements in adjacency matrix', {
	expect_error(get_hubs(cem0))
})
