context('adj_data')

data(expr0)
cem0 <- new_cem()
cem <- new_cem(expr0)
cem <- get_adj(cem, beta=7)
adj <- adj_data(cem)

test_that('adj_data throws an error when CEMiTool object has no expression data', {
	expect_error(adj_data(cem0) <- adj)
})

test_that('adj_data throws an error when object provided is not of class matrix', {
	expect_error(adj_data(cem0) <- "foo")
	expect_error(adj_data(cem0) <- 1)	
	expect_error(adj_data(cem0) <- NA)
	expect_error(adj_data(cem0) <- NULL)
	expect_error(adj_data(cem0) <- NaN)
})

test_that('adj_data throws an error when rownames in the object provided do not reflect those in expression', {
	bad_adj <- adj
	rownames(bad_adj) <- seq_len(nrow(bad_adj))
	expect_error(adj_data(cem) <- bad_adj)
})

