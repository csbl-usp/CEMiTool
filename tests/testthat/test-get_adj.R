context("get_adj")

data(expr0)
cem0 <- new_cem()
cem <- new_cem(expr0)

test_that("get_adj returns an error when no beta value is provided.", {
	expect_error(cem <- get_adj(cem))
})

test_that("get_adj returns an error when CEMiTool object has no expression data", {
	expect_error(cem0 <- get_adj(cem0, beta=7))	
})

test_that("get_adj returns an object of class matrix", {
	cem <- get_adj(cem, beta=7)
	expect_is(cem@adjacency, "matrix")
})
