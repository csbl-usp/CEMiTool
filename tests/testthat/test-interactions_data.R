context('interactions_data')

data(cem)
cem0 <- new_cem()
int_df <- read.delim(system.file("extdata", "interactions.tsv", 
								          package = "CEMiTool"))

test_that('interactions_data throws an error when there are no modules', {
	expect_error(interactions_data(cem0) <- int_df)
})

