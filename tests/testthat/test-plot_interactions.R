context('plot_interactions')

data(cem)
cem0 <- new_cem()
int_df <- read.delim(system.file("extdata", "interactions.tsv", 
								          package = "CEMiTool"))

test_that('plot_interactions throws an error when there are no modules', {
	expect_error(plot_interactions(cem0))
})

test_that('plot_interactions throws an error when interactions_data has not been run', {
	expect_error(plot_interactions(cem))
})



