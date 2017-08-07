context('mod_ora')

data(cem)

gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

test_that('mod_ora throws warning when there are no modules in CEMiTool object', {
	cem0 <- new_cem()
	expect_warning(mod_ora(cem0, gmt_in))
})

gmt_bad <- gmt_in[1:10, ]

test_that('mod_ora throws warning when enrichment is NULL', {
	expect_warning(cem <- mod_ora(cem, gmt_bad))
})
