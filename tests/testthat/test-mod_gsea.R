context('mod_gsea')

cem0 <- new_cem()
data(expr0)
data(sample_annot)
data(cem)

test_that('mod_gsea throws warning when there is no expression data', {
	expect_warning(mod_gsea(cem0))
})

expr_data(cem0) <- expr0

test_that('mod_gsea throws warning when there is no sample_annotation data', {
	expect_warning(mod_gsea(cem0))
})

sample_annotation(cem0) <- sample_annot

test_that('mod_gsea throws warning when there are no modules', {
	expect_warning(mod_gsea(cem0))
})

test_that('mod_gsea throws error when expression has more samples than annotation', {
	class1 <- unique(sample_annot$Class)[1]
	sample_annot_bad <- sample_annot[sample_annot$Class!=class1,]
	sample_annotation(cem) <- sample_annot_bad
	expect_error(mod_gsea(cem))
})

test_that('mod_gsea throws warning when expression has less samples than annotation', {
	expr_bad <- expr0[, 1:round(ncol(expr0)/2)]
	expr_data(cem) <- expr_bad
	expect_warning(mod_gsea(cem))	
})




