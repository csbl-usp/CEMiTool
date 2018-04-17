context('sample_annotation')

data(expr0)
data(sample_annot)

test_that('sample_annotation gives exactly the object received', {
	cem <- new_cem(expr0)
	sample_annotation(cem) <- sample_annot
	expect_identical(sample_annotation(cem), sample_annot)
})

test_that('sample_annotation returns errors when appropriate', {
	cem <- new_cem(expr0)
	expect_error(sample_annotation(cem, sample_name_column="foo") <- sample_annot)
	expect_error(sample_annotation(cem, class_column="foo") <- sample_annot)
})

test_that('sample_annotation returns a warning when there is only 1 sample in a group', {
	cem <- new_cem(expr0)
	one_row <- sample_annot[1,]
	bad_annot <- sample_annot[sample_annot$Class!="g0", ]
	bad_annot <- rbind(one_row, bad_annot)
	expect_warning(sample_annotation(cem) <- bad_annot)
})
