context('sample_annotation')

data(expr)
data(sample_annot)

test_that('sample_annotation gives exactly the object received', {
		      cem <- new_cem(expr)
			  	  sample_annotation(cem) <- sample_annot
			      expect_identical(sample_annotation(cem), sample_annot)
})
