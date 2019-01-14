context('new_cem')

data(expr0)
data(sample_annot)

cem0 <- new_cem()
cem1 <- new_cem(expr0)
cem2 <- new_cem(expr0, sample_annot)

test_that('objects constructed are of class CEMiTool', {
	expect_is(cem0, 'CEMiTool')
	expect_is(cem1, 'CEMiTool')
	expect_is(cem2, 'CEMiTool')
})

test_that('objects constructed receive the expected data', {
	expect_identical(expr0, expr_data(cem2, filter=FALSE, apply_vst=FALSE))
	expect_identical(sample_annot, sample_annotation(cem2))
})

test_that('new_cem returns errors when appropriate', {
	expect_error(cem <- new_cem(expr0, sample_annot, sample_name_column="foo"))
	expect_error(cem <- new_cem(expr0, sample_annot, class_column="bar"))
    expect_error(cem <- new_cem(expr0, sample_annot, filter="baz"))
    expect_error(cem <- new_cem(expr0, sample_annot, apply_vst="bum"))
})


