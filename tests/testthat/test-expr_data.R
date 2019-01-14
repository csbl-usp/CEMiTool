context('expr_data')

data(expr0)

test_that('expr_data gives exactly the object received', {
	cem <- new_cem()
	expr_data(cem) <- expr0
	expect_identical(expr_data(cem, filter=FALSE, apply_vst=FALSE), expr0)
})

test_that('expr_data is getting filter and apply_vst arguments from the cem object', {
    # filter = TRUE, apply_vst = TRUE
    cem <- new_cem(expr0, filter=TRUE, apply_vst=TRUE)
    expr <- expr_data(cem)
    expr_T_T <- filter_genes(expr0, apply_vst=TRUE)
    selected <- select_genes(expr_T_T)
    expect_identical(expr, expr_T_T[selected, ])

    # filter = TRUE, apply_vst = FALSE
    cem <- new_cem(expr0, filter=TRUE, apply_vst=TRUE)
    expr <- expr_data(cem)
    expr_T_F <- filter_genes(expr0, apply_vst=FALSE)
    selected <- select_genes(expr_T_F)
    expect_identical(expr, expr_T_F[selected, ])

    # filter = FALSE, apply_vst = FALSE
    cem <- new_cem(expr0, filter=FALSE, apply_vst=FALSE)
    expr <- expr_data(cem)
    expect_identical(expr, expr0)

    # filter = FALSE, apply_vst = TRUE
    cem <- new_cem(expr0, filter=FALSE, apply_vst=TRUE)
    expect_warning(expr <- expr_data(cem))
    expect_identical(expr, expr0)
})

test_that('expr_data prioritizes given arguments over those in cem object', {
    cem <- new_cem(expr0, filter=TRUE, apply_vst=TRUE)
    expr <- expr_data(cem, FALSE, FALSE)
    expect_identical(expr, expr0)
    expect_true(cem@parameters$filter)
    expect_true(cem@parameters$apply_vst)
})

test_that('expr_data throws errors for incorrect argument types', {
    expect_error(expr_data(cem, filter="foo"))
    expect_error(expr_data(cem, apply_vst=42))
})

