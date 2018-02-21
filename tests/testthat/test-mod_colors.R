context('mod_colors')

test_that('mod_colors input is a named character vector', {
	mod_cols <- mod_colors(cem)
	expect_is(mod_cols, "character")
	expect_named(mod_cols, mod_names(cem))
})

test_that('mod_colors returns the correct errors when given bad input', {
    mock_colors <- c("a", "b", "c", "d", "e")
    expect_error(mod_colors(cem) <- mock_colors)
    names(mock_colors) <- c("A", "B", "C", "D", "E")
    expect_error(mod_colors(cem) <- mock_colors)
})
