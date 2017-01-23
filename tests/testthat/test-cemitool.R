context('Cemitool methods')

data(exprs)
data(sample_annotation)

cem_base <- new('CEMiTool',
           expression=exprs,
           sample_annotation=sample_annotation,
           sample_name_column='Sample')

ppi <- system.file('extdata', 'interactions.tsv', package='CEMiTool')
ppi <- fread(ppi, data.table=F)

gmt <- system.file('extdata', 'pathways.gmt', package='CEMiTool')
gmt <- read_gmt(gmt)

test_that('all methods of signature CEMiTool returns CEMiTool objects', {
    cem <- filter_expr(cem_base)
    expect_equals(class(cem), 'CEMiTool') 
    
    cem <- find_modules(cem)
    expect_equals(class(cem), 'CEMiTool') 
    
    cem <- split_modules(cem)
    expect_equals(class(cem), 'CEMiTool')

    cem <- include_interactions(cem, ppi)
    expect_equals(class(cem), 'CEMiTool')

    cem <- mod_gsea(cem)
    expect_equals(class(cem), 'CEMiTool')

    cem <- mod_ora(cem) 
    expect_equals(class(cem), 'CEMiTool')

    cem <- plot_profile(cem)
    expect_equals(class(cem), 'CEMiTool')

    cem <- plot_gsea(cem)
    expect_equals(class(cem), 'CEMiTool')
    
    cem <- plot_ora(cem)
    expect_equals(class(cem), 'CEMiTool')
})
