context('Cemitool methods')

data(expr)
data(sample_annot)

cem_base <- new('CEMiTool',
           expression=expr,
           sample_name_column='Sample')

sample_annotation(cem_base) <- sample_annot

ppi <- system.file('extdata', 'interactions.tsv', package='CEMiTool')
ppi <- fread(ppi, data.table=FALSE)

gmt <- system.file('extdata', 'pathways.gmt', package='CEMiTool')
gmt <- read_gmt(gmt)

cem_filt <- filter_expr(cem_base, 1)
cem_fm <- find_modules(cem_filt)


test_that('all methods of signature CEMiTool returns CEMiTool objects', {
    expect_equal(class(cem_filt)[1], 'CEMiTool') 
    
    expect_equal(class(cem_fm)[1], 'CEMiTool') 
    
    cem <- split_modules(cem_fm)
    expect_equal(class(cem)[1], 'CEMiTool') 

    cem <- include_interactions(cem, ppi)
    expect_equal(class(cem)[1], 'CEMiTool') 

    cem <- mod_gsea(cem)
    expect_equal(class(cem)[1], 'CEMiTool') 

    cem <- mod_ora(cem, gmt) 
    expect_equal(class(cem)[1], 'CEMiTool') 

    cem <- plot_profile(cem)
    expect_equal(class(cem)[1], 'CEMiTool') 

    cem <- plot_gsea(cem)
    expect_equal(class(cem)[1], 'CEMiTool') 
    
    cem <- plot_ora(cem)
    expect_equal(class(cem)[1], 'CEMiTool') 
})

test_that('no expression data throws error', {
    expect_error(cemitool(), 'Must*')
})
