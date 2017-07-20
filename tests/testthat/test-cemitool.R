context('Cemitool methods')

data(cem)
data(expr)
data(sample_annot)

cem0 <- new_cem()

cem_base <- cem0
expr_data(cem_base) <- expr 

sample_annotation(cem_base) <- sample_annot

ppi <- system.file('extdata', 'interactions.tsv', package='CEMiTool')
ppi <- fread(ppi, data.table=FALSE)

gmt <- system.file('extdata', 'pathways.gmt', package='CEMiTool')
gmt <- clusterProfiler::read.gmt(gmt)

cem_filt <- filter_expr(cem_base, 1)

test_that('all methods of signature CEMiTool returns CEMiTool objects', {
    expect_is(cem_filt, 'CEMiTool') 
    
    expect_is(cem, 'CEMiTool') 
   
    interactions_data(cem) <- ppi
    expect_is(cem, 'CEMiTool') 

    cem <- mod_gsea(cem)
    expect_is(cem, 'CEMiTool') 

    cem <- mod_ora(cem, gmt) 
    expect_is(cem, 'CEMiTool') 

    cem <- plot_profile(cem)
    expect_is(cem, 'CEMiTool') 

    cem <- plot_gsea(cem)
    expect_is(cem, 'CEMiTool') 
    
    cem <- plot_ora(cem)
    expect_is(cem, 'CEMiTool') 

    cem <- plot_beta_r2(cem)
    expect_is(cem, 'CEMiTool') 

    cem <- plot_mean_k(cem)
    expect_is(cem, 'CEMiTool') 
})

test_that('expression file is getting the correct number of genes post-filtering',{
    expect_equal(nrow(expr_data(cem_filt)), length(cem_filt@selected_genes))
})

test_that('CEMiTool throws errors when there is missing data', {
    expect_error(cemitool(), "Please*")
    expect_error(filter_expr(cem0), "CEMiTool object has no expression file!")
    expect_error(find_modules(cem0), "CEMiTool object has no expression file!")
    expect_warning(mod_gsea(cem0), "CEMiTool object has no expression file!")
    expect_warning(mod_gsea(cem_filt), "No modules in CEMiTool object! Did you run find_modules()?")  
    expect_warning(mod_ora(cem0, gmt), "No modules in CEMiTool object! Did you run find_modules()?")
})

