context('Cemitool methods')

data(cem)
data(expr0)
data(sample_annot)

cem0 <- new_cem()

cem_base <- cem0
expr_data(cem_base) <- expr0 

sample_annotation(cem_base) <- sample_annot

ppi <- system.file('extdata', 'interactions.tsv', package='CEMiTool')
ppi <- fread(ppi, data.table=FALSE)

gmt <- system.file('extdata', 'pathways.gmt', package='CEMiTool')
gmt <- read_gmt(gmt)

cem_filt <- filter_expr(cem_base, 1)

cem_full <- cem

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

    cem_full <<- cem
})

test_that('CEMiTool throws errors when there is missing data', {
    expect_error(cemitool(), "Please*")
    expect_error(filter_expr(cem0), "CEMiTool object has no expression file!")
    expect_error(find_modules(cem0), "CEMiTool object has no expression file!")
    expect_warning(mod_gsea(cem0), "CEMiTool object has no expression file!")
    expect_warning(mod_gsea(cem_filt), "No modules in CEMiTool object! Did you run find_modules()?")  
    expect_warning(mod_ora(cem0, gmt), "No modules in CEMiTool object! Did you run find_modules()?")
})

test_that('Write files work as expected', {
    outputs <- c("enrichment_es.tsv", "enrichment_nes.tsv", "enrichment_padj.tsv",
                 "interactions.tsv", "module.tsv", "modules_genes.gmt",
                 "ora.tsv", "parameters.tsv", "selected_genes.txt",
                 "summary_eigengene.tsv", "summary_mean.tsv", "summary_median.tsv")
    expect_error(write_files(cem_full, force=TRUE), NA)
    expect_equal(outputs, dir('Tables'))
    unlink('Tables', recursive=TRUE)
})
