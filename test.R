library(devtools)
library(data.table)

load_all()

doParallel::registerDoParallel(cores=8)

# Load expression data
data(exprs)

# Load your sample annotation data
data(sample_annotation)

# Find the modules
gene_module <- find_modules(exprs) 

# If you desire, split the modules by positive and negative correlation
splitted_modules <- split_modules(exprs, gene_module)

# Take a look at the expression patterns of those modules.
# You can also use the gene_module variable here.
list_of_profiles <- plot_profile(exprs, splitted_modules)
print(list_of_profiles[[2]])

# Now that you got the modules, why don't you just give a look
# at the enrichment of those modules in your experimental classes
enrich <- mod_gsea(exprs, splitted_modules, sample_annotation)

# Heatmap of gene set enrichment analysis
plot_gsea(enrich)

# Performs over representation analysis
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_list <- read_gmt(gmt_fname)
gmt_in <- prepare_gmt(gmt_list)
ora_res <- mod_ora(splitted_modules, gmt_in)

# plot ora results
list_of_ora_results <- plot_ora(ora_res)
print(list_of_ora_results[[2]]$pl)

