library(devtools)
library(data.table)

load_all()

registerDoParallel(cores=8)

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

# Now that you got the modules, why don't you just give a look
# at the enrichment of those modules in your experimental classes
enrich <- mod_gsea(exprs, splitted_modules, annot)

# Heatmap of gene set enrichment analysis
plot_gsea(enrich)

# Performs over representation analysis
ora_res <- mod_ora(splitted_modules, "data/pathways.gmt")
