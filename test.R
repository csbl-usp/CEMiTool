library(devtools)
library(data.table)

load_all()

# Read your expression table
exprs <- fread('data/expression.txt', data.table=F)
rownames(exprs) <- exprs[,1]
exprs[,1] <- NULL
exprs <- head(exprs[order(apply(exprs, 1, var), 
                          decreasing=T),],4000)

# Read your sample annotation table
annot <- fread('data/annotation.txt', data.table=F)

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


