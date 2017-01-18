library(devtools)
library(data.table)

load_all()

doParallel::registerDoParallel(cores=8)

# Load expression data
data(exprs)

# Load your sample annotation data
data(sample_annotation)

# create a new CEMiTool object
cem_obj <- new("CEMiTool", expression=exprs, 
			   sample_annotation=sample_annotation,
			   sample_name_column="Sample")

# Find the modules
cem_obj@module <- find_modules(cem_obj@expression) 

# If you desire, split the modules by positive and negative correlation
cem_obj <- split_modules(cem_obj)

# Take a look at the expression patterns of those modules.
# You can also use the gene_module variable here.
cem_obj <- plot_profile(cem_obj)
print(cem_obj@profile_plot[[1]])

# Now that you got the modules, why don't you just give a look
# at the enrichment of those modules in your experimental classes
cem_obj <- mod_gsea(cem_obj)

# Heatmap of gene set enrichment analysis
cem_obj <- plot_gsea(cem_obj)
print(pl)

# Performs over representation analysis
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_list <- read_gmt(gmt_fname)
gmt_in <- prepare_gmt(gmt_list)
ora_res <- mod_ora(splitted_modules, gmt_in)

# plot ora results
list_of_ora_results <- plot_ora(ora_res)
print(list_of_ora_results[[2]]$pl)

# running cemitool
res <- cemitool(exprs, sample_annotation, gmt_in, plot=T, split_modules=T)

# generate report
generate_report(res)
