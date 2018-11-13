#' @importFrom data.table setDF data.table melt
#' @import limma
#' @importFrom gRbase combnPrim
#' @importFrom igraph graph_from_edgelist graph_from_data_frame simplify degree set_vertex_attr 
#' @importFrom igraph layout.drl communities cluster_fast_greedy cluster_edge_betweenness cluster_fast_greedy 
#' @importFrom igraph cluster_label_prop cluster_leading_eigen cluster_louvain cluster_optimal cluster_spinglass cluster_walktrap
#' @import GeneOverlap
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import fgsea
#' @importFrom plyr rbind.fill
#' @importFrom utils combn
#' @importFrom WGCNA cor
#' @importFrom ff as.ffdf ff ffdforder as.ff 
#' @importFrom ffbase ffdfrbind.fill merge.ffdf subset.ffdf duplicated.ffdf ffdfdply
#' @importFrom matrixStats rowSums2 rowMedians rowMeans2 rowSds
#' @importFrom RColorBrewer brewer.pal
NULL

#' Integrates CEMiTool analyses
#'
#' Returns the occurrence of edges between different analyses
#'
#' @param analyses List of objects of class \code{CEMiTool}
#' @param num_studies The minimum number of objects in the \code{analyses} list in
#' which an edge pair must be present to be selected (default = 0)
#' @param desired_table Character string indicating the type of output to be returned.
#' Default: 'adjacency'.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#'
#' @return Object of class \code{data.frame} containing edgelist describing common 
#' edges between the networks defined in module slots from the input objects
#'
#' @details The method assumes that all genes inside each module are connected to
#' every other gene from the same module. Argument desired_table must be one of 
#' \code{spearman} (returns Spearman's rho), \code{pearson} (Pearson's R), 
#' \code{b_correlations} (returns adjacency list defined in CEMiTool object),
#' \code{adjacency} (returns discretized edges)
#'
#' @examples
#' \dontrun{ 
#' # Run the cemitool function twice on expr dataset. Each time, one sample will be removed
#' data(expr0)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1, plot=FALSE) 
#' cem2 <- cemitool(dset2, plot=FALSE) 
#' cem_overlap_df <- cem_overlap(list(cem1, cem2))
#' }
#' @rdname cem_overlap
#' @export

cem_overlap <- function(analyses, num_studies = 0, desired_table = 'adjacency', verbose=TRUE){
    if(is.null(names(analyses))){
        names(analyses) <- paste0('cem', seq_along(analyses))
    }
    study_names <- names(analyses)
    if(verbose) message("Running overlap. May take a while.")
    analyses <- Filter(Negate(is.null), analyses)
    expr_dfs <- lapply(analyses, expr_data)
    if(desired_table == "b_correlations"){
        adj_list <- lapply(analyses, adj_data)
        names(adj_list) <- study_names
    }
    # combines all genes inside each module
    genes <- lapply(seq_along(analyses), function(index){
        cem_genes <- module_genes(analyses[[index]])
        cem_name <- names(analyses[index])
        mods <- split(cem_genes[,'genes'], cem_genes[,'modules'])
        mods_log <- sapply(mods, length) < 2
        mods <- mods[!mods_log]
        mods <- mods[names(mods) != 'Not.Correlated']
    })
    names(genes) <- study_names
    names(expr_dfs) <- study_names
    rm(analyses)
    gc()
    edgelist <- lapply(names(expr_dfs), function(cem_name){
        if(verbose) message(cem_name)
        if(desired_table %in% c('pearson', 'spearman')){
            per_mod <- lapply(names(genes[[cem_name]]), function(mod_name){
                if(verbose) message("Module ", mod_name, " of object ", cem_name)
                mod_expr <- expr_dfs[[cem_name]][genes[[cem_name]][[mod_name]], ]
                mod_expr <- t(mod_expr[ order(row.names(mod_expr)), ])
                mod_expr <- WGCNA::cor(mod_expr, method=desired_table, 
                                       use="pairwise.complete.obs")
                mod_expr[upper.tri(mod_expr, diag=TRUE)] <- NA
                mod_expr <- melt(mod_expr)
                names(mod_expr) <- c("gene1", "gene2", "value")
                mod_expr <- mod_expr[!is.na(mod_expr$value), ]
                mod_expr <- ff::as.ffdf(mod_expr)
                return(mod_expr)
            })
        }else if(desired_table == "adjacency"){
            per_mod <- 
                lapply(names(genes[[cem_name]]), function(mod_name){
                    if(verbose) message("Module ", mod_name, " of object ", cem_name)
                    
                    mod_outp <- 
                        do.call(rbind, lapply(gRbase::combnPrim(genes[[cem_name]][[mod_name]], 2, simplify = FALSE), sort))
                    mod_outp <- data.frame(mod_outp, TRUE)
                    colnames(mod_outp) <- c('gene1', 'gene2', cem_name)
                    rownames(mod_outp) <- NULL
                    mod_outp <- ff::as.ffdf(mod_outp)
                    return(mod_outp)
                })
        }else if(desired_table == "b_correlations"){
            per_mod <- 
                lapply(names(genes[[cem_name]]), function(mod_name){
                    if(verbose) message("Module ", mod_name, " of object ", cem_name)
                    adj_mod <- adj_list[[index]][genes[[index]][[mod_name]], genes[[index]][[mod_name]]]
                    adj_mod[upper.tri(adj_mod, diag=TRUE)] <- NA
                    adj_mod <- melt(adj_mod)
                    adj_mod <- adj_mod[!is.na(adj_mod$value), ]
                    colnames(adj_mod) <- c('gene1', 'gene2', cem_name)
                    adj_mod <- ff::as.ffdf(adj_mod)
                    return(adj_mod)
                })
        }
        edges <- do.call("ffdfrbind.fill", per_mod)
        rownames(edges) <- 1:nrow(edges)
        return(edges)
    })
    names(edgelist) <- study_names
    
    if(num_studies > 0){
        for(index in seq_along(edgelist)){
            original_study <- rep(names(edgelist)[index], nrow(edgelist[[index]]))
            edgelist[[index]]$object <- as.ff(as.factor(original_study))
        }
        
        full_edgelist <- do.call("ffdfrbind.fill", edgelist)
        full_edgelist$splitBy <- with(full_edgelist[, 1:3], 
                                      as.ff(as.factor(paste(gene1, gene2, sep="_"))), 
                                      by = 100000)
        
        tmp <- ffdfdply(full_edgelist, split=full_edgelist$splitBy, trace=TRUE, FUN=function(y){
            res <- y %>% group_by(splitBy) %>% mutate(count=n())
        })
        
        keep_edges <- subset(tmp, count >= num_studies)
        
        new_edgelist <- list()
        for(study in study_names){
            kept_edges <- subset(keep_edges, object == study)
            rownames(kept_edges) <- NULL
            kept_edges <- kept_edges[setdiff(colnames(kept_edges), c("object", "splitBy", "count"))]
            new_edgelist[[study]] <- kept_edges
        }
        
        edgelist <- new_edgelist    
    }
    
    # Merges all studies
    #out <- Reduce(function(...){merge(..., by=c('gene1', 'gene2'), all=TRUE)}, edgelist2)
    if(verbose) message("Merging...")
    out <- Reduce(function(...){outer_join_merge(...)}, edgelist)
    if(verbose) message("Merged!")
    #data.table::setDF(out)
    colnames(out)[!colnames(out) %in% c('gene1', 'gene2')] <- study_names
    # Sum of studies containing pair and order dataframe by sum of occurrences 
    presentin <- ncol(out[, study_names]) - matrixStats::rowSums2(is.na(out[, study_names]))
    out$edgeCount <- ff::ff(presentin)
    out$proportion <- ff::ff(presentin/length(study_names))
    # Keep only edges present in at least the number 
    #of cemitool objects specified in 'num_studies' variable
    ######## out <- subset(out, edgeCount >= num_studies)
    if(desired_table %in% c('spearman', 'pearson', 'b_correlations')){
        cor_median <- matrixStats::rowMedians(as.matrix(out[, study_names]), na.rm=TRUE)
        out$edgeCorMedian <- ff::ff(cor_median)
        cor_sd <- matrixStats::rowSds(as.matrix(out[, study_names]), na.rm=TRUE)
        out$edgeSd <- ff::ff(cor_sd)
        cor_mean <- matrixStats::rowMeans2(as.matrix(out[, study_names]), na.rm=TRUE)
        out$edgeCorMean <- ff::ff(cor_mean)
        idx <- ff::ffdforder(out[c("proportion", "edgeCorMedian")], decreasing=TRUE)
    }else{
        idx <- ff::ffdforder(out[c("proportion")], decreasing=TRUE)
    }
    out <- out[idx, ]
    return(out)
}


#' Calculate full outer join for two ffdf objects
#'
#' @param x An \code{ffdf} object
#' @param y An \code{ffdf} object
#'
#' @return A \code{ffdf} object containing the join result
#' @keywords internal
#'
outer_join_merge <- function(x, y){
    # do a left outer join
    leftjoin <- merge(x, y, by = c('gene1', 'gene2'), all.x = TRUE)
    list_index <- ncol(x) - 2
    names(leftjoin) <- c("gene1", "gene2", paste0("value", seq(1, list_index + 1)))
    
    # do a right outer join (it's just a left outer join with the objects swapped)
    rightjoin <- merge(y, x, by = c('gene1', 'gene2'), all.x = TRUE, suffixes=c(".y", ".x"))
    names(rightjoin) <- c("gene1", "gene2", paste0("value", list_index + 1), paste0("value", seq(1, list_index)))
    
    if(inherits(leftjoin, "ffdf") & inherits(rightjoin, "ffdf")){
        stacked <- ffbase::ffdfrbind.fill(leftjoin, rightjoin)
    }else{
        stacked <- plyr::rbind.fill(leftjoin, rightjoin)    
    }
    # remove duplicate rows
    stacked_names <- names(stacked)
    not_columns <- stacked_names[!stacked_names %in% c("gene1", "gene2")]
    stacked_genes <- stacked[setdiff(stacked_names, not_columns)]
    
    rownames(stacked) <- NULL
    rownames(stacked_genes) <- NULL
    uniques <- stacked[!duplicated.ffdf(stacked_genes), ]
    
    return(uniques)
}


#' Generates communities from edgelist 
#'
#' Returns communities from edgelist created by cemoverlap function.
#'
#' @param mod_intersection_df Module intersection dataframe obtained from cemoverlap function
#' @param presence_as_weights Logical. Should sums of node pair occurence be considered edge weights?
#' @param smallest_community Minimal number of genes in community (default:15)
#' @param method Character string denoting an \code{igraph} package clustering function. 
#' One of 'cluster_fast_greedy',  'cluster_edge_betweenness', 'cluster_fast_greedy',
#'  'cluster_label_prop', 'cluster_leading_eigen', 'cluster_louvain', 'cluster_optimal', 
#'  'cluster_spinglass' or 'cluster_walktrap'. Default: 'cluster_fast_greedy'
#'
#' @details Function takes edgelist as inputs and generates communities using functions 
#'    provided in igraph package (default:'cluster_fast_greedy')
#' 
#' @return A list containing the genes present in each community detected
#'
#' @export 
#' @examples
#' \dontrun{ 
#' # Run the cemitool function twice on expr dataset. Each time, one sample will be removed
#' data(expr0)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1, plot=FALSE) 
#' cem2 <- cemitool(dset2, plot=FALSE) 
#' cem_overlap_df <- cem_overlap(list(cem1, cem2))
#' comm_overlap_df <- overlap_community(cem_overlap_df)
#' }
overlap_community <- function(mod_intersection_df, presence_as_weights = FALSE,
                              smallest_community = 15, 
                              method = c('cluster_fast_greedy', 'cluster_edge_betweenness', 
                                         'cluster_label_prop', 'cluster_leading_eigen',
                                         'cluster_louvain', 'cluster_optimal', 
                                         'cluster_spinglass', 'cluster_walktrap')){
    
    method <- match.arg(method)
    edgemat <- as.matrix(mod_intersection_df[,c('gene1', 'gene2')])
    edgegraph <- igraph::graph_from_edgelist(edgemat, directed = FALSE)
    if(presence_as_weights){
        weight <- mod_intersection_df[, 'proportion']
        edgegraph$weight <- weight
    }
    commfunc <- get(method)
    comm <- commfunc(edgegraph)
    comm <- igraph::communities(comm)
    comm <- as.list(comm)
    names(comm) <- paste0('CM', seq_along(1:length(comm)))
    len_vec <- sapply(comm, length)
    names(comm) <- ifelse(len_vec >= smallest_community, 
                          names(comm), 
                          paste0(names(comm), '.SMALL'))
    comm <- comm[order(len_vec, decreasing = TRUE)]
    # comm <- comm[sapply(comm, length) >= smallest_community]
    return(comm)
}

#' Enriches communities 
#'
#' Returns enrichment of communities from edgelist created by cemoverlap function.
#'
#' @param community_list Community list obtained from overlap community function
#' @param analyses List of CEMiTool objects
#' @param comp_group Which group will be used as base for comparison. If 'none', 
#' then all combinations of comparisons will be made
#' @param subject_col Column containing subject information in sample annotation slot of CEMiTool objects 
#' @param run_fgsea Logical. Should fgsea be run?
#'
#' @details This function assumes that relevant modules for a comparison in a study will 
#' have a high proportion of differentially regulated genes in a certain direction. Base assumption is that
#' NON-relevant modules will be centered at zero.
#' 
#' @return A \code{data.frame} containing information of how much each comparison is enriched in each community 
#' in each CEMiTool object.
#' 
#' @examples 
#' \dontrun{
#' # Run the cemitool function twice on expr dataset. Each time, one sample will be removed
#' data(expr0)
#' data(sample_annot)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1, plot=FALSE) 
#' cem2 <- cemitool(dset2, plot=FALSE) 
#' cem_overlap_df <- cem_overlap(list(cem1, cem2))
#' 
#' comm_overlap <- overlap_community(cem_overlap_df)
#' 
#' samples1 <- names(expr_data(cem1))
#' samples2 <- names(expr_data(cem2))
#' 
#' sample_annotation(cem1) <- sample_annot[sample_annot$SampleName %in% samples1, ]
#' sample_annotation(cem2) <- sample_annot[sample_annot$SampleName %in% samples2, ]
#' 
#' mod_enrich <- enrich_mods(comm_overlap, list(cem1, cem2), comp_group='g0')
#' }
#' @rdname enrich_mods
#' @export 
enrich_mods <- function(community_list, analyses, 
                        comp_group = 'none', subject_col=NULL, 
                        run_fgsea = FALSE){
    if(is.null(names(analyses))){
        names(analyses) <- paste0('cem', seq_along(analyses))
    }
    
    # Setting
    mod_exp <- Map(function(cem, cem_name){
        annot <- sample_annotation(cem)
        annot <- annot[annot[[cem@sample_name_column]] %in% names(expr_data(cem)), ]
        contmat <- makeContMatrix(sample_annot = annot, class_column = cem@class_column,
                                  comp_group = comp_group, expr = expr_data(cem), 
                                  subject_col = subject_col) 
        toptables <- do.call(makeLimmaComp, contmat)
        iter_comp <- Map(function(comp, compname){
            if(run_fgsea){
                ranks <- rank(comp$logFC)
                names(ranks) <- comp$gene
                fgseaRes <- fgsea::fgsea(pathways = community_list, stats = ranks,
                                         minSize = 15, maxSize = 500, nperm = 1000)
                fgseaRes <- as.data.frame(fgseaRes)
                fgseaRes <- fgseaRes[order(fgseaRes$padj),]
            }else{
                fgseaRes <- data.frame()
            }
            percentage_in_comparison <- Map(function(mod, modname){
                comp <- subset(comp, gene %in% mod & P.Value < 0.05)
                per_up <- sum(comp$logFC > 0)/length(mod)
                per_do <- sum(comp$logFC < 0)/length(mod)
                fcs_pert <- (median(comp$logFC)/sd(comp$logFC)) * median(-log10(comp$P.Value))
                matched_gsea <- match(modname, fgseaRes$pathway)
                if(!is.na(matched_gsea)){
                    gsea_vec <- as.numeric(fgseaRes[matched_gsea, 2:7])
                }else{ 
                    gsea_vec <- rep(NA, 6)
                }
                scaled <- per_up - per_do
                scaled_percentage <- c(per_up, per_do, scaled, fcs_pert, compname, modname, 
                                       cem_name, paste0(compname, '_in_', cem_name), gsea_vec)
                return(scaled_percentage)
            }, community_list, names(community_list))
            percentage_in_comparison <- as.data.frame(do.call(rbind, percentage_in_comparison))
            colnames(percentage_in_comparison) <- c('percentage_up', 'percentage_down', 
                                                    'scaled_percentage', 'score_mod',
                                                    'comparison', 'community', 'cem_obj', 
                                                    'name_comp', 'pval', 'padj', 'ES', 'NES', 
                                                    'nMoreExtreme', 'size')
            return(percentage_in_comparison)
        }, toptables, names(toptables))
        p_incomp_df <- do.call(rbind, iter_comp)
        return(p_incomp_df)
    }, analyses, names(analyses))
    mod_exp_df <- as.data.frame(do.call(rbind, mod_exp))
    ncolumns <- c('percentage_up', 'percentage_down', 'scaled_percentage', 'score_mod',
                  'pval', 'padj', 'ES', 'NES', 'nMoreExtreme', 'size')
    for(i in ncolumns){
        if(i %in% colnames(mod_exp_df)){
            mod_exp_df[[i]] <- as.numeric(mod_exp_df[[i]])  
        }
    }
    mod_exp_df <- Filter(function(x) !all(is.na(x)), mod_exp_df)
    #rownames(mod_exp_df) <- NULL
    return(mod_exp_df)
}

#' Make a LIMMA contrast matrix
#'
#' @param sample_annot Sample annotation \code{data.frame}
#' @param class_column Character string indicating the sample 
#' grouping column. Default: "Class"
#' @param which_groups Optional character vector to use if selecting groups 
#' to be compared (e.g. disease1, disease3, but not disease2)
#' @param comp_group Optional character string indicating group to compare groups against. 
#' If "none" (the default), compares all groups against each other. 
#' @param subject_col Optional character string indicating a column containing 
#' subject information 
#' @param expr Optional expression matrix for direct compatibility with makeLimmaComp 
#'
#' @return A list containing LIMMA contrast and design matrices
#' @keywords internal 
#'
#' @examples
#' \dontrun{
#' # Create a mock expressionset to test functions 
#' mockexpset <- matrix(rnorm(15000), ncol = 15)
#' colnames(mockexpset) <- paste0('GSM', seq(1, 15))
#' set.seed(100)
#' rownames(mockexpset) <- apply(replicate(n = 1000, sample(letters, 8)),
#'     2, function(x) paste(x, collapse = '')) 
#' sample_annot <- data.frame(geo_accession = colnames(mockexpset), 
#'                   group = c(rep('skin_healthy', 5), rep('degree1', 5), rep('degree2', 5)), 
#'                   subject = paste0('S', c(rep(1:5),rep(1:5), rep(1:5))))
#' # Testing function below 
#' # Paired
#' cont_mat1 <- makeContMatrix(class_column = 'group', which_groups = c('degree1', 'degree2'), 
#'                  comp_group = 'skin_healthy', subject_col= 'subject', 
#'                  sample_annot = sample_annot, expr = mockexpset)
#' limma_result1 <- do.call(makeLimmaComp, cont_mat1)
#' # Unpaired
#' cont_mat2 <- makeContMatrix(class_column = 'group', which_groups = c('degree1', 'degree2'), 
#'                  comp_group = 'skin_healthy', sample_annot = sample_annot, expr = mockexpset)
#' limma_result2 <- do.call(makeLimmaComp, cont_mat2)
#' }
makeContMatrix <- function(sample_annot, class_column = 'Class', which_groups = 'all',
                           comp_group = 'none', subject_col=NULL, expr){
    
    groups <- as.character(unique(sample_annot[[class_column]]))
    stop_if(is.null(groups), 
                       "No groups in the class column of the sample annotation data")
    
    if(which_groups == 'all' && comp_group == 'none'){
        comparisons <- utils::combn(groups, 2, simplify = FALSE)
    }else if(which_groups == 'all' && comp_group != 'none'){
        if(!comp_group %in% groups){
            stop('Comparison group is not present in groups')
        }
        groups <- groups[!groups %in% comp_group]
        comparisons <- lapply(groups, c, comp_group)
    }else if(which_groups != 'all' && comp_group == 'none'){
        comparisons <- utils::combn(which_groups, 2, simplify = FALSE)
    }else if(which_groups != 'all' && comp_group != 'none'){
        comparisons <- as.list(paste(which_groups, comp_group))
        comparisons <- sapply(comparisons, 
                              function(x) strsplit(x, ' ', fixed = TRUE))
    }
    design_levels <- factor(sample_annot[[class_column]])
    if(!is.null(subject_col)){
        subject_levels <- factor(sample_annot[[subject_col]])
        design <- model.matrix(~0 + design_levels + subject_levels) 
        colnames(design) <- gsub('design_levels', '', colnames(design))
        colnames(design) <- gsub('subject_levels', '', colnames(design))
    }else{
        design <- model.matrix(~0 + design_levels)
        colnames(design) <- gsub('design_levels', '', colnames(design))
    }
    comp <- lapply(comparisons, function(x) 
        paste0(x[1],'_vs_',x[2],'=',x[1],'-',x[2]))
    comp$levels <- design
    contrasts_names <- lapply(comparisons, function(x)
        paste0(x[1], '_vs_',x[2]))
    cont_matrix <- do.call(limma::makeContrasts, comp)
    colnames(cont_matrix) <- contrasts_names
    if(missing(expr)){
        warning('An expression matrix must be submitted for direct compatibility with makeLimmaComp() function. Alternatively, a matrix can be added to the list object returned by this function.')
        ret_list <- list(design = design, cont.matrix = cont_matrix)
    }else{
        ret_list <- list(design = design, cont.matrix = cont_matrix, expr = expr)
    }
    return(ret_list)
}

#' Make LIMMA comparisons
#'
#' @param expr Expression
#' @param design Design matrix
#' @param cont.matrix Contrast matrix
#'
#' @return A list with one data.frame per comparison
#' @keywords internal 
#' @examples 
#' \dontrun{
#' # Create a mock expressionset to test functions 
#'       mockexpset <- matrix(rnorm(15000), ncol = 15)
#'       colnames(mockexpset) <- paste0('GSM', seq(1, 15))
#'       set.seed(100)
#'       rownames(mockexpset) <- apply(replicate(n = 1000, sample(letters, 8)),
#'                                   2, function(x) paste(x, collapse = '')) 
#'       sample_annot <- data.frame(geo_accession = colnames(mockexpset), 
#'                                group = c(rep('skin_healthy', 5), rep('degree1', 5), 
#'                                          rep('degree2', 5)), 
#'                                subject = paste0('S', c(rep(1:5),rep(1:5), rep(1:5))))
#' # Testing function below 
#'       # Paired
#'       cont_mat1 <- makeContMatrix(class_column = 'group', which_groups = c('degree1', 'degree2'), 
#'                                  comp_group = 'skin_healthy', subject_col= 'subject', 
#'                                  sample_annot = sample_annot, expr = mockexpset)
#'       limma_result1 <- do.call(makeLimmaComp, cont_mat1)
#'       # Unpaired
#'       cont_mat2 <- makeContMatrix(class_column = 'group', which_groups = c('degree1', 'degree2'), 
#'                                  comp_group = 'skin_healthy', sample_annot = sample_annot, 
#'                                  expr = mockexpset)
#'       limma_result2 <- do.call(makeLimmaComp, cont_mat2)
#' }
makeLimmaComp <- function(expr, design, cont.matrix){
    fit <- lmFit(expr, design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    coeffs <- colnames(cont.matrix)
    tops <- lapply(coeffs, function(x) 
        topTable(fit2, number = Inf, coef = x, adjust.method = 'BH'))
    names(tops) <- coeffs
    tops <- lapply(tops, function(x) cbind(gene = rownames(x), x))
    tops <- lapply(tops, function(x) { rownames(x) <- NULL; x})
    nmes <- names(tops)
    toptbs <- lapply(setNames(nmes, nmes), function(x) {cbind(tops[[x]], comparison = rep(x, nrow(tops[[x]])))})
    return(toptbs)
}


#' Plot module edge co-membership
#'
#' @param cem_overlap_df Output of function \code{cem_overlap} 
#'
#' @return A plot containing the number of edges present in each fraction of studies
#' @export
#'
#' @examples
#' \dontrun{ 
#' # Run the cemitool function five times on expr0 dataset. Each time, 10 samples will be removed.
#' data(expr0)
#' data(sample_annot)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(11)
#' dset2 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(12)
#' dset3 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(13)
#' dset4 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(14)
#' dset5 <- expr0[,-sample(1:ncol(expr0), 10)]
#' 
#' cem1 <- cemitool(dset1, sample_annot, plot=FALSE) 
#' cem2 <- cemitool(dset2, sample_annot, plot=FALSE)
#' cem3 <- cemitool(dset3, sample_annot, plot=FALSE) 
#' cem4 <- cemitool(dset4, sample_annot, plot=FALSE) 
#' cem5 <- cemitool(dset5, sample_annot, plot=FALSE) 
#' 
#' cem_overlap_df <- cem_overlap(list(cem1, cem2, cem3, cem4, cem5))
#' plot_comembership(cem_overlap_df)
#' }
plot_comembership <- function(cem_overlap_df){
    
    x <- table(cem_overlap_df$edgeCount[])
    x <- as.data.frame(x)
    x$Var1 <- as.factor(round(as.numeric(as.character(x$Var1)), digits = 3))
    
    ggplot(x, aes(x=Var1, y=Freq, group=1)) +
        geom_line(color="darkgrey") +
        geom_text(aes(label=Freq)) +
        labs(x="Number of CEMiTool objects", y="Number of edges", title="Edge co-membership") +
        theme(axis.text=element_text(size=12), plot.title=element_text(hjust=0.5)) +
        scale_x_discrete(limits=levels(x$Var1))
}

#' Plot graph of consensus modules
#'
#' @param cem_overlap_df Output \code{data.frame} from function 
#' \code{cem_overlap}
#' @param comm_overlap_df Output \code{data.frame} from function 
#' \code{overlap_community}
#' @param study_num Minimum number of studies an edge must be present in for it 
#' to be included
#' @param num_sd_cut Number of standard deviations an edge's mean must be above in order
#' to be included in the final plot (Default: 2)
#'
#' @return A plot containing the consensus modules of the studies
#' @export
#'
#' @examples
#' \dontrun{
#' # Run the cemitool function twice on expr dataset. Each time, one sample will be removed
#' data(expr0)
#' data(sample_annot)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1, sample_annot, plot=FALSE) 
#' cem2 <- cemitool(dset2, sample_annot, plot=FALSE) 
#' cem_overlap_df <- cem_overlap(list(cem1, cem2))
#' comm_overlap_df <- overlap_community(cem_overlap_df)
#' plot_consensus(cem_overlap_df, comm_overlap_df, study_num=2)
#' }
plot_consensus <- function(cem_overlap_df, comm_overlap_df, study_num, num_sd_cut=2){
    overlap_df <- cem_overlap_df
    
    cem_num <- grep("cem", names(overlap_df), value=TRUE)
    
    stop_if(study_num > length(cem_num), "Variable study_num cannot exceed number of studies")
    stop_if(study_num == 1, "Cannot do consensus of one study")
    overlap_df <- overlap_df[overlap_df$edgeCount >= study_num, ]
    
    if("edgeCorMean" %in% names(overlap_df)){
        means <- overlap_df$edgeCorMean[]
        overlap_df <- overlap_df[abs(means) >= num_sd_cut*sd(means), ]    
    }
    
    ig_obj <- igraph::graph_from_data_frame(overlap_df, directed=FALSE)
    ig_obj <- igraph::simplify(ig_obj)
    degrees <- igraph::degree(ig_obj, normalized=FALSE)
    ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
    
    plotcord <- data.frame(igraph::layout.drl(ig_obj, options=list(simmer.attraction=0, 
                                                           simmer.temperature=400)))
    
    net_obj <- intergraph::asNetwork(ig_obj)
    edglist <- network::as.matrix.network.edgelist(net_obj)
    
    edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
    colnames(edges) <- c("X1","Y1","X2","Y2")
    
    plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
    plotcord$Names <- network::get.vertex.attribute(net_obj, "vertex.names")
    
    comm_list <- lapply(names(comm_overlap_df), function(x){
        genes <- comm_overlap_df[[x]]
        names(genes) <- rep(x, length(genes))
        genes
    })
    comms <- unlist(comm_list)
    comms <- flip_vector(comms)
    
    plotcord$Communities <- comms[plotcord$Names]
    
    mycolors <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(30) 
    
    ggplot() + 
        geom_segment(aes(x=X1, y=Y1, xend=X2, yend=Y2), 
                     data=edges, size = 0.5, alpha=0.5, colour="#DDDDDD") + 
        geom_point(aes(X1, X2, size=plotcord$Degree, alpha=0.9, 
                       color=Communities), data=plotcord) + 
        scale_color_manual(values=mycolors) +
        scale_size_continuous(name="Degree") +
        ggplot2::theme_bw(base_size = 12, base_family = "") +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white",
                                                                colour = NA),
                       panel.border = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank()) + 
        guides(alpha=FALSE)
}

#' Calculate module overlap statistics
#'
#' @param analyses List of \code{CEMiTool} objects
#' @param comp_group Character string indicating the group to be compared against
#' @param subject_col Optional character string indicating a column in the \code{CEMiTool} 
#' objects' sample annotation slot object containing subject information 
#' @param ... Additional parameters
#'
#' @return A list containing overlap statistics, node and study information and module activity
#' @export
#'
#' @examples
#' \dontrun{
#' # Run the cemitool function twice on expr dataset. Each time, one sample will be removed
#' data(expr0)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1, sample_annot, plot=FALSE) 
#' cem2 <- cemitool(dset2, sample_annot, plot=FALSE)
#' mod_stats <- stat_overlap_mods(analyses=list(cem1, cem2), comp_group="g0")
#' }
stat_overlap_mods <- function(analyses, comp_group, subject_col=NULL, ...){
    if(is.null(names(analyses))){
        names(analyses) <- paste0('cem', seq_along(analyses)) 
    }
    stopifnot(sum(duplicated(names(analyses))) == 0)
    lapply(names(analyses), function(name) {
        if(length(gsea_data(analyses[[name]])) == 0){
            stop("CEMiTool object '", name, 
                 "' has no enrichment slot! Please run function mod_gsea.")
        }
    })
    
    df_output <- mod_compare(analyses, ...)
    
    # Compute node information and study information
    info_mod <- mod_info(analyses, df_output)
    
    # Determine module activity based on limma
    mod_mean <- mod_activity(analyses, comp_group=comp_group, subject_col=subject_col)
    
    # Return three dataframes with aggregation of modules.
    meta_info <- list(module_comparison = df_output, metric_df = info_mod, module_info = mod_mean)
    return(meta_info)
}

#' Overlap statistics
#'
#' @param analyses List of \code{CEMiTool} objects
#' @param p_thresh p-value threshold
#' @param fdr_thresh Threshold for FDR 
#' @param jac_thresh Threshold for Jaccard
#'
#' @return A \code{data.frame} containing the results
#' @keywords internal
#'
mod_compare <- function(analyses, p_thresh = 1, fdr_thresh = 1, jac_thresh = 0){
    names(analyses) <- paste0("cem_", seq_along(analyses))
    new_mods_per_cem <- lapply(analyses, function(cem){
        tmpmod <- subset(module_genes(cem), modules != 'Not.Correlated')
        spmod <- split(tmpmod$genes, tmpmod$modules)
        spmod
    })
    new_mods <- unlist(new_mods_per_cem, recursive = FALSE)
    universe <- length(unique(Reduce(union, new_mods)))
    first_lev <- Map(function(fmod_name, firstlev){
        second_lev <- Map(function(smod_name, seclev){
            gene_ovlp <- GeneOverlap::newGeneOverlap(firstlev, seclev, genome.size = universe)
            gene_ovlp <- GeneOverlap::testGeneOverlap(gene_ovlp)
            jac <- slot(gene_ovlp, 'Jaccard')
            fis <- slot(gene_ovlp, 'pval')
            # fmod_len <- paste0(fmod_name, '//', length(firstlev))
            # smod_len <- paste0(smod_name, '//', length(seclev))
            mod_ord <- sort(c(fmod_name, smod_name)) 
            #df_row <- c(mod_ord, jac, fis)
            df_row <- data.frame(mod_ord[1], mod_ord[2], jac, fis, stringsAsFactors = FALSE)
            return(df_row)
        }, names(new_mods), new_mods)
        do.call(rbind, second_lev)
    }, names(new_mods), new_mods)
    df_output <- as.data.frame(do.call(rbind, first_lev))
    rownames(df_output) <- NULL
    colnames(df_output) <- c('first_mod', 'second_mod', 'Jaccard', 'Fisherp')
    
    # Remove interactions from the same modules, duplicates and filter out significant interactions
    df_output <- subset(df_output, first_mod != second_mod)
    remove_pair <- duplicated(paste(df_output$first_mod, df_output$second_mod, sep = '.'))
    df_output <- df_output[!remove_pair,]
    df_output$Jaccard <- as.numeric(df_output$Jaccard)
    df_output$Fisherp <- as.numeric(df_output$Fisherp)
    df_output$fdr <- p.adjust(df_output$Fisherp)
    df_output$logp <- -log10(df_output$Fisherp) 
    df_output$logfdr <- -log10(df_output$fdr) 
    df_output <- df_output[order(df_output$Jaccard, decreasing = TRUE),]
    df_output <- subset(df_output, Jaccard > jac_thresh & Fisherp < p_thresh & fdr < fdr_thresh)
    return(df_output)
}

#' Compute node information and study information
#'
#' @param analyses List of \code{CEMiTool} objects
#' @param df_output Output of function \code{mod_compare}
#' @param gsea_metric GSEA metric to be used. One of "nes", 
#' "es" or "pval". Default: "nes".
#'
#' @return A \code{data.frame} containing the results
#' @keywords internal
mod_info <- function(analyses, df_output, gsea_metric="nes"){
    if(is.null(names(analyses))){
        names(analyses) <- paste0("cem_", seq_along(analyses))    
    }
    info_mod <- Map(function(cem_name, cem){
        # tmpmod <- subset(cem@module, modules != 'Not.Correlated')
        scores <- Map(function(scrname, scr){
            scr <- scr %>%
                gather(key = class, value = value, -pathway) %>%
                mutate(cem_name = cem_name) %>%
                mutate(metric = scrname) %>%
                filter(pathway != 'Not.Correlated') %>%
                unite(module, cem_name, pathway, sep = '.') %>%
                select(module, class, value, metric) %>%
                as.data.frame()
            return(scr)}, names(cem@enrichment), cem@enrichment)
        scores <- do.call(rbind, c(scores, make.row.names = FALSE))
        # Subset one metric. Loop was kept in case more metrics are needed.
        scores <- subset(scores, metric == gsea_metric)
        scores <- scores[,!colnames(scores) == 'metric']
        cemdf <- cem@module %>%
            filter(modules != 'Not.Correlated') %>%
            group_by(modules) %>%
            summarise(mod_length = n()) %>%
            ungroup() %>%
            mutate(cem_name = cem_name) %>%
            unite(module, cem_name, modules, sep = '.')
        cemdf <- merge(cemdf, scores, by = 'module')
        return(cemdf)
    }, names(analyses), analyses)
    info_mod <- do.call(rbind, c(info_mod, make.row.names = FALSE))
    info_mod <- info_mod %>%
        spread(key = class, value = value)%>%
        filter(module %in% c(df_output$first_mod, df_output$second_mod))
    info_mod[is.na(info_mod)] <- 0
}

#' # Determine module activity based on limma
#'
#' @param analyses List of \code{CEMiTool} objects
#' @param comp_group Character string indicating the group to be compared against
#' @param subject_col Optional character string indicating a column in the \code{CEMiTool} objects'
#' sample annotation slot object containing subject information
#'
#' @return A \code{data.frame} containing the results
#' @keywords internal
mod_activity <- function(analyses, comp_group, subject_col){
    mod_mean <- Map(function(cem_name, cem){
        # Group and subject are hardcoded
        annot <- sample_annotation(cem)
        annot <- annot[annot[[cem@sample_name_column]] %in% names(expr_data(cem)), ]
        contmat <- makeContMatrix(sample_annot = annot, 
                                  class_column = cem@class_column,
                                  comp_group = comp_group, expr = expr_data(cem),
                                  subject_col = subject_col) 
        toptables <- do.call(makeLimmaComp, contmat)
        toptables <- do.call(rbind, c(toptables, make.row.names = FALSE))
        # Now determine activity of module
        cemsp <- cem@module %>%
            filter(modules != 'Not.Correlated')
        cemsp <- split(cemsp$genes, cemsp$modules)
        cem_actv <- Map(function(modname, mod){
            tmptop <- toptables %>%
                filter(gene %in% mod) %>%
                group_by(comparison) %>%
                summarise(fc_median = median(logFC)) %>%
                ungroup() %>%
                mutate(module = paste0(cem_name, '.', modname)) %>%
                select(comparison, module, fc_median) %>%
                gather(key = parameter, value = value, fc_median)
            tmptop 
        }, names(cemsp), cemsp)
        cem_actv <- do.call(rbind, c(cem_actv, make.row.names = FALSE))
        cem_actv <- cem_actv %>%
            unite(new_col, comparison, parameter, sep = '.') %>%
            spread(new_col, value)
        cem_actv
    }, names(analyses), analyses)
    mod_mean <- plyr::rbind.fill(mod_mean)
    return(mod_mean)
}

#' Plot study module similarity
#'
#' @param mod_stats List output from function \code{stat_overlap_mods}
#' @param weight_col Character string denoting the weighting column for module similarity.
#' One of "Jaccard", "Fisherp", "fdr", "logp" or "logfdr". Default: "logfdr". 
#'
#' @return A plot showing the similarity between study modules
#' @export
#'
#' @examples
#' \dontrun{
#' # Run the cemitool function five times on expr0 dataset. Each time, 10 samples will be removed.
#' data(expr0)
#' data(sample_annot)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(11)
#' dset2 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(12)
#' dset3 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(13)
#' dset4 <- expr0[,-sample(1:ncol(expr0), 10)]
#' set.seed(14)
#' dset5 <- expr0[,-sample(1:ncol(expr0), 10)]
#' 
#' cem1 <- cemitool(dset1, sample_annot, plot=FALSE) 
#' cem2 <- cemitool(dset2, sample_annot, plot=FALSE)
#' cem3 <- cemitool(dset3, sample_annot, plot=FALSE) 
#' cem4 <- cemitool(dset4, sample_annot, plot=FALSE) 
#' cem5 <- cemitool(dset5, sample_annot, plot=FALSE) 
#' mod_stats <- stat_overlap_mods(list(cem1, cem2, cem3, cem4, cem5), comp_group="g0")
#' plot_similarity(mod_stats)
#' }
plot_similarity <- function(mod_stats, weight_col="logfdr"){
    df_output <- mod_stats[[1]]
    ig_obj <- igraph::graph_from_data_frame(df_output, directed=FALSE)
    ig_obj <- igraph::simplify(ig_obj)
    degrees <- igraph::degree(ig_obj, normalized=FALSE)
    ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
    ig_obj <- igraph::set.edge.attribute(ig_obj, "weight", value=df_output[[weight_col]])
    
    plotcord <- data.frame(igraph::layout.drl(ig_obj, options=list(simmer.attraction=0, 
                                                                   simmer.temperature=400)))
    net_obj <- intergraph::asNetwork(ig_obj)
    edglist <- network::as.matrix.network.edgelist(net_obj)
    
    edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
    edges$Weight <- network::get.edge.attribute(net_obj, "weight")
    colnames(edges) <- c("X1","Y1","X2","Y2", "Weight")
    
    plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
    plotcord$Names <- network::get.vertex.attribute(net_obj, "vertex.names")
    
    ggplot() + 
        geom_segment(aes(x=X1, y=Y1, xend=X2, yend=Y2, size = Weight), 
                     data=edges, alpha=0.5, colour="#DDDDDD") + 
        geom_point(aes(X1, X2, size=plotcord$Degree, alpha=0.9), data=plotcord) + 
        geom_text(aes(x=X1, y=X2, label=Names), hjust=0, vjust=0, data=plotcord) +
        scale_color_brewer(palette="Set1") +
        ggplot2::theme_bw(base_size = 12, base_family = "") +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white",
                                                                colour = NA),
                       panel.border = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank()) + 
        guides(alpha=FALSE)
}

