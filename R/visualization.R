#' @import ggplot2
#' @importFrom sna gplot.layout.fruchtermanreingold
#' @importFrom data.table melt
#' @importFrom ggrepel geom_label_repel
#' @importFrom igraph degree set_vertex_attr
#' @import intergraph
#' @importFrom scales squish
#' @import stringr
#' @importFrom network as.matrix.network.adjacency as.matrix.network.edgelist get.vertex.attribute
#' @import grid
NULL

#' Expression profile visualization
#'
#' Creates a plot with module gene expression profiles along samples
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param order_by_class Logical. Only used if a sample
#'           annotation file is present. Whether or not to order by the
#'           class column in the sample annotation file (as defined by the
#'           class_column slot in \code{cem}).
#' @param center_func Character string indicating the centrality measure to show in
#'        the plot. Either 'mean' (the default) or 'median'.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with profile plots
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot module gene expression profiles
#' cem <- plot_profile(cem)
#' # Check resulting plot
#' show_plot(cem, "profile")
#'
#' @rdname plot_profile
#' @export
setGeneric('plot_profile', function(cem, ...) {
    standardGeneric('plot_profile')
})

#' @rdname plot_profile
setMethod('plot_profile', signature('CEMiTool'),
    function(cem, order_by_class=TRUE, center_func='mean') {
        if(!tolower(center_func) %in% c("mean", "median")){
            stop("Invalid center_func type. Valid values are 'mean' and 'median'")
        }
        modules <- unique(cem@module$modules)
        if(is.null(modules)){
               stop("No modules in this CEMiTool object.")
        }
        #vars <- mget(ls())
        #vars$modules <- NULL
        #cem <- get_args(cem=cem, vars=vars)

        modules <- modules[order(as.numeric(stringr::str_extract(modules, "\\d+")))]
        expr <- expr_data(cem, filter=cem@parameters$filter,
                          apply_vst=cem@parameters$apply_vst)
        annot <- sample_annotation(cem)
        sample_name_column <- cem@sample_name_column
        class_column <- cem@class_column
        mod_cols <- mod_colors(cem)
        plots <- lapply(modules, function(mod){
            # subsets from expr all genes inside module mod
            genes <- cem@module[cem@module[,'modules']==mod, 'genes']
            expr[, 'id'] <- rownames(expr)
            mod_expr <- data.table::melt(expr[genes,], 'id',
                              variable.name='sample',
                              value.name='expression')

            # initialize plot base layer
            g <- ggplot(mod_expr, aes_(x=~sample, y=~expression))

            # adds different background colours if annot is provided
            if (nrow(annot)!=0) {
                if (order_by_class) {
                    # sorts data.frame by class name
                    annot <- annot[order(annot[, class_column]),]
                }
                annot[, sample_name_column] <- factor(annot[, sample_name_column],
                                                      levels=annot[, sample_name_column])
                mod_expr[, 'sample'] <- factor(mod_expr[, 'sample'],
                                               levels=annot[, sample_name_column])

                # y positioning of background tiles
                y_pos <- mean(mod_expr[, 'expression'])

                # reinitialize base layer adding background tiles
                g <- ggplot(mod_expr, aes_(x=~sample, y=~expression)) +
                     geom_tile(data=annot, alpha=0.3, height=Inf,
                               aes(x=get(sample_name_column), y=y_pos,
                               fill=as.factor(get(class_column))))
            }

            # adding lines
            g <- g + geom_line(aes_(group=~id), alpha=0.2, colour=mod_cols[mod]) +
                stat_summary(aes(group=1), size=1, fun.y=get(tolower(center_func)), geom='line')

            # custom theme
            g <- g + theme(plot.title=element_text(lineheight=0.8,
                                                   face='bold',
                                                   colour='black',
                                                   size=15),
                           axis.title=element_text(face='bold',
                                                   colour='black',
                                                   size=15),
                           axis.text.y=element_text(angle=0,
                                                    vjust=0.5,
                                                    size=8),
                           axis.text.x=element_text(angle=90,
                                                    vjust=0.5,
                                                    size=6),
                           panel.grid=element_blank(),
                           legend.title=element_blank(),
                           legend.text=element_text(size = 8),
                           legend.background=element_rect(fill='gray90',
                                                          size=0.5,
                                                          linetype='dotted'),
                           legend.position='bottom'
                           )
            # title
            g <- g + ggtitle(mod)

            return(g)
        })
        names(plots) <- modules
        cem@profile_plot <- plots
        return(cem)
    })

#' ORA visualization
#'
#' Creates a bar plot with the results of module overrepresentation analysis
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param n number of modules to show
#' @param pv_cut p-value significance cutoff. Default is 0.05.
#' @param ... parameters to plot_ora_single
#'
#' @return Object of class \code{CEMiTool} with ORA plots
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Read example gmt file
#' gmt <- read_gmt(system.file('extdata', 'pathways.gmt',
#'                    package='CEMiTool'))
#' # Run overrepresentation analysis
#' cem <- mod_ora(cem, gmt)
#' # Plot module gene expression profiles
#' cem <- plot_ora(cem)
#' # Check resulting plot
#' show_plot(cem, "ora")
#'
#' @rdname plot_ora
#' @export
setGeneric('plot_ora', function(cem, ...) {
    standardGeneric('plot_ora')
})

#' @rdname plot_ora
setMethod('plot_ora', signature('CEMiTool'),
    function(cem, n=10, pv_cut=0.05, ...){
        if(length(unique(cem@module$modules)) == 0){
            stop("No modules in CEMiTool object! Did you run find_modules()?")
        }
        if(nrow(cem@ora) == 0){
            stop("No ORA data! Did you run mod_ora()?")
        }

        #cem <- get_args(cem=cem, vars=mget(ls()))
        ora_splitted <- split(cem@ora, cem@ora$Module)
        mod_cols <- mod_colors(cem)
        res <- lapply(ora_splitted, function(x){
                      plot_ora_single(head(x, n=n),
                      pv_cut=pv_cut,
                      graph_color=mod_cols[unique(x$Module)],
                      title=unique(x$Module),
                      ...)
        })
        modules <- names(res)
        modules <- modules[order(as.numeric(stringr::str_extract(modules, "\\d+")))]
        cem@barplot_ora <- res[modules]
        return(cem)
    })

#' ORA visualization for one module
#'
#' @keywords internal
#'
#' @param es a data.frame from ora function containing only one module
#' @param ordr_by column to order the data.frame
#' @param max_length max length of a gene set name
#' @param pv_cut p-value cuttoff
#' @param graph_color color of bars
#' @param title title of the graph
#'
#' @return a list with ggplot2 object and the number of significant gene sets
plot_ora_single <- function(es, ordr_by='p.adjust', max_length=50, pv_cut=0.05,
                            graph_color="#4169E1", title="Over Representation Analysis"){

    comsub <- function(x){
        #split the first and last element by character
        d_x <- strsplit(x[c(1, length(x))], "")
        #search for the first not common element, and so, get the last matching one
        der_com <- match(FALSE, do.call("==", d_x))-1
        return(substr(x, 1, der_com + 1))
    }

    es[, "GeneSet"] <- es[, "ID"]

    # limits name length
    ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length) # overflow
    ovf_data <- es[ovf_rows, "GeneSet"]
    test <- strtrim(ovf_data, max_length)
    dupes <- duplicated(test) | duplicated(test, fromLast=TRUE)
    if(sum(dupes) > 0){
        test[dupes] <- ovf_data[dupes]
        test[dupes] <- comsub(test[dupes])
        max_length <- max(nchar(test))
    }

    es[ovf_rows, "GeneSet"] <-  paste0(strtrim(test, max_length), "...")
    es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)

    # order bars
    lvls <- es[order(es[, ordr_by], decreasing=TRUE), "GeneSet"]
    es[, "GeneSet"] <- factor(es[, "GeneSet"], levels=lvls)

    es[, "alpha"] <- 1
    es[es[, ordr_by] > pv_cut, "alpha"] <- 0

    # Avoid 0's
    es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
    my_squish <- function(...){
        return(scales::squish(..., only.finite=FALSE))
    }

    # plot
    y_axis <- paste('-log10(', ordr_by, ')')
    pl <- ggplot(es, aes_string(x="GeneSet", y=y_axis, alpha="alpha", fill=y_axis)) +
        geom_bar(stat="identity") +
        theme(axis.text=element_text(size=8), legend.title=element_blank()) +
        coord_flip() +
        scale_alpha(range=c(0.4, 1), guide="none") +
        labs(y="-log10(adjusted p-value)", title=title, x="") +
        geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") +
        scale_fill_gradient(low="gray", high=graph_color, limits=c(2, 5), oob=my_squish)
    res <- list('pl'=pl, numsig=sum(es[, ordr_by] < pv_cut, na.rm=TRUE))
    return(res)
}

#' GSEA visualization
#'
#' Creates a heatmap with the results of gene set enrichment analysis (GSEA) of co-expression modules
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param pv_cut P-value cut-off. Default \code{0.05}
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with GSEA plots
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get example sample annotation file
#' # Run GSEA on network modules
#' cem <- mod_gsea(cem)
#' # Plot GSEA results
#' cem <- plot_gsea(cem)
#' # Check resulting plot
#' show_plot(cem, "gsea")
#'
#' @rdname plot_gsea
#' @export
setGeneric('plot_gsea', function(cem, ...) {
    standardGeneric('plot_gsea')
})

#' @rdname plot_gsea
setMethod('plot_gsea', signature('CEMiTool'),
    function(cem, pv_cut=0.05) {
        if(length(unique(cem@module$modules)) == 0){
            stop("No modules in CEMiTool object! Did you run find_modules()?")
        }
        if(length(cem@enrichment) == 0){
            stop("No GSEA data! Did you run mod_gsea()?")
        }
        if(all(unlist(lapply(cem@enrichment, nrow))) == 0){
            warning("No modules were enriched for any classes. Unable to plot enrichment.")
            return(cem)
        }
        #cem <- get_args(cem, vars=mget(ls()))

        stats <- names(cem@enrichment)
        enrichment <- lapply(cem@enrichment, function(stat){
                                 stat[is.na(stat)] <- 0
                                 rownames(stat) <- stat[,1]
                                 stat[,1] <- NULL
                                 return(stat)
        })
        names(enrichment) <- stats

        pval <- enrichment[['padj']]
        nes <- enrichment[['nes']]

        pval <- pval[rowSums(pval < pv_cut) >= 1, , drop=FALSE]
        nes <- nes[rownames(pval), , drop=FALSE]

        # check if there is any signif. module
        if(nrow(nes) < 0){
            stop("No significant modules found!")
        }

        custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                        "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                        "#D6604D", "#B2182B", "#67001F")
        custom_pal <- colorRampPalette(custom_pal)(200)

        nes <- as.matrix(nes)
        pval <- as.matrix(pval)
        nes[which(pval > pv_cut, arr.ind=TRUE)] <- 0

        if(nrow(nes) > 2){
            row_order <- rownames(nes)[hclust(dist(nes))$order]
        } else {
            row_order <- rownames(nes)
        }

        nes_melted <- reshape2::melt(nes)
        colnames(nes_melted) <- c("Module", "Class", "NES")
        nes_melted$Module <- factor(nes_melted$Module, levels=row_order)
        nes_melted$Class <- as.character(nes_melted$Class)
        max_abs_nes <- max(abs(nes_melted$NES))
        res <- ggplot(nes_melted, aes_(x=~Class, y=~Module, size=~abs(NES), fill=~NES)) +
            geom_point(color = "white", shape=21) +
            scale_fill_gradientn(colours=custom_pal, space = "Lab",
                                 limits=c(-max_abs_nes, max_abs_nes)) +
            scale_size(range=c(0,30), limits=c(0, NA)) +
            guides(size="none") +
            theme_minimal() +
            theme(panel.grid.major = element_blank()) +
            scale_x_discrete(position = "top")
        res_list <- list(enrichment_plot=res)
        cem@enrichment_plot <- res_list

        return(cem)
    })

#' Network visualization
#'
#' Creates a graph based on interactions provided
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param n number of nodes to label
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with profile plots
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get example gene interactions data
#' int <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
#' int_df <- read.delim(int)
#' # Include interaction data into CEMiTool object
#' interactions_data(cem) <- int_df
#' # Plot resulting networks
#' cem <- plot_interactions(cem)
#' # Check resulting plot
#' show_plot(cem, "interaction")
#'
#' @rdname plot_interactions
#' @export
setGeneric('plot_interactions', function(cem, ...) {
    standardGeneric('plot_interactions')
})

#' @rdname plot_interactions
setMethod('plot_interactions', signature('CEMiTool'),
    function(cem, n=10, ...) {
        if(length(unique(cem@module$modules)) == 0){
            stop("No modules in CEMiTool object! Did you run find_modules()?")
        }
        if(length(interactions_data(cem)) == 0){
            stop("No interactions information! Did you run interactions_data()?")
        }
        #cem <- get_args(cem, vars=mget(ls()))
        mod_cols <- mod_colors(cem)
        mod_names <- names(cem@interactions)
        mod_names <- mod_names[which(mod_names!="Not.Correlated")]
        hubs <- lapply(get_hubs(cem), names)
        zero_ints <- character()
        zero_ints <- lapply(names(cem@interactions), function(mod){
                        degree <- igraph::degree(cem@interactions[[mod]], normalized=FALSE)
                        if(length(degree) == 0) {
                            zero_ints <- append(zero_ints, mod)
                        }
                     })
        zero_ints <- unlist(zero_ints)
        if(!is.null(zero_ints)){
            mod_names <- mod_names[which(!(mod_names %in% zero_ints))]
        }
        if(length(mod_names) == 0){
            warning("There are no interactions in the given modules. Please check interactions file.")
            return(cem)
        }
        res <- lapply(mod_names, function(name){
                   plot_interaction(ig_obj=cem@interactions[[name]],
                                    n=n, color=mod_cols[name], name=name,
                                    mod_genes=module_genes(cem, name)$genes,
                                    coexp_hubs=hubs[[name]])
               })
        names(res) <- mod_names
        mod_names_ordered <- mod_names[order(as.numeric(stringr::str_extract(mod_names, "\\d+")))]
        cem@interaction_plot <- res[mod_names_ordered]
        return(cem)
    })

plot_interaction <- function(ig_obj, n, color, name, mod_genes, coexp_hubs){
    degrees <- igraph::degree(ig_obj, normalized=FALSE)
    ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
    max_n <- min(n, length(degrees))
    net_obj <- intergraph::asNetwork(ig_obj)
    m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
    # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
    plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
    # or get it them from Kamada-Kawai's algorithm:
    # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
    colnames(plotcord) <- c("X1","X2")
    edglist <- network::as.matrix.network.edgelist(net_obj)
    edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
    plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
    plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
    plotcord[, "shouldLabel"] <- FALSE
    plotcord[, "Hub"] <- ""
    int_hubs <- names(sort(degrees, decreasing=TRUE))[1:max_n]
    int_bool <- plotcord[, "vertex.names"] %in% int_hubs
    plotcord[which(int_bool), "Hub"] <- "Interaction"
    sel_vertex <- int_hubs
    if(!missing(coexp_hubs)){
        coexp_bool <- plotcord[, "vertex.names"] %in% coexp_hubs
        coexp_and_int <- coexp_bool & int_bool
        plotcord[which(coexp_bool), "Hub"] <- "Co-expression"
        plotcord[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
        sel_vertex <- c(sel_vertex, coexp_hubs)
    }

    colnames(edges) <-  c("X1","Y1","X2","Y2")
    #edges$midX  <- (edges$X1 + edges$X2) / 2
    #edges$midY  <- (edges$Y1 + edges$Y2) / 2
    plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
    plotcord$Degree_cut <- cut(plotcord$Degree, breaks=3, labels=FALSE)
    plotcord$in_mod <- TRUE
    #mod_genes <- cem@module[cem@module$modules==name,]$genes
    not_in <- setdiff(plotcord[,"vertex.names"], mod_genes)
    plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE

    pl <- ggplot(plotcord)  +
        geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
                     size = 0.5, alpha=0.5, colour="#DDDDDD") +
        geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
        geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),
                         box.padding=unit(1, "lines"),
                         data=function(x){x[x$shouldLabel, ]}) +
        scale_colour_manual(values=c("Co-expression" = "#005E87",
                                     "Interaction" = "#540814",
                                     "Co-expression + Interaction" = "#736E0B")) +
        labs(title=name) +
        ggplot2::theme_bw(base_size = 12, base_family = "") +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white",
                                                                colour = NA),
                       panel.border = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())

    return(pl)
}

#' Soft-threshold beta selection graph
#'
#' Creates a graph showing each possible soft-threshold value and its corresponding R squared value
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param plot_title title of the graph
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with beta x R squared plot
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot scale-free model fit as a function of the soft-thresholding beta parameter choice
#' cem <- plot_beta_r2(cem)
#' # Check resulting plot
#' show_plot(cem, "beta_r2")
#'
#' @rdname plot_beta_r2
#' @export
setGeneric('plot_beta_r2', function(cem, ...){
    standardGeneric('plot_beta_r2')
})

#' @rdname plot_beta_r2
setMethod('plot_beta_r2', signature('CEMiTool'),
    function(cem, plot_title="Scale independence (beta selection)"){
          if(nrow(cem@fit_indices) == 0){
            stop("No fit indices data! Did you run find_modules()?")
        }
        #cem <- get_args(cem, vars=mget(ls()))

        fit <- cem@fit_indices
        fit$new_fit <- -sign(fit[, 3])*fit[, 2]
        beta_power <- cem@parameters$beta

        pl <- ggplot(fit, aes_(x=~Power, y=~new_fit)) +
              geom_line(color="darkgrey") +
              geom_point(size=1.5) +
              theme(axis.text=element_text(size=12),
                    plot.title=element_text(hjust=0.5)) +
              scale_y_continuous(breaks=seq(round(min(fit$new_fit), 1), 1, by=0.2),
                                 limits=c(NA, 1)) +
              labs(y="Scale-free topology model fit, R squared", title=plot_title,
                   x="Soft-threshold beta")

        if(!is.na(beta_power)){
            pl <- pl + annotate(geom="text", label=beta_power,
                                x=beta_power,
                                y=fit[fit$Power == beta_power, "new_fit"] + 0.1,
                                color="red", size=7)
        }

        res_list <- list(beta_r2_plot=pl)
        cem@beta_r2_plot <- res_list
        return(cem)
    })

#' Network mean connectivity
#'
#' Creates a graph showing the mean connectivity of genes in the network
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param title title of the graph
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with connectivity plot
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot scale-free model fit as a function of the soft-thresholding beta parameter choice
#' cem <- plot_mean_k(cem)
#' # Check resulting plot
#' show_plot(cem, "mean_k")
#'
#' @rdname plot_mean_k
#' @export
setGeneric('plot_mean_k', function(cem, ...){
    standardGeneric('plot_mean_k')
})

#' @rdname plot_mean_k
setMethod('plot_mean_k', signature('CEMiTool'),
    function(cem, title="Mean connectivity"){
        if(nrow(cem@fit_indices) == 0){
              stop("No fit indices data! Did you run find_modules()?")
        }
        # cem <- get_args(cem, vars=mget(ls()))
        fit <- cem@fit_indices

        pl <- ggplot(fit, aes_(x=~Power, y=~mean.k.)) +
            geom_line(color="darkgrey") +
            geom_point(size=1.5) +
            theme(axis.text=element_text(size=12), plot.title=element_text(hjust=0.5)) +
            labs(y="Mean connectivity", title=title, x="Soft-threshold beta")

          res_list <- list(mean_k_plot=pl)
          cem@mean_k_plot <- res_list
        return(cem)
    })


#' Retrieve CEMiTool object plots
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param value A character string containing the name of the plot to be shown.
#' One of "profile", "gsea", "ora", "interaction", "beta_r2", "mean_k", "sample_tree",
#' "mean_var", "hist", "qq".
#'
#' @return A plot corresponding to a CEMiTool analysis
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot beta x R squared graph
#' cem <- plot_beta_r2(cem)
#' # Check plot
#' show_plot(cem, "beta_r2")
#' @rdname show_plot
#' @export
setGeneric('show_plot', function(cem, value) {
    standardGeneric('show_plot')
})

#' @rdname show_plot
setMethod('show_plot', signature('CEMiTool'),
    function(cem, value=c("profile", "gsea", "ora", "interaction",
                          "beta_r2", "mean_k", "sample_tree", "mean_var",
                          "hist", "qq")) {
        value <- match.arg(value)
        if(value!="sample_tree"){
            x_plot <- switch(value,
                             profile=cem@profile_plot,
                             gsea=cem@enrichment_plot,
                             ora=cem@barplot_ora,
                             interaction=cem@interaction_plot,
                             beta_r2=cem@beta_r2_plot,
                             mean_k=cem@mean_k_plot,
                             mean_var=cem@mean_var_plot,
                             hist=cem@hist_plot,
                             qq=cem@qq_plot)
            return(x_plot)
        }else{
            grid::grid.draw(cem@sample_tree_plot)
        }
    })

#' @title
#' Save CEMiTool object plots
#'
#' @description
#' Save plots into the directory specified by the \code{directory} argument.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param value A character string containing the name of the plot to be saved.
#' @param directory Directory into which the files will be saved.
#' @param force If the directory exists, execution will not stop.
#' @param ... Optional parameters
#' One of "all", "profile", "gsea", "ora", "interaction", "beta_r2", "mean_k",
#' "sample_tree", "mean_var", "hist", "qq".
#'
#' @return A pdf file or files with the desired plot(s)
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot beta x R squared graph
#' cem <- plot_beta_r2(cem)
#' # Save plot
#' \dontrun{save_plots(cem, value="beta_r2", directory="./Plots")}
#' @rdname save_plots
#' @export
setGeneric('save_plots', function(cem, ...) {
    standardGeneric('save_plots')
})

#' @rdname save_plots
setMethod('save_plots', signature('CEMiTool'),
    function(cem, value=c("all", "profile", "gsea", "ora",
                          "interaction", "beta_r2", "mean_k",
                          "sample_tree", "mean_var", "hist", "qq"),
                  force=FALSE, directory="./Plots") {
        if(dir.exists(directory)){
            if(!force){
                stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
            }
        }else{
            dir.create(directory, recursive=TRUE)
        }
        value <- match.arg(value)
        if(value == "all"){
            plots <- list(cem@profile_plot, cem@enrichment_plot, cem@barplot_ora,
                          cem@interaction_plot, cem@beta_r2_plot, cem@mean_k_plot,
                          cem@sample_tree_plot, cem@mean_var_plot, cem@hist_plot,
                          cem@qq_plot)
            all_plots<- c("profile", "gsea", "ora", "interaction", "beta_r2", "mean_k",
                                    "sample_tree", "mean_var", "hist", "qq")
            names(plots) <- all_plots

            plots <- Filter(function(x) length(x) >= 1, plots)
            if(length(plots) < length(all_plots)){
                message("Some plots have not been defined. Please run the appropriate plot functions. Saving available plots.")
            }
            lapply(names(plots), function(pl){
                pdf(file=file.path(directory, paste0(pl, ".pdf")))
                if(pl=="sample_tree"){
                    if(!is.null(nrow(plots[[pl]]))){
                        grid::grid.draw(plots[[pl]])
                        dev.off()
                    }
                }else{
                    print(plots[[pl]])
                    dev.off()
                }
            })
        }else if(value=="sample_tree"){
            if(!is.null(nrow(cem@sample_tree_plot))){
                pdf(file.path(directory, paste0(value, ".pdf")))
                grid::grid.draw(cem@sample_tree_plot)
                dev.off()
            }else{
                stop("Plot not available! Please use the appropriate plotting function on the CEMiTool object.")
            }
        }else{
            x_plot <- switch(value,
                             profile=cem@profile_plot,
                             gsea=cem@enrichment_plot,
                             ora=cem@barplot_ora,
                             interaction=cem@interaction_plot,
                             beta_r2=cem@beta_r2_plot,
                             mean_k=cem@mean_k_plot,
                             mean_var=cem@mean_var_plot,
                             hist=cem@hist_plot,
                             qq=cem@qq_plot)
            pdf(file.path(directory, paste0(value, ".pdf")))
            print(x_plot)
            dev.off()
        }
      })
