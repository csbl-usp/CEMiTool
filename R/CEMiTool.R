 
#!/usr/bin/env Rscript

base.path <- getwd()

" CEMiTool - Co-Expression Modules Identification Tool - Co-Expression Module analysis made easy

Just submit your data and CEMiTool will identify Co-Expression gene Modules and send results to your e-mail 

(Full article: 

Usage: CEMiTool.R [-d DIR] EXPRSFILE [-i INTERACT -o OUT -p CORPVAL -q ORAPVAL -n PERM -m MINGEN -c COR -t ANNOT -s GMT -e GENNUM --min-module=MINMOD --set-beta=BETA --diss-thres=THRESH --samples-column=SAMPLECOL --gene-column=GENECOL -f -b -u -l -v --merge-bool]

Input:
EXPRSFILE           a normalized expression file .txt format

Output:
A label_phi.pdf file (for all accounts, label is the name you gave the task and phi is the
resulting parameter, phi) with an R2 x curve, showing which value CEMiTool chose for the 
co-expression module analyses (the value cut by the red line). 
A label_Modules.pdf file with the expression profiles of the genes in each module, separated 
into the positively and negatively correlated gene submodules A and B, and a network image 
of your modules.
A label_Genes_in_Modules.txt file, with information pertaining the modules the genes are 
in.
A label_PPIs.txt file, with information about interactions between each gene pair connected
in the network.

Options:
-d DIR --directory=DIR                set working directory [default: ./]
-i INTERACT --interact=INTERACT       optional gene interaction file
-o OUT --output=OUT                   name output file [default: results]
-p CORPVAL --pvalue=CORPVAL           p-value for gene-gene correlation [default: 0.05]
-q ORAPVAL --pvalue-ora=ORAPVAL       p-value for over representation analysis [default: 0.05]
-n PERM --permutations=PERM           number of permutations [default: 1000]
-m MINGEN --mingenes=MINGEN           minimal number of genes per SUBMODULE [default: 20]
-c COR --correlation=COR              selected correlation method [default: spearman]
-t ANNOT --template=ANNOT             optional template file
-s GMT --gmtfile=GMT                  optional gmt file (geneSets)
-e GENNUM --genenum=GENNUM            number of genes to filter [default: 4000]
--min-module=MINMOD                   minimum MODULE size [default: 30]
--set-beta=BETA                       override selected beta value
--diss-thres=THRESH                   module merging correlation threshold for eigengene similarity [default: 0.8]
--gene-column=GENECOL	              the column name containing gene symbols in expression file [default: Symbol]
--samples-column=SAMPLECOL            the column name containing sample names in template file [default: SampleName]
-f --filter                           filter expression file by 60% most expressed and 'genenum' most variant genes [default: FALSE]
-b --beta-bool                        change WGCNA beta
-u --unsign                           use unsigned networks
-l --split                            split modules
-v --verbose                          display progress messages
--merge-bool                          merge related modules based on eigengene similarity

Authors:

Gustavo R Ferreira - gustavo.rodrigues.ferreira at usp.br
Pedro S T Russo - pedro.russo at usp.br
Matheus C Burger - burger at usp.br
Thiago D C Hirata - thiagodch at gmail.com
Helder I Nakaya - hnakaya at usp.br

More information:
www.csbiology.com
Computational Systems Biology Laboratory
University of Sao Paulo, Brazil

" -> doc

########################################################################
if (!interactive() && !exists('SOURCE')) {
    # Get and check arguments.
    suppressMessages(library(docopt))
    arg <- docopt(doc, version="1.0.2\n", strict=T)
    #arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))
    
    
    if (!file.exists(arg$exprsfile)) stop("No such file: ", arg$exprsfile)
    
    if(!is.null(arg$interact)){
        if (!file.exists(arg$interact)){
            stop("No such file: ", arg$interact)
        }
    }
    if (!grepl('\\.txt$', arg$exprsfile) && !grepl('\\.tsv$', arg$exprsfile)) stop("Bad format: ", arg$exprsfile)
    if(!is.null(arg$interact)){
        if (!grepl('\\.txt$', arg$interact)) stop("Bad format: ", arg$interact)
    }
    if (!grepl('^[[:digit:]]+$', arg$permutations)) stop("'n' must be an integer")
    if (!grepl('^[[:digit:]]+$', arg$mingenes)) stop("'m' must be an integer")
    
    if(!is.null(arg$template)){
        if (!file.exists(arg$template)){
            stop("No such file: ", arg$template)
        }
    }
    if(!is.null(arg$gmtfile)){
        if (!file.exists(arg$gmtfile)){
            stop("No such file: ", arg$gmtfile)
        }
    }
}   

set_dir     <- arg[["directory"]] #Ex: /Users/helder/Desktop/CEMtool/Test
exprs_f     <- arg[["exprsfile"]] #Ex: Expression_dataset2.txt
name_out    <- arg[["output"]] #Ex: Lepto
cutPvalue   <- as.numeric(arg[["pvalue"]]) #Cutoff permutation-based p-value for gene-gene correlation. Default = 0.05
nPerm       <- as.numeric(arg[["permutations"]]) #Number of permutations. Default = 1000.
MinSize     <- as.numeric(arg[["mingenes"]]) #Minimal number of genes in sub-modules. Default = 20
#Default = c(topo.colors(16),rainbow(16),heat.colors(20),
#            topo.colors(16),rainbow(16),heat.colors(20))
corr.method <- tolower(arg[["correlation"]])
interact <- arg[["interact"]]
template_f  <- arg[["template"]]
gmt_f       <- arg[["gmtfile"]]
filt        <- arg[["filter"]]
genenum     <- as.numeric(arg[["genenum"]])
beta_bool   <- arg[["beta_bool"]]
unsign      <- arg[["unsign"]]
split       <- arg[["split"]]
verbose     <- arg[["verbose"]]
gene_column <- arg[["gene_column"]]
samples_column <- arg[["samples_column"]]
minModuleSize <- as.numeric(arg[["min_module"]])
set_beta    <- arg[["set_beta"]]
merge_bool  <- arg[["merge_bool"]]
MEDissThres <- as.numeric(arg[["diss_thres"]])

# Transform set_beta into numeric if not NULL
if(!is.null(set_beta)) set_beta <- as.numeric(set_beta)
# Add paths to filename variables
if(!is.null(exprs_f)){
    exprs_f <- normalizePath(exprs_f)
}
if(!is.null(interact)){
    interact <- normalizePath(interact)
}
if(!is.null(template_f)){
    template_f <- normalizePath(template_f)
}
if(!is.null(gmt_f)){
    gmt_f <- normalizePath(gmt_f)
}
# Check variable values
print(c(set_dir, exprs_f, name_out, cutPvalue, nPerm, MinSize, corr.method, interact, template_f, gmt_f, filt, genenum, beta_bool, unsign, split, verbose))

suppressMessages({
    library("ggplot2")
    sink("/dev/null")
    library(WGCNA)
    sink()
    library("networkD3")
    library("reshape2") #for melt function
    library("matrixStats") #for rowCounts function
    library("igraph")
    require("cowplot") #for plot_grid function
    library("scales")
    library("org.Hs.eg.db")
    library("corrplot")
    library("mGSZ")
    library("dplyr")
    library("data.table")
    library("stringr") # for axis label wrapping (ORA) 
    library("gplots")  # for printing the parameters table
    library("rmarkdown")
    library("htmlwidgets")
    library("ff")
    library("foreach")
    library("doParallel")
    library("iterators")
    library("BiocParallel")
    if(!require("devtools")){
        install.packages("devtools")
        library("devtools")
    }
    if(!require("plotflow")){
        install_github("trinker/plotflow")
        library("plotflow")
    }
    if(!require("fgsea")){
        install_github("ctlab/fgsea")
        library("fgsea")
    }
})
allowWGCNAThreads()


#' Determines soft-threshold and creates co-expression modules.
#' 
#' @param data An expression data.frame.
#' @param file_name A string for the output file name.
#' @param corr.method A string for the correlation method 
#' @return List of WGCNA results.
goWGCNA <- function(data, file_name, corr.method = "spearman"){
    # This function takes the expression data,
    # determines the soft-threshold, merges based on eigengene similarity
    # and finds the co-expression modules
    D <- t(data)
    names(D) <- rownames(data)
    rownames(D) <- names(data)
    # Define a range of soft-thresholding candidates
    powers <- c(c(1:10), seq(12, 20, 2))
    
    # Selecting type of network to be used
    networkType <- ifelse(unsign == TRUE, "unsigned", "signed") # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
    
    #Use Spearman as default. Cite paper:
    #Evaluation of Gene Association Methods for 
    #Coexpression Network Construction and Biological Knowledge Discovery
    
    # Automatic selection of soft-thresholding power beta
    if(verbose) message("Selecting Beta") # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
    beta <- pickSoftThreshold(D, powerVector = powers, verbose = 5, 
                              corOptions = list(use = "p", method = corr.method), networkType = networkType, moreNetworkConcepts = TRUE)
    fit <- -sign(beta$fitIndices[, 3])*beta$fitIndices[, 2]
    k   <- beta$fitIndices[, 5]
    At  <- powers[length(powers)] - powers[1]
    A   <- 0.0
    for(cont in 2:length(fit)){
        A <- A + 0.5*(fit[cont] + fit[cont-1])*(powers[cont] - powers[cont-1])
    }
    #Area under the curve/threshold
    phi <- A/At 
    eps <- 0.1
    st  <- c(NA, NA)
    
    # Selecting smallest beta value in Cauchy sequence range (CEMiTool beta)
    for(count in (1:(length(fit)-2))){
        if(fit[count] >= 0.8){
            d <- c(abs(fit[count] - fit[count+1]), abs(fit[count] - fit[count+2]), abs(fit[count+1] - fit[count+2]))
            if(max(d) < eps){
                j <- which.max(k[count:count+2]) + count - 1
                st <- c(fit[j], powers[j])
                break
            }
        }
    }
    
    # Choosing between set beta, WGCNA beta and CEMiTool beta
    if(!is.null(set_beta)){ # THIS IS BAD AND YOU SHOULD FEEL BAD: referecen to global variable
        st[2] <- set_beta
        print("Using set beta:")
        print(st[2])
    }else if(!as.logical(beta_bool)){ # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
        st[2] <- beta$powerEstimate 
        print("Using WGCNA beta:")
        print(st[2])
    }else{
        print("Using CEMiTool beta:")
        print(st[2])
    }
    
    # Get network connectivities
    ourK <- NA
    if(!is.na(st[2])){
        print(st)
        if(st[2] <= 10){
            ourK <- beta$fitIndices[(st[2]) , 5]
        }else{
            line <- (st[2] + 10)/2
            ourK <- beta$fitIndices[line, 5]
        }
    }
    ourBeta <- as.integer(st[2])
    ourR2 <- st[1]
    
    # Calculating adjacency matrix
    D2 <- apply(D, 2, as.numeric)
    ourAdj <- adjacency(D2, power = as.numeric(ourBeta), type = networkType)
    
    #Calculating Topological Overlap Matrix
    ourTOM <- TOMsimilarity(ourAdj)
    # Determining TOM based distance measure 
    ourDiss <- 1 - ourTOM
    # Clustering
    ourTree <- hclust(as.dist(ourDiss), method = "average")
    # Cutting tree to determine modules
    ourMods <- cutreeDynamic(dendro = ourTree, distM = ourDiss, deepSplit = 2, 
                             pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
    ourTable <- table(ourMods)
    
    # Determining module colors
    ourColors <- labels2colors(ourMods)
    
    # Create pdf output for unmerged modules
    file_out <- paste(unlist(strsplit(file_name, ".txt")), "_", phi, "_", st[2], sep = "")
    pdf(file = paste0(file_out, ".pdf"),onefile=TRUE)
    ## Plot beta x R2
    plot(powers, fit, type = "n", xlab = "Soft Threshold Beta", 
         ylab = "signed R2", main = file_name, ylim = c(-0.1,1))
    text(powers, fit, labels = powers, col = "red")
    abline(h = st[1], col = "red")
    par(mfrow = c(2,1))
    ## Plot dendrogram
    plotDendroAndColors(ourTree, ourColors, "Dynamic Tree Cut", 
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                        main = "Gene dendrogram and module colors (recalculated beta)")
    dev.off()
    
    # Merging similar modules
    if(merge_bool){
        # Calculating eigengenes
        MEList <- moduleEigengenes(D, colors = ourColors)
        MEs <- MEList$eigengenes
        # Calculating dissimilarity of module eigengenes
        MEDiss <- 1 - cor(MEs)
        # Clustering module eigengenes
        METree <- hclust(as.dist(MEDiss), method = "average")
        # Plot cluster tree
        #sizeGrWindow(7, 6)
        #plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
        
        # Setting cut height
        MEDissThres <- 1 - MEDissThres
        # Merging modules
        moduleMerge <- mergeCloseModules(D, ourColors, cutHeight = MEDissThres, verbose = 3)
        # The merged module colors
        ourColors <- moduleMerge$colors
        # Eigengenes of the new merged modules
        mergedMEs <- moduleMerge$newMEs
        
        # Create pdf output for merged modules
        file_out <- paste(unlist(strsplit(file_name, ".txt")), "_", phi, "_", st[2], "_merged", sep = "")
        pdf(file = paste0(file_out, ".pdf"),onefile=TRUE)
        ## Plot beta x R2
        plot(powers, fit, type = "n", xlab = "Soft Threshold Beta", 
             ylab = "signed R2", main = file_name, ylim = c(-0.1,1))
        text(powers, fit, labels = powers, col = "red")
        abline(h = st[1], col = "red")
        par(mfrow = c(2,1))
        ## Plot dendrogram
        plotDendroAndColors(ourTree, ourColors, "Dynamic Tree Cut", 
                            dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                            main = "Gene dendrogram and module colors (recalculated beta)")
        dev.off()
    }
    
    # 
    # # Create pdf output 
    # file_out <- paste(unlist(strsplit(file_name, ".txt")), "_", phi, "_", st[2], sep = "")
    # pdf(file = paste0(file_out, ".pdf"),onefile=TRUE)
    # ## Plot beta x R2
    # plot(powers, fit, type = "n", xlab = "Soft Threshold Beta", 
    #      ylab = "signed R2", main = file_name, ylim = c(-0.1,1))
    # text(powers, fit, labels = powers, col = "red")
    # abline(h = st[1], col = "red")
    # par(mfrow = c(2,1))
    # ## Plot dendrogram
    # plotDendroAndColors(ourTree, ourColors, "Dynamic Tree Cut", 
    #                     dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
    #                     main = "Gene dendrogram and module colors (recalculated beta)")
    # dev.off()
    
    data2 <- cbind(data, ourColors, ourMods)
    params <- c(phi, ourBeta)
    names(params) <- c("phi", "ourBeta")
    res <- list(data=data2, parameters=params)
    
    params <- list(exprs_f, set_dir, name_out, cutPvalue, nPerm, MinSize, corr.method, interact, template_f, gmt_f, filt, genenum, ourBeta, phi, ourR2, ourK, ncol(exp.df))
    params <- unlist(lapply(params, function(x) ifelse(is.null(x), NA, x)))
    
    analysis.res <- data.frame(t(params))
    colnames(analysis.res) <- c("exprs_f", "set_dir", "name_out", "cutPvalue", "nPerm", "MinSize", "corr.method", "interact", "template_f", "gmt_f", "filt", "genenum", "ourBeta", "phi", "ourR2", "ourK", "colnumber")
    #write.table(analysis.res, file=paste0(name_out, "_analysis_res.txt"), sep="\t", quote=FALSE, row.names=FALSE)
    
    return(list(res, analysis.res))
}

SplitModules <- function(CorMatrix, Gpos, cutPvalue, AorB, mod, geneS, corr.method = "spearman", label) {
    # This function takes the expression data for a single module and refines it. The steps are:
    # - Separating into 2 submodules
    # - Evaluating edges through a permutation test
    # - Cutting edges that do not contribute to module conectedness
    
    Mpos <- CorMatrix[Gpos, Gpos]
    Mpos[lower.tri(Mpos, diag=TRUE)] <- NA
    
    posL <- melt(Mpos)
    posL[,1] <- as.character(posL[,1])
    posL[,2] <- as.character(posL[,2])
    posL <- posL[!is.na(posL$value), ]
    
    #posL[1:2] <- t( apply(posL[1:2], 1, sort))
    #posL[1:2]
    
    #posL      <- posL[!duplicated(posL[1:2]),]
    posL <- posL[order(posL[, 3], decreasing = TRUE), ]
    posL <- posL[which(posL[, 3] > 0), ]
    Gpos2 <- unique(c(posL[, 1], posL[, 2]))
    PPI_wanted_pos <- posL[, c(1,2)]
    
    
    if(split){ # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
        NewM_M <- ifelse(AorB == 2, paste0("M", mod, ".A"), paste0("M", mod, ".B"))
    }else{
        NewM_M <- paste0("M.", mod)
    }
    
    NewM_color <- ifelse(AorB == 2, all_colors[(mod*2)], all_colors[(mod*2)-AorB])
    
    out_wgcna[Gpos2,"NewModule"] <- NewM_M
    out_wgcna[Gpos2,"NewColor"]  <- NewM_color
    
    results_pos <- PPI_wanted_pos
    results_pos[, "Correlation"] <- CorMatrix[as.matrix(PPI_wanted_pos)]
    
    #no_cores <- detectCores() - 1
    #registerDoParallel(no_cores)
    
    
    tgeneSGpos2 <- t(geneS[as.character(Gpos2),])
    # for (i in 1:nPerm) {
    #     if(i %% 25 == 0) print(i)
    #     random_pos <- geneS[as.character(Gpos2), sample(ncol(geneS), ncol(geneS))]
    #     CorMatrixB_pos <- cor(t(random_pos), tgeneSGpos2, use = "everything", method = corr.method)
    #     results_pos <- cbind(results_pos, CorMatrixB_pos[as.matrix(PPI_wanted_pos)])
    # }
    # colnames(results_pos) <- append(colnames(results_pos)[1:3], 1:nPerm)
    
    if(verbose) message("Doing permutations") # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
    
    results_perm <- foreach(icount(nPerm)) %do% { # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
        random_pos <- geneS[as.character(Gpos2), sample(ncol(geneS), ncol(geneS))]
        CorMatrixB_pos <- cor(t(random_pos), tgeneSGpos2, use = "everything", method = corr.method)
        #results_pos <- cbind(results_pos, CorMatrixB_pos[as.matrix(PPI_wanted_pos)])
        CorMatrixB_pos[as.matrix(PPI_wanted_pos)]
    }
    results_perm <- as.data.frame(results_perm)
    results_pos <- cbind(results_pos, results_perm)
    colnames(results_pos) <- append(colnames(results_pos)[1:3], 1:nPerm) # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
    stopImplicitCluster()
    
    #times <- results_pos[,4:ncol(results_pos)] < results_pos[,3]
    
    #times2 <- ff(results_pos[,4:ncol(results_pos)] - results_pos[,3] < 0)
    
    times <- list()
    corrcol <- results_pos[,3]
    
    for(col in 4:ncol(results_pos)){
        #print(col)
        comp <- as.ff(results_pos[,col] - corrcol < 0)
        times[[as.character(col)]] <- as.data.frame(as.ffdf(as.matrix(comp)))
    }
    
    times.df <- as.data.frame(times)
    
    #t <- list()
    #for(cl in 1:ncol(m3)){
    #    print(cl)
    #    comparacao <- as.ff(m3[, cl] - x < 0)
    #    t[[as.character(cl)]] <- as.data.frame(as.ffdf(as.matrix(comparacao)))
    #}
    
    
    #save(file="saves.RData", list=ls())
    
    #message("Counting number of TRUEs...")
    times_bool <- rowCounts(as.matrix(times.df))
    #message("Finished counting")
    
    PPI_wanted_pos[, "Correlation"] <- CorMatrix[as.matrix(PPI_wanted_pos)]
    PPI_wanted_pos[, "P.value"] <- 1-times_bool/nPerm # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
    #	message("Function SplitModules:: rows with p-value above ", cutPvalue, " : ", sum(PPI_wanted_pos[, "P.value"] >= cutPvalue), " (", sum(PPI_wanted_pos[, "P.value"] >= cutPvalue)/nrow(PPI_wanted_pos), ")")
    #	message("Function SplitModules:: rows with p-value below ", cutPvalue, " : ", sum(PPI_wanted_pos[, "P.value"] < cutPvalue), " (", sum(PPI_wanted_pos[, "P.value"] < cutPvalue)/nrow(PPI_wanted_pos), ")")
    PPI_wanted_pos <- PPI_wanted_pos[which(PPI_wanted_pos[, "P.value"] < cutPvalue), ]
    PPI_wanted_pos <- PPI_wanted_pos[order(PPI_wanted_pos[, "Correlation"], decreasing = TRUE), ]
    Gpos3 <- unique(c(PPI_wanted_pos[, 1], PPI_wanted_pos[, 2])) 
    
    PPI_wanted_final <- GetPPIs(PPI_wanted_pos, Gpos3) 
    
    column.name <- 'module' #paste0("CEMiTool_", label)
    
    PPI_wanted_final[, column.name] <- NewM_M
    colnames(PPI_wanted_final)[1:2] <- c("Gene1", "Gene2")
    return(list(g = Gpos2, newC = NewM_color, newM = NewM_M, ppi = PPI_wanted_final))
}

GetPPIs <- function (PPI_wanted_pos, Gpos3) {
    a <- 1
    b <- nrow(PPI_wanted_pos)
    while(b > a+1){
        m <- floor((b - a)/2) + a
        g <- graph.edgelist(as.matrix(PPI_wanted_pos[1:m, 1:2]), directed = FALSE)
        ifelse(length(which(Gpos3 %in% V(g)$name == FALSE))==0 && clusters(g)$no == 1, b <- m, a <- m)
    }
    ppiX <- PPI_wanted_pos[1:m+1, ]
    return (ppiX)
}

PlotLines <- function(ExpX, rects) {
    # This function creates the expression profile line graphs
    notExp <- c((ncol(ExpX) - 4):ncol(ExpX))
    corX   <- paste("#", ExpX[,"NewColor"], "66", sep = "") 
    modX   <- ExpX[, "NewModule"]
    title_pos <- paste("Module:", modX[1]) 
    ExpX   <- ExpX[saida$g, -notExp]
    Gmean  <- apply(ExpX, 2, mean)
    time <- c(1:(ncol(ExpX)))
    
    mmm <- as.data.frame(cbind(time, Gmean)) # THIS IS BAD AND YOU SHOULD FEEL BAD: bad variable name
    mmm <- cbind(mmm,rep("mean", nrow(mmm)))
    colnames(mmm)[3] <- "group"
    xx <- t(ExpX) # THIS IS BAD AND YOU SHOULD FEEL BAD: bad variable name
    xx_long <- melt(xx,id = colnames(xx)) # THIS IS BAD AND YOU SHOULD FEEL BAD: bad variable name
    coltime <- rep(time, nrow(ExpX))
    xx_long <- cbind(xx_long, coltime) 
    colnames(xx_long)[2] <- "gene"
    Lineplot <- ggplot() + 
        geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.2) +
        geom_line(data=xx_long, aes(x=coltime, y=value, group = gene), colour = unique(corX)) +
        geom_line(data = mmm, aes(x = time, y = Gmean, group = group)) +
        ylab(label = "Intensity") + 
        xlab("Samples") +
        ggtitle(title_pos) +
        scale_x_discrete(limits = rownames(mmm)) + 
        theme(plot.title = element_text(lineheight  = .8, face = "bold", colour = "black", size = 15),
              axis.title.y = element_text(face  = "bold", colour = "black", size = 15),
              axis.title.x = element_text(face  = "bold", colour = "black", size = 15),
              axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 8),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 6),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 8),
              legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
              legend.position="bottom") +
        scale_fill_discrete(breaks=levels(rects[, "col"])) +
        coord_cartesian(xlim=c(1.5,15.5)) +
        coord_cartesian(xlim=c(mmm$time[1]+0.5, mmm$time[nrow(mmm)]-0.5 )) +
        guides(fill=guide_legend(ncol=2,byrow=TRUE)) #Determine column breaks for the legend
    return(Lineplot)
}

CreateGGplot2Graph <- function (dataX, geneX, sampleNames, plot.text = TRUE) {
    # This function creates the ggplot2 network graph
    dataX$compare <- apply(dataX, 1, function(x){paste(sort(c(x[1], x[2])), collapse="_")} )
    dataX$origin <- "Corr"    
    
    g <- graph.edgelist(as.matrix(dataX[, 1:2]), directed = FALSE)
    
    E(g)$weights <- 1- dataX[,"Correlation"]
    
    
    if(!is.null(interact)){ # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
        interactions.1  <- paste0(interact.df[, "Gene1"], "_", interact.df[, "Gene2"])
        interactions.1 <- interactions.1[which(interactions.1 %in% dataX$compare)]
        interactions.2 <- paste0(interact.df[, "Gene2"], "_", interact.df[, "Gene1"])
        interactions.2 <- interactions.2[which(interactions.2 %in% dataX$compare)]
        interactions <- c(interactions.1, interactions.2)
        
        dataX$origin[dataX$compare %in% interactions] <- "Interact"
        
        intersection.df <- subset(dataX, dataX$compare %in% interactions)
        
        
        
        #linhas que estão no arquivo de interação mas que não estão no dataX (mas são ambos genes do msm módulo)
        interactions.1  <- paste0(interact.df[, "Gene1"], "_", interact.df[, "Gene2"])
        not.interactions.1 <- interactions.1[which(!(interactions.1 %in% dataX$compare))]
        interactions.2 <- paste0(interact.df[, "Gene2"], "_", interact.df[, "Gene1"])
        not.interactions.2 <- interactions.2[!which(!(interactions.2 %in% dataX$compare))]
        
        not.interactions <- as.data.frame(c(not.interactions.1, not.interactions.2), stringsAsFactors=FALSE)
        names(not.interactions)[1] <- "compare"
        not.interactions$Gene1 <- lapply(not.interactions$compare, function(x) strsplit(x, "_")[[1]][1])
        not.interactions$Gene2 <- lapply(not.interactions$compare, function(x) strsplit(x, "_")[[1]][2])
        not.interactions$Correlation <- 1
        not.interactions$P.value <- 0
        not.interactions[, 6] <- unique(saida$ppi[, 5])
        names(not.interactions)[6] <- paste("CEMiTool", name_out, sep="_")
        not.interactions$origin <- "Interact_not_Corr"
        not.interactions <- not.interactions[, c(2, 3, 4, 5, 6, 1, 7)]
        not.interactions <- not.interactions[which(not.interactions$Gene1 %in% row.names(geneX) && not.interactions$Gene2 %in% row.names(geneX)) ,]
        colnames(not.interactions)[1] <- "X"
        colnames(not.interactions)[2] <- "Y"
        
        int.not.int <- rbind(intersection.df, not.interactions)
    }
    
    
    g <- mst(g)
    
    dataX <- get.data.frame(g, what="edges")
    
    layout.used <- layout.graphopt(g)
    colnames(layout.used) <- c("x", "y")  
    #layout.used[, "nodes"] <- rownames(layout.used) # for pca
    
    dataX$compare <- apply(dataX, 1, function(x){paste(sort(c(x[1], x[2])), collapse="_")} )
    dataX$origin <- "Corr"
    
    my.df <- as.data.frame(layout.used) # for layout (igraph)
    my.df[,"nodes"] <- unlist(vertex.attributes(g)) # for layout (igraph)
    my.df[,"group"] <- geneX[unlist(my.df$nodes), "NewModule"]
    my.df[,"CorGrupo"] <- paste("#",geneX[unlist(my.df$nodes), "NewColor"],sep = "")
    
    dataX$from.x <- my.df[, "x"][match(dataX[,1], my.df$nodes)] 
    dataX$from.y <- my.df[, "y"][match(dataX[,1], my.df$nodes)]
    dataX$to.x   <- my.df[, "x"][match(dataX[,2], my.df$nodes)] 
    dataX$to.y   <- my.df[, "y"][match(dataX[,2], my.df$nodes)]
    
    
    if(!is.null(interact)){ # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
        int.not.int$from.x <- my.df[, "x"][match(int.not.int[,1], my.df$nodes)] 
        int.not.int$from.y <- my.df[, "y"][match(int.not.int[,1], my.df$nodes)]
        int.not.int$to.x   <- my.df[, "x"][match(int.not.int[,2], my.df$nodes)] 
        int.not.int$to.y   <- my.df[, "y"][match(int.not.int[,2], my.df$nodes)]
    }
    
    Edge_width   <- (5+2*log(abs(dataX[, 3])) - min(5 + 2*log(abs(dataX[, 3]))))/diff(range(5 + 2*log(abs(dataX[, 3])))) + 0.1
    
    colnames(dataX)[1:2] <- c("X", "Y")
    
    if(!is.null(interact)){
        colnames(int.not.int)[1:2] <- c("X", "Y")
    }
    
    dataX$weights      <- ifelse(is.null(interact), Edge_width, dataX$weights) # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
    Numbers <- table(my.df$group)
    Names   <- rownames(Numbers)
    rects   <- data.frame()
    
    ########################################NETWORK D3###############################################
    
    #         geneList <- as.data.frame(sort(unique(c(dataX[, "X"], dataX[, "Y"]))))
    #         names(geneList) <- "names"
    #         geneList[,"numbers"] <- seq(0, nrow(geneList)-1)
    #         geneList[, "color"] <- unique(my.df[,"CorGrupo"])
    #         
    #         geneLinks <- dataX[,1:3]
    #         geneLinks <- merge(x=geneLinks, y=geneList, by.x="X", by.y="names")
    #         geneLinks <- merge(x=geneLinks, y=geneList, by.x="Y", by.y="names")
    #         geneLinks <- geneLinks[,c(4, 6, 2, 1, 3, 5)]
    #         names(geneLinks) <- c("source", "target", "X", "Y", "Correlation", "color")
    #         
    #         moduleColor <- unique(geneLinks[,"color"])
    #JScolor <- paste0('d3.scale.ordinal().range(["', moduleColor, '"])')
    
    #         my.graph <- forceNetwork(
    #             Links=geneLinks, 
    #             Nodes=geneList,
    #             Source="source", 
    #             Target="target",
    #             Value="Correlation", 
    #             NodeID="names",
    #             Group="color", 
    #             opacity=1,
    #             linkDistance=JS("function(d){return d.value*10}"),
    #             colourScale=JS(JScolor),
    #             linkWidth=0.5,
    #             linkColour="grey",
    #             zoom=T)
    #         
    #         this.module <- unique(geneX[unlist(my.df$nodes), "NewModule"])
    #         saveWidget(my.graph, paste0("Module_", this.module, "_graph.html"), selfcontained = TRUE)
    
    ########################################NETWORK D3###############################################
    
    for (grupo in 1:length(Names))
    {
        if (Numbers[grupo] != "0") {
            minX  <- min(my.df[which(my.df$group == Names[grupo]), "x"])-(abs(min(my.df[which(my.df$group == Names[grupo]), "x"])/10))
            maxX  <- max(my.df[which(my.df$group == Names[grupo]), "x"])+(abs(max(my.df[which(my.df$group == Names[grupo]), "x"])/10))
            minY  <- min(my.df[which(my.df$group == Names[grupo]), "y"])-(abs(min(my.df[which(my.df$group == Names[grupo]), "y"])/10))
            maxY  <- max(my.df[which(my.df$group == Names[grupo]), "y"])+(abs(max(my.df[which(my.df$group == Names[grupo]), "y"])/10))
            corG  <- unique(my.df[which(my.df$group == Names[grupo]),5])
            row   <- data.frame(Names[grupo], minX, maxX, minY, maxY, corG)
            rects <- rbind(rects,row)
        }
    }
    c <- geom_rect(data = rects, aes(xmin = minX, xmax = maxX, ymin = minY, ymax = maxY), 
                   fill = "blue", alpha = 0, linetype = "dashed", colour = corG)
    d <- annotate(	"text",   x = (rects[,"maxX"]+rects[, "minX"])/2,
                   y = rects[,"maxY"]+(abs(rects[, "maxY"]/10)),
                   label = paste(rects[,1]," (", length(V(g))," nodes;", length(E(g)), " edges)", sep=""))
    
    
    mod.col <- unique(my.df$CorGrupo)
    
    graphplot <- ggplot() +
        geom_segment(data = dataX, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, 
                                       width = 0.5), color = "black", alpha = 0.5)
    
    
    if(!is.null(interact)){ # THIS IS BAD AND YOU SHOULD FEEL BAD: reference to global variable
        if(nrow(int.not.int) > 0){
            graphplot <- graphplot + geom_segment(data = int.not.int, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, 
                                                                          width = 0.5), color = mod.col, alpha = 0.5)
        }
    }
    
    graphplot <- graphplot +
        geom_point(data = my.df,aes(x = x, y = y), size = 8, colour = "white") + 
        geom_point(data = my.df,aes(x = x, y = y), size = 6, colour = my.df$CorGrupo, alpha = 0.5)
    
    if(plot.text){
        graphplot <- graphplot + geom_text(data = my.df, aes(x = x, y = y, label = nodes), size = 1.5, colour = "black") 
    }
    
    graphplot <- graphplot + 
        scale_x_continuous(expand = c(0,3)) +
        scale_y_continuous(expand = c(0,(abs(rects[,"maxY"]/9)))) + 
        theme_bw() +  c + d +
        theme(axis.text.x      = element_blank(), 
              axis.text.y      = element_blank(), 
              axis.ticks       = element_blank(),  
              axis.title.x     = element_blank(), 
              axis.title.y     = element_blank(), 
              panel.background = element_blank(), 
              panel.border     = element_blank(), 
              panel.grid.major = element_blank(),  
              panel.grid.minor = element_blank(),  
              legend.position  = "none", 
              plot.background  = element_blank())
    
    return(graphplot)
    
}

RankE <- function(D){
    ME <- moduleEigengenes(D, rep("blue", dim(D)[2]))
    L <- signedKME(D, ME$eigengenes)
    colnames(L) <- NULL
    return(L)
}

drop.col <- function(x, cols){
    return(x[, -which(colnames(x) %in% cols)])
}

filter.prop <- function(x, prop, fun=mean){
    rows <- floor(prop * nrow(x))
    return(filter.rows(x, rows, fun))
}

filter.rows <- function(x, rows, fun=sd){
    rows <- min(rows, nrow(x))
    val <- apply(x, 1, fun)
    sel.rows <- order(val, decreasing=TRUE)[1:rows]
    return(x[sel.rows, ])
}

just.filter <- function(x, genenum){
    if(floor(nrow(x)*.7) > genenum){
        y <- filter.prop(x, 0.7, fun=mean)
    }else{
        y <- x
    }
    res <- filter.rows(y, genenum, fun=sd)
    return(res)
}

createRects2 <- function(exp.df, template.df=NULL){
    if(!is.null(template.df)){
        rects2  <- data.frame(xstart = seq(0.5,ncol(exp.df)-0.5,1), 
                              xend = seq(1.5,ncol(exp.df)+0.5,1), col = factor(template.df[,"Class"], levels=unique(template.df[, "Class"])))
        classNames   <- rownames(table(template.df[,"Class"]))
        Zexp.df <- t(scale(t(exp.df),center=TRUE,scale=TRUE))
        Zmeans  <- Zexp.df
        Zmeans  <- data.frame(matrix(ncol = 1, nrow = nrow(Zexp.df)))
        rownames(Zmeans) <- rownames(Zexp.df)
        for (groupX in 1:length(classNames)){
            meanTemp <- apply(Zexp.df[, which(template.df[,"Class"]==classNames[groupX])], 1, mean)
            Zmeans <- cbind(Zmeans,meanTemp)
            colnames(Zmeans)[ncol(Zmeans)] <- classNames[groupX]
        }
        Zmeans <- Zmeans[,-1]
        return(list(Zmeans, rects2))
    }else{
        rects2  <- data.frame(xstart = seq(0.5,ncol(exp.df)-0.5,1), 
                              xend = seq(1.5,ncol(exp.df)+0.5,1), col = "no_class")
    }
    return(list(rects2))
}

########################################################################

all_colors <- rep(rainbow(16, s = 1, v = 0.7 ), 20)

#Remove "#"
all_colors <- substr(all_colors, 2, 7)
plot.text <- TRUE


setwd(set_dir)

if(verbose) message("Reading expression file")


exp.df <- read.delim(exprs_f, sep = "\t", header = TRUE, row.names = gene_column)
exp.gsea <- exp.df
#Remove ".CEL" from sample names
# sampleNames <- toupper(names(exp.df))
# sampleNames <- gsub("[_,.,;,-,,|,/,\\].*", "", sampleNames) 
# names(exp.df) <- sampleNames

################################

#Filter 70% most expressed genes and <genenum> most variant
if (filt == TRUE){
    if(verbose) message("Filtering")
    exp.df <- just.filter(exp.df, genenum)
}


if (!is.null(template_f)){
    if(verbose) message("Reading template file")
    
    template.df  <- as.data.frame(read.delim(template_f, sep = "\t", header = TRUE, row.names = samples_column))
    
    if(!any(names(template.df) == "Class")){
        stop("Template file does not contain a grouping column named 'Class'")
    }
    
    if (length(names(exp.df)) > length(row.names(template.df))){
        stop("Expression file contains more samples than template file")
    }else if (length(names(exp.df)) < length(row.names(template.df))){
        template.df <- template.df[colnames(exp.df),]
        warning("Template file contains more samples than expression file; cutting template file")
        #exp.df  <- exp.df[,rownames(template.df)]
    }
    
    if(!all((sort(names(exp.df)) == sort(rownames(template.df))))){
        rownames(template.df) <- make.names(rownames(template.df))
        if(!all((sort(names(exp.df)) == sort(rownames(template.df))))){
            stop("Template sample names different from expression file sample names")
        }else{
            message("Sample names started with numbers. Added syntactically valid names.")
        }
    }
    
    z.and.rects2 <- createRects2(exp.df, template.df)
    Zmeans <- z.and.rects2[[1]]
    rects2 <- z.and.rects2[[2]]
    
}else{
    ###Caso contrario (nao tenha template):
    # rects2  <- data.frame(xstart = seq(0.5,ncol(exp.df)-0.5,1), 
    #                       xend = seq(1.5,ncol(exp.df)+0.5,1), col = "no_class") 
    z.and.rects2 <- createRects2(exp.df)
    rects2 <- z.and.rects2[[1]]
}

write.gmt <- function(all.gmts, fname){
    # This function takes a list of modules and their genes and converts it into a gmt file
    fname <- paste0(fname, ".gmt")
    out <- file(fname, "w")
    alllines <- sapply(names(all.gmts), function(x){
        paste(c(x, "",all.gmts[[x]]), collapse="\t")
        
    })
    writeLines(alllines, con=out)
    close(out)
}

doFastGSEA <- function(exp.gsea, template.df, GS, ranks = F){
    # This function takes the expression and template files and the genes/module list and runs
    # fgsea to find the enrichment of each module's genes across the classes
    Temp <- template.df
    Temp$Class <- as.character(Temp$Class)
    classes <- as.vector(unique(Temp[, "Class"]))
    
    # Transforms into Z-scores
    Zexp.gsea <- data.frame(t(scale(t(exp.gsea), center=TRUE, scale=TRUE)))
    
    gseaList <- list()
    
    for(j in 1:length(classes)){
        print(classes[j])
        curr_class <- classes[j]
        class_samples <- rownames(subset(Temp, subset=Temp$Class==curr_class))
        print(j)
        print(class_samples)
        if(ranks){
            geneList <- rank(apply(exp.gsea[, class_samples], 1, mean))
            geneList <- sort(geneList, decreasing = T)
        }else{
            geneList <- apply(Zexp.gsea[, class_samples], 1, mean)
            geneList <- sort(geneList, decreasing = T)
        }
        register(SerialParam())
        bpparameters <- bpparam()
        message("bpparam()")
        message(str(bpparameters))
        if(verbose) message("Running fgsea")
        fgseaRes <- fgsea(pathways = GS, 
                          stats = geneList,
                          minSize=15,
                          maxSize=500,
                          nperm=10000,
                          nproc=0)
        lead.edge <- fgseaRes[["leadingEdge"]]
        lead.edge <- lapply(lead.edge, function(x){ 
            x <- paste(x, collapse=",")
        })
        lead.edge <- unlist(lead.edge)
        fgseaRes[["lead.edge"]] <- lead.edge
        fgseaRes[["leadingEdge"]] <- NULL
        gseaList[[classes[j]]] <- setDF(fgseaRes)
    }
    
    list.es <- lapply(gseaList, "[", c("pathway", "ES"))
    list.es <- lapply(names(list.es), function(x) setNames(list.es[[x]], paste0(x, "_", names(list.es[[x]]))))
    es.combined <- Reduce(function(x, y) merge(x, y, all=T, by.x=grep("pathway", colnames(x)), by.y=grep("pathway", colnames(y))), list.es, accumulate=F)
    names(es.combined)[1] <- "pathway"
    
    write.table(es.combined, file=paste0(name_out, "_Enrichment.txt"), sep="\t", row.names = F)
    
    list.pval <- lapply(gseaList, "[", c("pathway", "padj"))
    list.pval <- lapply(names(list.pval), function(x) setNames(list.pval[[x]], paste0(x, "_", names(list.pval[[x]]))))
    pval.combined <- Reduce(function(x, y) merge(x, y, all=T, by.x=grep("pathway", colnames(x)), by.y=grep("pathway", colnames(y))), list.pval, accumulate=F)
    names(pval.combined)[1] <- "pathway"
    
    write.table(pval.combined, file=paste0(name_out, "_Enrichment_Pvalue.txt"), sep="\t", row.names = F)
    
    list.nes <- lapply(gseaList, "[", c("pathway", "NES"))
    list.nes <- lapply(names(list.nes), function(x) setNames(list.nes[[x]], paste0(x, "_", names(list.nes[[x]]))))
    nes.combined <- Reduce(function(x, y) merge(x, y, all=T, by.x=grep("pathway", colnames(x)), by.y=grep("pathway", colnames(y))), list.nes, accumulate=F)
    names(nes.combined)[1] <- "pathway"
    
    write.table(nes.combined, file=paste0(name_out, "_Enrichment_NES.txt"), sep="\t", row.names = F)
    gsea.res <- list(es.combined, pval.combined, nes.combined)
    names(gsea.res) <- c("ES", "padj", "NES")
    return(gsea.res)
}    

doCorrplot <- function(gsea.res, pv.cut=0.05){
    
    res.es <- gsea.res[[1]]
    row.names(res.es) <- res.es$pathway
    res.es$pathway <- NULL
    
    res.pval <- gsea.res[[2]]
    row.names(res.pval) <- res.pval$pathway
    res.pval$pathway <- NULL
    
    res.nes <- gsea.res[[3]]
    row.names(res.nes) <- res.nes$pathway
    res.nes$pathway <- NULL
    
    res.pval[is.na(res.pval)] <- 0
    res.es[is.na(res.es)] <- 0
    res.nes[is.na(res.nes)] <- 0
    
    res.pval <- res.pval[rowSums(res.pval < pv.cut) >= 1, ]
    #RES.ES <- RES.ES[rownames(RES.pval),]
    res.nes <- res.nes[rownames(res.pval),]
    
    if(nrow(res.nes) < 0){
        warning("No significant modules found!")
    }
    
    pdf(file=paste0(name_out, "_Enrichment.pdf"), width=max(2, ncol(res.es)/nrow(res.es)), max(2, height=nrow(res.es)/ncol(res.es)))
    corrColors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                    "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
    corrColors <- colorRampPalette(corrColors)(200)
    #es.mat <- as.matrix(RES.ES)
    nes.mat <- as.matrix(res.nes)
    pv.mat <- as.matrix(res.pval)
    nes.mat[which(pv.mat > pv.cut, arr.ind=T)] <- 0
    
    order <- row.names(nes.mat)[hclust(dist(nes.mat))$order]
    
    corrplot(nes.mat[order, ], p.mat=pv.mat[order, ], col=corrColors, is.corr=FALSE, addgrid.col="white", insig="blank",
             pch.cex=0.5, pch.col="black", tl.col="black", tl.cex=0.5, cl.cex=0.4, cl.ratio=0.5,
             cl.pos="r", cl.align.text="l", mar=c(0,0,0,0), sig.level=pv.cut)
    dev.off()
}    

get.parameters <- function(pathway.genes, significant.genes, all.genes){
    params <- list()
    significant.genes <- unique(significant.genes)
    pathway.genes <- unique(pathway.genes)
    params[['intersection']] <- intersect(significant.genes, pathway.genes)
    params[['white.balls.drawn']] <- length(params[['intersection']])
    params[['white.balls.in.urn']] <- length(pathway.genes)
    params[['total.balls.in.urn']] <- all.genes
    params[['black.balls.in.urn']] <- params[['total.balls.in.urn']] - params[['white.balls.in.urn']]
    params[['balls.pulled.from.urn']] <- length(significant.genes)
    params[['black.balls.drawn']] <- params[['balls.pulled.from.urn']] - params[['white.balls.drawn']]
    params[['white.balls.left']] <- params[['white.balls.in.urn']] - params[['white.balls.drawn']]
    params[['black.balls.left']] <- params[['black.balls.in.urn']] - params[['black.balls.drawn']]
    params[['confusion.matrix']] <- matrix(c(params[['white.balls.drawn']], params[['white.balls.left']],
                                             params[['black.balls.drawn']], params[['black.balls.left']]),
                                           ncol=2, nrow=2, byrow=T)
    return(params)
} 

#'Do hypergeometric test
#' adapted from https://stat.ethz.ch/pipermail/bioconductor/2013-September/054904.html
#'
#' @return p-value
hyperg.test <- function(...){
    params <- list(...)
    p.val <- phyper(q=params[['white.balls.drawn']]-1, m=params[['white.balls.in.urn']], 
                    n=params[['black.balls.in.urn']], k=params[['balls.pulled.from.urn']], lower.tail=FALSE)
    return(p.val)
}

fisher.exact.test <- function(...){
    params <- list(...)
    res.test <- fisher.test(params[['confusion.matrix']], alternative="greater")
    res <- c(res.test$p.value, res.test$estimate)
    names(res) <- c('p.value', 'odds.ratio')
    return(res)
}

do.enrichment <- function(pathway.genes, significant.genes, all.genes, method=hyperg.test){
    params <- get.parameters(pathway.genes, significant.genes, all.genes)
    return(do.call(method, params))
}

names.in.warning <- function(x, max.names=10){
    names.out <- x
    if(length(names.out) > max.names){
        names.out <- c(names.out[1:10], "...")
    }
    names.p <- paste(names.out, collapse=", ")
    return(names.p)
}

symbol2eid <- function(genes, env=org.Hs.egSYMBOL2EG, remove.NAs=TRUE){
    ids <- mget(genes, env, ifnotfound=NA)
    lens <- sapply(ids, length)
    multiple.ids <- names(lens)[which(lens > 1)]
    if(length(multiple.ids) > 1){
        names.p <- names.in.warning(multiple.ids)
        warning(length(multiple.ids), "Some genes have more than one entrez gene id, considering only first.\nGenes: ", names.p)
    }
    without.ids <- names(ids)[sapply(ids, function(x) all(is.na(x)))]
    res <- sapply(ids, "[", 1)
    if(length(without.ids) > 1){
        names.p <- names.in.warning(without.ids)
        warning(length(without.ids), " genes does not have an entrez gene id.\nGenes: ", names.p)
        
        if(remove.NAs){
            res <- res[-which(is.na(res))]
        }
    }
    return(res)
}

apply.ora <- function(id, genesets, significant.genes, universe){
    col.names <- c('id', 'name', 'hyper', 'hyper.adj', 'fisher', 'fisher.adj', 'odds.ratio', 
                   'gene.set.length', 'query.length', 'intersection')
    row.list <- list('id'=id)
    genes.in.set <- genesets[[id]]
    params <- get.parameters(genes.in.set, significant.genes, universe)
    row.list[['intersection']] <- paste(params[['intersection']], collapse=", ")
    row.list[c('fisher', 'odds.ratio')] <- do.call(fisher.exact.test, params)
    row.list[['hyper']] <- do.call(hyperg.test, params)
    row.list[['hyper.adj']] <- NA
    row.list[['fisher.adj']] <- NA
    row.list[['query.length']] <- length(significant.genes)
    row.list[['gene.set.length']] <- length(genes.in.set)
    #row.list[['name']] <- as.character(get(str_match(id, "\\d+"), KEGGPATHID2NAME))
    row.list[['name']] <- NA
    row.list <- as.data.frame(row.list[col.names])
    return(row.list)
}

read.gmt <- function(fname){
    res <- list(genes=list(), desc=list())
    gmt <- file(fname)
    gmt.lines <- readLines(gmt)
    close(gmt)
    gmt.list <- lapply(gmt.lines, function(x) unlist(strsplit(x, split="\t")))
    gmt.names <- sapply(gmt.list, '[', 1)
    gmt.desc <- lapply(gmt.list, '[', 2)
    gmt.genes <- lapply(gmt.list, function(x){x[3:length(x)]})
    names(gmt.desc) <- names(gmt.genes) <- gmt.names
    res[['genes']] <- gmt.genes
    res[['desc']] <- gmt.desc
    return(res)
}

do.gmt <- function(genes, genesets, universe, adjust.method="fdr"){
    if(missing(universe)){
        universe <- length(org.Hs.egSYMBOL)
        message("Parameter universe missing. Using human genes as background. #genes:", universe)
    }
    res.df <- do.call(rbind, lapply(names(genesets), apply.ora, genesets=genesets, 
                                    significant.genes=genes, universe=universe))
    res.df[, 'hyper.adj'] <- p.adjust(res.df[, 'hyper'], method=adjust.method)
    res.df[, 'fisher.adj'] <- p.adjust(res.df[, 'hyper'], method=adjust.method)
    res.df <- res.df[order(res.df[, 'fisher.adj']), ]
    rownames(res.df) <- NULL
    return(res.df)
}

my.squish <- function(...){
    return(squish(..., only.finite=FALSE))
}

#ordr.by <- 'hyper.adj'
#es <- res.reactome[res.reactome[, ordr.by] < 0.01, ]
#ggsave(plot=pl, filename="teste.svg", height=height, width=width)
plot.ora <- function(es, ordr.by='hyper.adj', maxlength=50, pv.cut=0.01, graphColor, title){
    # paste id and name
    ids <- as.character(es[, "id"])
    ids[is.na(ids)] <- ""
    if(!all(is.na(es[, "name"]))){
        name <- as.character(es[, "name"])
        name[is.na(name)] <- ""
        es[, "GeneSet"] <- paste(ids, name, sep="_")
    }else{
        es[, "GeneSet"] <- ids
    }
    
    # limits name length
    es[which(nchar(es[, "GeneSet"])>maxlength), "GeneSet"] <- paste0(strtrim(es[which(nchar(es[, "GeneSet"])>maxlength), "GeneSet"], maxlength), "...")
    es[, "GeneSet"] <- str_wrap(es[,"GeneSet"], width = 20)
    
    # order bars
    lvls <- es[order(es[, ordr.by], decreasing=TRUE), "GeneSet"]
    es[, "GeneSet"] <- factor(es[, "GeneSet"], levels=lvls)
    
    es[, "alpha"] <- 1
    es[es[, ordr.by] > pv.cut, "alpha"] <- 0
    
    # Avoid 0's
    es[es[, ordr.by] > 0.8, ordr.by] <- 0.8
    
    # plot
    pl <- ggplot(es, aes_string(x="GeneSet", y=paste('-log10(', ordr.by, ')'), alpha="alpha", fill=paste('-log10(', ordr.by, ')'))) + 
        geom_bar(stat="identity") +
        #scale_y_continuous(limits=c(0,5), oob=my.squish) + 
        #scale_x_discrete(labels = str_wrap(es[,"GeneSet"], width = 20)) + 
        theme(axis.text=element_text(size=8), legend.title=element_blank()) +
        coord_flip() + scale_alpha(range=c(0.4, 1), guide="none") +
        labs(y="-log10(adjusted p-value)", title=title, x="") +
        geom_hline(yintercept=-log10(pv.cut), colour="grey", linetype="longdash") + 
        scale_fill_gradient(low="gray", high=paste0("#",graphColor), limits=c(2, 5), oob=my.squish)
    height <- .12*nrow(es)
    width <- .08*max(nchar(es[, 'name'])) + 5
    res <- list('pl'=pl, 'height'=height, 'width'=width, numsig=sum(es[, ordr.by] < pv.cut, na.rm=TRUE))
    return(res)
}

gset <- NA
if (!is.null(gmt_f)){
    
    gmt <- read.gmt(gmt_f)
    genesets <- gmt[['genes']]
    names.universe <- row.names(exp.df)
    universe <- length(names.universe)
    
    # filter genes not in microarray
    genesets.filt <- lapply(genesets, function(x){x[x %in% names.universe]})
    
    if(all(sapply(genesets.filt, length)==0)){
        warning("There is no intersection between GMT and genes on expression file.")
    }
    
    # filter geneset with less than 5 genes
    gset <- genesets.filt[which(sapply(genesets.filt, length)>5)]
}

sample.names <- colnames(exp.df)

# Run WGCNA
if(verbose) message("WGCNA")
wgcna_data_params <- goWGCNA(exp.df, name_out, corr.method = corr.method)
wgcna_results  <- wgcna_data_params[[1]][["data"]]
out_wgcna <- data.frame(wgcna_results)
out_wgcna[, "NewColor"] <- "NA"
out_wgcna[, "NewModule"] <- "NA"
all_ppi <- data.frame()
all_enrich <- data.frame() # enrichment results
modNumbers <- table(wgcna_results[, ncol(wgcna_results)] )
modNames   <- rownames(modNumbers)
if (modNames[1] == "0"){
    if(modNumbers['0'] == nrow(exp.df)){
        writeLines("Error! We could not specify the parameter Beta. No modules found", con = paste0(name_out, "_error.txt"), sep = "\n", useBytes = FALSE)
        warning("We could not specify the parameter Beta. No modules found! Exiting ...")
        quit()
    }
    modNumbers <- as.table(modNumbers[-1 , drop=FALSE])
    modNames   <- rownames(modNumbers)
}

finalBeta <- wgcna_data_params[[1]][['parameters']]['ourBeta']
finalPhi <- round(wgcna_data_params[[1]][['parameters']]['phi'], digits = 3)
parameter.names <- c("Directory", "Exprs.File", "Output.Name", "Cutoff.Pvalue", "No.Permutations", "Min.Genes", "Corr.Method", "Interaction.File", "Template.File", "Gmt.File", "Filter", "Num.Filt", "Beta", "Phi") 
inputs <- c("set_dir", "exprs_f", "name_out", "cutPvalue", "nPerm", "MinSize", "corr.method", "interact", "template_f", "gmt_f", "filt", "genenum", "finalBeta", "finalPhi")

rm.rows <- which(sapply(inputs, exists) & sapply(inputs, function(x) is.null(get(x))))
if(length(rm.rows) > 0){
    inputs <- inputs[-rm.rows]
    parameter.names <- parameter.names[-rm.rows]
}

parameters <- sapply(inputs, get)

# Create modules pdf
modules.pdf <- paste0(name_out, "_Modules_tmp.pdf")

if(verbose) message("Creating modules")
pdf(file = modules.pdf) #, onefile = TRUE)

results.table <- data.frame("Names" = parameter.names,  "Values" = parameters, row.names=NULL)

if (!is.null(interact)){
    # Read interaction file
    interact.df  <- as.data.frame(read.table(interact, sep = "\t", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA"))
    interact.df[, "module"] <- NA 
    names(interact.df)[1:2] <- c("Gene1", "Gene2")
    interact.df$origin <- "Interact"
}

eigen.list <- data.frame()
mean.list <- data.frame()

for (mod in 1:length(modNames)){
    if(verbose) message(paste("Creating module", mod, "of", length(modNames)))
    
    geneS <- wgcna_results[which(wgcna_results[, ncol(wgcna_results)] == modNames[mod], arr.ind = TRUE), sample.names]
    C <- cor(t(geneS), use = "everything", method = corr.method)
    signs <- sign(range(C))
    if((signs[1] != signs[2]) && split==TRUE ){	# somente se houver correlacoes pos e neg
        DM <- as.dist(1 - C)
        H <- hclust(DM)
        k <- cutree(H, 2)
        Gpos <- rownames(C)[which(k == 1)] 
        Gneg <- rownames(C)[which(k == 2)]
        
        colorSeq <- c(paste0("M", mod, ".A"), paste0("M", mod, ".B"))
        dynamicColors <- labels2colors(k, colorSeq=colorSeq)
        
    }else{
        Gpos <- rownames(C) 
        Gneg <- c()
        
        if(split==TRUE){
            colorSeq <- paste0("M", mod, ".A")
        }else{
            colorSeq <- paste0("M.", mod)
        }
        
        dynamicColors <- rep(colorSeq, length(rownames(C)))
    }
    MEs <- moduleEigengenes(t(geneS), colors=dynamicColors)$eigengenes
    eigen.list <- rbind(eigen.list, t(MEs))
    
    if(verbose) message("Splitting module")
    if (length(Gpos) >= MinSize) {
        AorB <- 2
        saida <- SplitModules(C, Gpos, cutPvalue, AorB, mod, geneS = geneS, corr.method = corr.method, label = name_out)
        out_wgcna[saida$g,"NewColor"]  <- saida$newC
        out_wgcna[saida$g,"NewModule"] <- saida$newM
        out_wgcna[saida$g, "Membership"] <- rank(RankE(t(exp.df[saida$g, ])))
        all_ppi <- rbind(all_ppi,saida$ppi)
        ExpX   <- out_wgcna[saida$g, ]
        
        only.exp <- ExpX[, colnames(ExpX) %in% names(exp.df)]
        mean.exp <- t(apply(only.exp, 2, mean))
        
        if(split==TRUE){
            rownames(mean.exp) <- paste0("M", mod, ".A")
        }else{
            rownames(mean.exp) <- paste0("M.", mod)
        }
        mean.list <- rbind(mean.list, mean.exp)
        
        if(verbose) message("Plotting expression profile")
        plot1 <- PlotLines(ExpX,rects2)
        if(verbose) message("Plotting network graph")
        plot2 <- CreateGGplot2Graph(saida$ppi, ExpX, sampleNames, plot.text = plot.text)
        print(plot_grid(plot1, plot2, labels = c("a", "b"), ncol = 1, nrow = 2))
    } 
    if (length(Gneg) >= MinSize) {
        AorB <- 1
        saida <- SplitModules(C, Gneg, cutPvalue, AorB, mod, geneS = geneS, corr.method = corr.method, label = name_out)
        out_wgcna[saida$g,"NewColor"]  <- saida$newC
        out_wgcna[saida$g,"NewModule"] <- saida$newM
        out_wgcna[saida$g, "Membership"] <- rank(RankE(t(exp.df[saida$g, ])))
        all_ppi <- rbind(all_ppi,saida$ppi)
        ExpX   <- out_wgcna[saida$g, ]
        
        only.exp <- ExpX[, colnames(ExpX) %in% names(exp.df)]
        mean.exp <- t(apply(only.exp, 2, mean))
        
        if(split==TRUE){
            rownames(mean.exp) <- paste0("M", mod, ".B")
        }
        mean.list <- rbind(mean.list, mean.exp)
        
        if(verbose) message("Plotting expression profile")
        plot1 <- PlotLines(ExpX,rects2)
        if(verbose) message("Plotting network graph")
        plot2 <- CreateGGplot2Graph(saida$ppi, ExpX, sampleNames, plot.text = plot.text)
        print(plot_grid(plot1, plot2, labels = c("a", "b"), ncol = 1, nrow = 2))
    }
}

names(eigen.list) <- names(exp.df)
rownames(eigen.list) <- unlist(lapply(rownames(eigen.list), function(x){strsplit(x, "E")[[1]][2]}))
eigen.list <- eigen.list[rownames(mean.list),]

write.table(eigen.list, file = paste(name_out, "_eigenlist.txt", sep = ""), sep = "\t", row.names = T, quote = FALSE)
write.table(mean.list, file = paste(name_out, "_meanlist.txt", sep = ""), sep = "\t", row.names = F, quote = FALSE)
#plot3 <- CreateGGplot2GraphFULL(all_ppi, out_wgcna)
#print(plot3)

dev.off()

if (!is.null(interact)){
    if(verbose) message("Creating interactions file")
    #interact.df  <- as.data.frame(read.table(interact, sep = "\t", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA"))
    #interact.df[, "Module"] <- NA 
    same.mod <- which(out_wgcna[interact.df[, 1], "NewModule"] == out_wgcna[interact.df[, 2], "NewModule"])
    interact.df[same.mod, "module"] <- out_wgcna[interact.df[same.mod, 1], "NewModule"]
    interact.df <- interact.df[!(is.na(interact.df[,"module"]) | interact.df[,"module"]=="NA"),]
    all_ppi <- interact.df
}

write.table(all_ppi, file = paste(name_out, "_PPIs.txt", sep = ""), sep = "\t", row.names = F, quote = FALSE)


# Rename "GenesInModules" table column names and write it
new.names <- c(paste0("WGCNA_", name_out), 'module',#paste0("CEMiTool_", name_out), 
               paste0("CEMiTool_", name_out, "_membership"), paste0("WGCNA_", name_out, "_colors"))
names(new.names) <- c("ourMods", "NewModule", "Membership", "ourColors")
mtch <- match(colnames(out_wgcna), names(new.names))
colnames(out_wgcna)[colnames(out_wgcna) %in% names(new.names)] <- na.exclude(new.names[mtch])
out_wgcna <- drop.col(out_wgcna, sample.names)
write.table(out_wgcna, file = paste(name_out, "_GenesInModules.txt", sep = ""), sep="\t", col.names = NA, quote=FALSE)

###Transform each submodule into gene sets
subMs   <- rownames(table(out_wgcna[,'module'])) #paste0("CEMiTool_", name_out)])) 

GS <- list()

for (mX in 1:length(subMs)){
    Mgenes <- rownames(out_wgcna[which(out_wgcna[,'module']==subMs[mX]),]) 
    GS[[subMs[mX]]] <- as.character(Mgenes)
}

# Create gmt file with gene modules 
write.gmt(GS, name_out)

#################################fastGSEA#########################################
if(!is.null(template_f)){
    gsea.res <- doFastGSEA(exp.gsea, template.df, GS, ranks = F)  
    doCorrplot(gsea.res)
}
#################################fastGSEA#########################################

####################################mGSZ###########################################
#Vamos usar o Zmeans como ranks para o mGSZ ou o weighted Kolmogorov Smirnov
##http://ekhidna.biocenter.helsinki.fi/downloads/pashupati/mGSZ.html
##weighted Kolmogorov Smirnov
##Os submodulos serao os gene.sets que serao testados contra o Zmeans
##saida eh um corrplot com os gene.sets nas linhas e os ranks nas colunas

# if(!is.null(template_f)){
#     library(mGSZ)
#     doGSZ(GS)
# }

####################################mGSZ###########################################

# Enrichment Analysis (ORA)
sig_count_total <- NA
if(!is.null(gmt_f) & length(gset) > 0){
    if(verbose) message("Doing ORA")
    
    gmt.genes <- character()
    for(i in seq(1:length(gset))){
        gmt.genes <- as.character(c(gmt.genes, gset[[i]]))
    }
    gmt.genes <- unique(gmt.genes)
    
    if(sum(rownames(exp.df) %in% gmt.genes) > 0){
        
        pdf(paste(name_out, "_ORA.pdf", sep = ""))
        mod2genes <- split(rownames(out_wgcna), out_wgcna[, 'module'])#paste0("CEMiTool_", name_out)])
        
        if("NA" %in% names(mod2genes)){
            mod2genes <- mod2genes[-which(names(mod2genes) == "NA")]
        }
        mod2genes <- mod2genes[order(as.numeric(str_extract(names(mod2genes), "\\d+")))]
        
        res.gmt <- list()
        res.gmt$sig_count <- 0
        
        res.gmt <- lapply(names(mod2genes), function(name){
            x <- mod2genes[[name]]
            res.gmt <- as.data.frame(c(res.gmt, do.gmt(x, gset, universe, adjust.method="fdr")))
            res.gmt <- res.gmt[order(res.gmt[, "fisher"]), ]
            graphColor <- as.character(out_wgcna[match(name, out_wgcna[, 'module']), 'NewColor']) #paste0("CEMiTool_", name_out)]), "NewColor"])
            plot.res <- plot.ora(res.gmt[1:10, ], graphColor=graphColor, title=name)
            if(plot.res[["numsig"]] > 0){
                print(plot.res[['pl']])
                res.gmt$sig_count <- res.gmt$sig_count + 1
            }
            return(cbind(NewModule=name, res.gmt))
        })
        
        sig_count_total <- lapply(res.gmt, function(x){
            unique(x$sig_count)
        })
        sig_count_total <- sum(unlist(sig_count_total))
        
        names(res.gmt) <- names(mod2genes)
        dev.off()
        
        res.gmt <- do.call(rbind, res.gmt)
        write.table(res.gmt, file = paste(name_out, "_ORA.txt", sep = ""), sep = "\t", row.names = F, quote = FALSE, col.names=TRUE)
        
        #Create module x pathway table
        new.res.gmt <- res.gmt[, c("NewModule", "id", "hyper.adj")]
        new.res.gmt <- dcast(res.gmt, NewModule ~ id, value.var = "hyper.adj")  #Define NewModule as identifier, id as columns and hyper.adj as values
        write.table(new.res.gmt, file = paste(name_out, "_ORA_PathsInModules.txt", sep = ""), sep = "\t", row.names = F, quote = FALSE, col.names=TRUE)
    }else{
        message("Expression file gene symbols not in gmt file. Skipping ORA.")
        sig_count_total <- NA
    }
}


# Create analysis.res parameter table
analysis.res <- wgcna_data_params[[2]]
col <- 'module'#paste("CEMiTool_", name_out, sep="")
analysis.res$ourNGzero <- sum(out_wgcna[, col] == "NA")
analysis.res$ourNMods  <- length(unique(out_wgcna[,col][which(out_wgcna[, col] != "NA")]))
analysis.res$ourGM <- mean(table(out_wgcna[,col][which(out_wgcna[, col] != "NA")]))
analysis.res$sig_count <- sig_count_total

if(verbose) message("Outputting parameters")
write.table(analysis.res, file = paste(name_out, "_analysis_res.txt", sep = ""), sep="\t", col.names = NA, quote=FALSE)

# parameters
#parameters <- data.frame(name = c("Label: ", "p-value cutoff: ", "number of permutations: ", 
#				   		 			"minimum module size: ", "correlation method: "), value = args[-c(1,2)] )

#write.table(parameters, file = paste(name_out, "_params.txt", sep = ""), sep = "\t", row.names = F, quote = FALSE)

render(file.path(base.path, "template_modules.Rmd"), envir = globalenv(), output_file = paste0(name_out, "_params.pdf"), output_dir = getwd())

#save(file="saves.RData", list=ls())

gowgcna_plot_fname <- paste(name_out, "_", wgcna_data_params[[1]][["parameters"]]["phi"], "_", wgcna_data_params[[1]][["parameters"]]["ourBeta"], ".pdf", sep = "")
system(paste0("pdftk A=", name_out, "_params.pdf"," B=", modules.pdf, " C=", gowgcna_plot_fname, " cat A1 C B2-end output ", name_out, "_Modules.pdf"))

if(verbose) message("Finished CEMiTool")
#load(file="saves.RData") 
