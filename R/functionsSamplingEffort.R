### Source file with functions
# Source it in R: source(functionsSamplingEffort.R)

# Function for selecting random rows/samples within a site ---------------

randomRows <- function(data, n) {
    ## data: dataframe with site (1st col), sample (2nd), and species abundances
    ## n: number of samples to sample within a site

    return(data[sample(nrow(data), n), ])
}



# Function for running DCA on reduced datasets -----------------------------

dcaRedData <- function(data, iter, n) {
    ## data:        data frame containing samples nested within sites/stations
    # First col: site/station (named Sta), second col: sample nr, then species
    ## iter:        number of random sampling iterations
    ## n            number of samples to randomly select per site/station
    
    require(plyr)
    require(reshape2)
    require(vegan)
    
    dcaList <- list()
    richnessTot <- NULL
    # richnessSite <- NULL
    richnessSite <- matrix(nrow = nrow(data),
                           ncol = iter)
    
    for(i in 1:iter)
    {
        # Random selection of a given number n of samples within each site
        sel.i <- ddply(data, .(Sta), randomRows, n)
        sel.i <- sel.i[, -2] # Removing sample id
        
        # Unifying data into abundance of species per site
        siteMelt.i <- melt(sel.i, "Sta")
        siteCast.i <- dcast(siteMelt.i, formula = Sta ~ variable, sum)
        site.i <- siteCast.i[, -1]
        
        # Weighting of abundance: R=32
        siteR32.i <- as.data.frame(site.i^(log(32)/log(max(site.i))))
        
        # Total and site richness
        richTot.i <- sum(colSums(siteR32.i) > 0)
        richSite.i <- as.vector(rowSums(siteR32.i > 0))
        
        dcaSites.i <- decorana(siteR32.i)
        
        axesSites.i <- data.frame(dca1 = scores(dcaSites.i, choices = 1, disp = "sites"),
                                  dca2 = scores(dcaSites.i, choices = 2, disp = "sites"))
        
        dcaList[[length(dcaList) + 1]] <- axesSites.i
        richnessTot[length(richnessTot) + 1] <- richTot.i
        # richnessSite[length(richnessSite) + 1] <- richSite.i
        richnessSite[, i] <- richSite.i
        
    }
    out <- list(dcas = dcaList,
                richnessTot = richnessTot,
                richnessSite = richnessSite)    
    out
}



# Function to extract DCA axes and store them in a matrix -----------------

extractAxes <- function(lst, ax, noSites, iter) {
    
    ## lst: list of data frames
    ## ax: axis to extract (1 or 2)
    ## noSites: number of sites in study
    ## iter: number of randomly selected reduced datasets
    
    axes <- matrix(data = unlist(lapply(lst, function(x) x[, ax])),
                   nrow = noSites,
                   ncol = iter,
                   dimnames = list(1:noSites, paste("rand", 1:iter, sep = "")))
}


# Function for checking and changing the sign of DCA axes ------------------------------

signCheck <- function(source, target) {
    ## source: vector of axes scores with known direction
    ## target: dataframe with columns that contain axes scores for which the sign are checked
    
    check <- cor(cbind(source, target))
    for(i in 1:ncol(target)) {
        check.i <- check[1, i]
        # print(sign(check.i))
        ifelse(sign(check.i == 1), next, target[, i] <- -target[, i])
    }
    target
}



# Function conducting procrustes analysis ---------------------------------

procrAnalysis <- function(source, first, second) {
    ## source: data frame/matrix containing ordination axes of the full data set
    ## first: data frame/matrix containing first axes of reduced data sets
    ## second: data frame/matrix containing second axes of reduced data sets
    # coeffs <- rep(NA, ncol(first))
    procrRes <- list()
    
    for(i in 1:ncol(first)) {
        target.i <- cbind(first[, i], second[, i])
        procr.i <- protest(X = source, Y = target.i,
                           scores = "sites",
                           symmetric = T,
                           permutations = 999)
        # coeffs[i] <- procr.i$scale
        procrRes[[length(procrRes) + 1]] <- procr.i
    }
    procrRes
}


# Function selecting procrustes coefficients and then plotting ------------

# procrPlots <- function(coeffs) {
#     ### Makes a 3x2 plot showing procrustes plots for the 0, 20, 40, 60, 80, and 100 
#     ### percentiles of procrustes coefficients
#     ## coeffs: vector of procrustes coefficients
#     ## 
#     index <- rep(NA, 6)
#     perc <- unname(round(quantile(coeffs, probs = seq(0, 1, 0.2)), 2))
#     for(i in seq_len(perc)) {
#         index[i] <- match(perc[i], round(coeffs, 2))
#     }
#     
#     pdf(file = paste("output/", "procr41Plots.pdf", sep = ""),
#         height = 10, width = 10)
#     
# #     svg(file = paste(avhandling, "fig2_procr41Plots.svg", sep = ""),
# #         height = 10, width = 10,
# #         onefile = TRUE)
#     
#     par(mfrow = c(3, 2))
#     
#     for(i in index)
#     {
#         plot(procr41Ls[[i]],
#              xlab = "DCA1", ylab = "DCA2", 
#              main = paste("r = ", round(procr41[i], 2), sep = ""),
#              cex.main = 1.5,
#              xlim = c(-0.3, 0.6), ylim = c(-0.25, 0.25))
#         text(procr41Ls[[i]],
#              display = "rotated",
#              labels = 1:28,
#              cex = 1.25,
#              pos = 2,
#              offset = 0.25)
#     }
#     dev.off()
# }
# 
