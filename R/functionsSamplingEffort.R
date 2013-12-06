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
    richnessSite <- NULL
    
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
        richnessSite[length(richnessSite) + 1] <- richSite.i
        
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
    coeffs <- rep(NA, ncol(first))
    
    for(i in 1:ncol(first)) {
        target.i <- cbind(first[, i], second[, i])
        procr.i <- protest(X = source, Y = target.i,
                           scores = "sites",
                           symmetric = T,
                           permutations = 999)
        coeffs[i] <- procr.i$scale
    }
    coeffs
}