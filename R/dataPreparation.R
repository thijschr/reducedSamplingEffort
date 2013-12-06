### Data preparation script


# Abundance weighting for samples within-sites data -----------------------

### Weight abundance to grabs (R=16, gives an R = 40 for sites, similar to the intermediate 
## range in GA paper)
## Data set with all grabs (samples)
allStaRep16 <- allStaRep
allStaRep16[, -c(1:2)] <- as.data.frame(allStaRep16[, -c(1:2)]^(log(16)/log(max(allStaRep16[, -c(1:2)]))))
allStaRep16$Sta <- as.factor(allStaRep16$Sta)
allStaRep16$Gr <- as.factor(allStaRep16$Gr)



# Abundance weighting for abundance-by-site data --------------------------

allSta32 <- allSta
allSta32[, -1] <- as.data.frame(allSta32[, -1]^(log(32)/log(max(allSta32[, -1]))))
allSta32$Sta <- as.factor(allSta32$Sta)
