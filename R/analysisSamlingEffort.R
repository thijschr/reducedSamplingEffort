### Analysis of sampling effort


# Loading data ------------------------------------------------------------

allStaRep <- read.csv("data/abundanceByStationAndSample.csv")
allSta <- read.csv("data/abundanceByStation.csv")



# Sourcing script with functions ------------------------------------------

source("R/functionsSamplingEffort.R")



# Sourcing data preparation script ----------------------------------------

source("R/dataPreparation.R")



# ANALYSIS OF INTRA-SITE VARIABILITY --------------------------------------

# ORDINATION OF FULL DATASET BY GRAB --------------------------------------

# decorana() call and extraction of axes scores ---------------------------
require(vegan)

dca4StaRep <- decorana(allStaRep16[, -c(1:2)])

axes4StaRep <- data.frame(dca1 = scores(dca4StaRep, choices = 1, display = "sites"),
                          dca2 = scores(dca4StaRep, choices = 2, display = "sites"))

## Adding Site/Station as grouping factor
axes4StaRep$Sta <- allStaRep16$Sta


# Extracting mean and range of DCA site-grab scores -----------------------

axes4StaRepMeanRange <- data.frame(
    d1 = aggregate(x = axes4StaRep$dca1,
                   by = list(axes4StaRep$Sta),
                   FUN = mean)[, 2],
    d1Lo = aggregate(x = axes4StaRep$dca1,
                     by = list(axes4StaRep$Sta),
                     FUN = range)[, 2][, 1],
    d1Hi = aggregate(x = axes4StaRep$dca1,
                     by = list(axes4StaRep$Sta),
                     FUN = range)[, 2][, 2],
    d2 = aggregate(x = axes4StaRep$dca2,
                   by = list(axes4StaRep$Sta),
                   FUN = mean)[, 2],
    d2Lo = aggregate(x = axes4StaRep$dca2,
                     by = list(axes4StaRep$Sta),
                     FUN = range)[, 2][, 1],
    d2Hi = aggregate(x = axes4StaRep$dca2,
                     by = list(axes4StaRep$Sta),
                     FUN = range)[, 2][, 2])

axes4StaRepMeanRange$deltaD1 <- with(axes4StaRepMeanRange, d1Hi - d1Lo)
axes4StaRepMeanRange$deltaD2 <- with(axes4StaRepMeanRange, d2Hi - d2Lo)



# PLOTTING INTRA-SITE VARIABILITY -----------------------------------------

require(extrafont)
loadfonts(device = "postscript")

postscript(file = paste("output/", "intraSiteVarPlot.eps", sep = ""),
           height = 6.25, width = 6.83,
           family = "Arial", paper = "special",
           onefile = FALSE, horizontal = FALSE)

# pdf(file = paste("output/", "intraSiteVarPlot.pdf", sep = ""),
#     height = 5.5, width = 6)

# svg(file = paste(avhandling, "fig1_intraSiteVarPlot.svg", sep = ""),
#     height = 5.5, width = 6,
#     onefile = TRUE)

# par(mfrow = c(1, 2))

## Complete DCA diagram
plot(x = axes4StaRepMeanRange$d1,
     y = axes4StaRepMeanRange$d2,
     xlab = "DCA1",
     ylab = "DCA2",
     xlim = c(-1.00, 1.75),
     ylim = c(-0.75, 1.25),
     main = "Intra-site variability",
     type = "n")

# DCA1 variability
segments(x0 = axes4StaRepMeanRange$d1Lo,
         y0 = axes4StaRepMeanRange$d2,
         x1 = axes4StaRepMeanRange$d1Hi,
         y1 = axes4StaRepMeanRange$d2)

# DCA2 variability
segments(x0 = axes4StaRepMeanRange$d1,
         y0 = axes4StaRepMeanRange$d2Lo,
         x1 = axes4StaRepMeanRange$d1,
         y1 = axes4StaRepMeanRange$d2Hi)

text(x = axes4StaRepMeanRange$d1,
     y = axes4StaRepMeanRange$d2,
     labels = 1:28) 
dev.off()



# WILCOXON-MANN-WHITNEY: test of DCA1 vs. DCA2 variability -----------------

with(axes4StaRepMeanRange, wilcox.test(deltaD1, deltaD2,
                                       paired = T))


# ORDINATION OF FULL DATASET BY SITE --------------------------------------

dcaSta <- decorana(allSta32[, -1])

axesSta <- data.frame(dca1 = scores(dcaSta, choices = 1, display = "sites"),
                      dca2 = scores(dcaSta, choices = 2, display = "sites"))


# DCA STABILITY FOR FULL VS REDUCED DATA SETS -----------------------------

# Running DCA, extracting axes and richness data for reduced datasets -----

dcas3gr <- dcaRedData(data = allStaRep, iter = 99, n = 3)
dcas2gr <- dcaRedData(data = allStaRep, iter = 99, n = 2)
dcas1gr <- dcaRedData(data = allStaRep, iter = 99, n = 1)


# Storing DCA axes for reduced datasets in data frames --------------------
dca1_3gr <- as.data.frame(extractAxes(lst = dcas3gr[[1]], ax = 1,
                                      noSites = 28, iter = 99))
dca2_3gr <- as.data.frame(extractAxes(lst = dcas3gr[[1]], ax = 2,
                                      noSites = 28, iter = 99))
dca1_2gr <- as.data.frame(extractAxes(lst = dcas2gr[[1]], ax = 1,
                                      noSites = 28, iter = 99))
dca2_2gr <- as.data.frame(extractAxes(lst = dcas2gr[[1]], ax = 2,
                                      noSites = 28, iter = 99))
dca1_1gr <- as.data.frame(extractAxes(lst = dcas1gr[[1]], ax = 1,
                                      noSites = 28, iter = 99))
dca2_1gr <- as.data.frame(extractAxes(lst = dcas1gr[[1]], ax = 2,
                                      noSites = 28, iter = 99))



# Checking and changing (if necessary) the sign of axes -------------------

## All axes need to point in the same direction when running procrustes analysis

## Changing sign (if necessary) by calling signCheck()
dca1_3gr <- signCheck(source = axesSta[, 1], target = dca1_3gr)
dca2_3gr <- signCheck(source = axesSta[, 2], target = dca2_3gr)
dca1_2gr <- signCheck(source = axesSta[, 1], target = dca1_2gr)
dca2_2gr <- signCheck(source = axesSta[, 2], target = dca2_2gr)
dca1_1gr <- signCheck(source = axesSta[, 1], target = dca1_1gr)
dca2_1gr <- signCheck(source = axesSta[, 2], target = dca2_1gr)



# Procrustes analysis -----------------------------------------------------

## Calling procrAnalysis()

procr3 <- procrAnalysis(target = axesSta, first = dca1_3gr, second = dca2_3gr)
procr2 <- procrAnalysis(target = axesSta, first = dca1_2gr, second = dca2_2gr)
procr1 <- procrAnalysis(target = axesSta, first = dca1_1gr, second = dca2_1gr)



# Correlation between corresponding DCA axes ------------------------------

## Vectors to hold correlation coefficients for each of the nine iterations
corr43_ax1 <- rep(NA, 99); corr43_ax2 <- rep(NA, 99)
corr42_ax1 <- rep(NA, 99); corr42_ax2 <- rep(NA, 99)
corr41_ax1 <- rep(NA, 99); corr41_ax2 <- rep(NA, 99)

for(i in seq(length(corr43_ax1)))
{
    # 4 VS 3
    corr43_ax1.i <- cor.test(axesSta[, 1], dca1_3gr[, i], method = "k")
    corr43_ax2.i <- cor.test(axesSta[, 2], dca2_3gr[, i], method = "k")
    
    ## 4 VS 2
    corr42_ax1.i <- cor.test(axesSta[, 1], dca1_2gr[, i], method = "k")
    corr42_ax2.i <- cor.test(axesSta[, 2], dca2_2gr[, i], method = "k")
    
    ## 4 VS 1
    corr41_ax1.i <- cor.test(axesSta[, 1], dca1_1gr[, i], method = "k")
    corr41_ax2.i <- cor.test(axesSta[, 2], dca2_1gr[, i], method = "k")
    
    # Store results
    corr43_ax1[i] <- abs(corr43_ax1.i$estimate); corr43_ax2[i] <- abs(corr43_ax2.i$estimate)
    corr42_ax1[i] <- abs(corr42_ax1.i$estimate); corr42_ax2[i] <- abs(corr42_ax2.i$estimate)
    corr41_ax1[i] <- abs(corr41_ax1.i$estimate); corr41_ax2[i] <- abs(corr41_ax2.i$estimate)
}



# Congruence and total species richness results ------------------------------------------------------

congruence <- data.frame(samples = c("3", "2", "1"),
                         procrMean = c(mean(unlist(lapply(procr3, function(x) x$scale))),
                                       mean(unlist(lapply(procr2, function(x) x$scale))),
                                       mean(unlist(lapply(procr1, function(x) x$scale)))),
                         procrSD = c(sd(unlist(lapply(procr3, function(x) x$scale))),
                                     sd(unlist(lapply(procr2, function(x) x$scale))),
                                     sd(unlist(lapply(procr1, function(x) x$scale)))),
                         procrMin = c(min(unlist(lapply(procr3, function(x) x$scale))),
                                      min(unlist(lapply(procr2, function(x) x$scale))),
                                      min(unlist(lapply(procr1, function(x) x$scale)))),
                         procrMax = c(max(unlist(lapply(procr3, function(x) x$scale))),
                                      max(unlist(lapply(procr2, function(x) x$scale))),
                                      max(unlist(lapply(procr1, function(x) x$scale)))),
                         tau1Mean = c(mean(corr43_ax1), 
                                      mean(corr42_ax1), 
                                      mean(corr41_ax1)),
                         tau1SD = c(sd(corr43_ax1), 
                                    sd(corr42_ax1), 
                                    sd(corr41_ax1)),
                         tau1Min = c(min(corr43_ax1), 
                                     min(corr42_ax1), 
                                     min(corr41_ax1)),
                         tau1Max = c(max(corr43_ax1), 
                                     max(corr42_ax1), 
                                     max(corr41_ax1)),
                         tau2Mean = c(mean(corr43_ax2), 
                                      mean(corr42_ax2), 
                                      mean(corr41_ax2)),
                         tau2SD = c(sd(corr43_ax2), 
                                    sd(corr42_ax2), 
                                    sd(corr41_ax2)),
                         tau2Min = c(min(corr43_ax2), 
                                     min(corr42_ax2), 
                                     min(corr41_ax2)),
                         tau2Max = c(max(corr43_ax2), 
                                     max(corr42_ax2), 
                                     max(corr41_ax2)),
                         spMean = c(mean(dcas3gr$richnessTot), 
                                    mean(dcas2gr$richnessTot), 
                                    mean(dcas1gr$richnessTot)),
                         spSD = c(sd(dcas3gr$richnessTot), 
                                  sd(dcas2gr$richnessTot), 
                                  sd(dcas1gr$richnessTot)),
                         spMin = c(min(dcas3gr$richnessTot), 
                                   min(dcas2gr$richnessTot), 
                                   min(dcas1gr$richnessTot)),
                         spMax = c(max(dcas3gr$richnessTot), 
                                   max(dcas2gr$richnessTot), 
                                   max(dcas1gr$richnessTot)))

                         

write.csv(congruence, paste("output/", "congruence.csv", sep = ""), row.names = F)



# Changes in species richness per site ------------------------------------

speciesRichSite <- data.frame(site = 1:28, 
                              sp4 = as.vector(rowSums(allSta[, -1] > 0)),
                              sp3 = apply(dcas3gr$richnessSite, 1, mean),
                              sp3sd = apply(dcas3gr$richnessSite, 1, sd),
                              sp3lo = apply(dcas3gr$richnessSite, 1, range)[1, ],
                              sp3hi = apply(dcas3gr$richnessSite, 1, range)[2, ],
                              sp2 = apply(dcas2gr$richnessSite, 1, mean),
                              sp2sd = apply(dcas2gr$richnessSite, 1, sd),
                              sp2lo = apply(dcas2gr$richnessSite, 1, range)[1, ],
                              sp2hi = apply(dcas2gr$richnessSite, 1, range)[2, ],
                              sp1 = apply(dcas1gr$richnessSite, 1, mean),
                              sp1sd = apply(dcas1gr$richnessSite, 1, sd),
                              sp1lo = apply(dcas1gr$richnessSite, 1, range)[1, ],
                              sp1hi = apply(dcas1gr$richnessSite, 1, range)[2, ])

write.csv(speciesRichSite, paste("output/", "speciesRichSite.csv", sep = ""), row.names = F)



# PROCRUSTES PLOTS --------------------------------------------------------


# Selection of iterations to plot, based on percentiles -------------------

indices <- percIndices(procr1)



# Plotting Procrustes plots to file ---------------------------------------

procrPlots(procr1, indices)



# Embedding fonts in figures ----------------------------------------------

embed_fonts(file = "output/intraSiteVarPlot.eps", 
            outfile = "output/intraSiteVarPlot_embed.eps",
            options = "-dEPSCrop")

embed_fonts(file = "output/procr41PLots.eps", 
           outfile = "output/procr41PLots_embed.eps",
           options = "-dEPSCrop")
