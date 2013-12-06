### Analysis of sampling effort


# Loading data ------------------------------------------------------------

allStaRep <- read.csv("data/abundanceByStationAndSample.csv")
allSta <- read.csv("data/abundanceByStation.csv")



# Sourcing script with functions ------------------------------------------

source("R/functionsSamplingEffort.R")


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

pdf(file = paste("output/", "intraSiteVarPlot.pdf", sep = ""),
    height = 5.5, width = 6)

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

procr3 <- procrAnalysis(source = axesSta, first = dca1_3gr, second = dca2_3gr)
procr2 <- procrAnalysis(source = axesSta, first = dca1_2gr, second = dca2_2gr)
procr1 <- procrAnalysis(source = axesSta, first = dca1_1gr, second = dca2_1gr)
