source("main.R")

##############################
## Heat plots theoretical masses

theoreticalLevelplotMatrix <- function(dt, breakrange) {
    matrix(c(
        hist(dt[, monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "metabolite", monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "tripep+1", monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "tripep", monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "dipep+1", monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "dipep", monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "single+1", monomass],
             plot = FALSE, breaks = breakrange)$counts,
        hist(dt[type == "single", monomass],
             plot = FALSE, breaks = breakrange)$counts),
        ncol = 8, byrow = FALSE)
}

levelplotLevels <- c("Total",
                     "Metabolite",
                     "Tripeptide + PTM",
                     "Tripeptide",
                     "Dipeptide + PTM",
                     "Dipeptide",
                     "Single amino acid + PTM",
                     "Single amino acid")


breakrange <- seq(270, 1300, 10)
theoheat10 <- levelplot(
    theoreticalLevelplotMatrix(peptides, breakrange),
    col.regions = hsv(h=0.8, s=seq(0,1,length.out=1000), v=1),
    xlab = list("m/z", cex = 2), ylab = NULL,
    main = list("Theoretical masses", cex = 2.5),
    scales = list(x = list(at = which(breakrange %in% c(300, 600, 900, 1200)),
                           labels = c(300, 600, 900, 1200)),
                  y = list(at = seq(1,8),
                           labels = levelplotLevels),
                  cex = 2),
    aspect = c(1/3),
    colorkey = list(labels=list(cex=1.5)))
pdf("theoreticalheat10.pdf", height=3, width=9)
theoheat10
dev.off()
jpeg("theoreticalheat10.jpeg", width=1024, height=384)
theoheat10
dev.off()

##############################
## Heat plots experimental masses


IDinBins <- function (colNumber, colName) {
    identifications[nReporter > 2 & get(colNumber) > 0,
                    length(unique(get(colName))),
                    by = MZInterval]
}
IDTotinBins <- function () {
    identifications[nReporter > 2  &
                    (singleAAN > 0 | singleAA1PTMN > 0 |
                     twoAAN > 0 | twoAA1PTMN > 0 |
                     threeAAN > 0 | threeAA1PTMN > 0 |
                     metabolitesN > 0),
                    sum(length(unique(na.omit(singleAA))),
                        length(unique(na.omit(singleAA1mod))),
                        length(unique(na.omit(twoAA))),
                        length(unique(na.omit(twoAA1mod))),
                        length(unique(na.omit(threeAA))),
                        length(unique(na.omit(threeAA1mod))),
                        length(unique(na.omit(metabolites)))),
                    by = MZInterval]
}

breakrange <- seq(270, 1300, 10)
identifications[, MZInterval := findInterval(MZ, breakrange)]
binned <- rev(list(single = IDinBins("singleAAN", "singleAA"),
               singlePTM = IDinBins("singleAA1PTMN", "singleAA1mod"),
               two = IDinBins("twoAAN", "twoAA"),
               twoPTM = IDinBins("twoAA1PTMN", "twoAA1mod"),
               three = IDinBins("threeAAN", "threeAA"),
               threePTM = IDinBins("threeAA1PTMN", "threeAA1mod"),
               metab = IDinBins("metabolitesN", "metabolites"),
               total = IDTotinBins()))
plotmat <- matrix(0, nrow = length(breakrange), ncol = 8)
for (i in seq(8))
    plotmat[binned[[i]]$MZInterval, i] <- binned[[i]]$V1
expheat10 <- levelplot(
    plotmat,
    col.regions = hsv(h=0.875, s=seq(0,1,length.out=1000), v=1),
    xlab = list("m/z", cex = 2), ylab = NULL,
    main = list("Identified masses", cex = 2.5),
    scales = list(x = list(at = which(breakrange %in% c(300, 600, 900, 1200)),
                           labels = c(300, 600, 900, 1200)),
                  y = list(at = seq(1,8),
                           labels = levelplotLevels),
                  cex = 2),
    aspect = c(1/3),
    colorkey = list(width = 2.5, labels=list(cex=1.5)))
pdf("experimentalheat10.pdf", height=3, width=9)
expheat10
dev.off()


jpeg("heat10.jpeg", width=1024, height=768)
grid.arrange(theoheat10, expheat10, nrow=2)
dev.off()

