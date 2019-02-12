library(stringr)
library(ggplot2)
library(gridExtra)
library(lattice)

## define colours for bar chart and scatterplot
## mycols <- c("Single amino acid" = "#0000ff",
##            "Single amino acid + PTM" = "#7f7fff",
##            "Dipeptide" = "#ff0000",
##            "Dipeptide + PTM" = "#ff7f7f",
##            "Tripeptide" = "#ffa500",
##            "Tripeptide + PTM" = "#ffc04c",
##            "Metabolite" = "#00ff00",
##            "None" = "#cccccc")
mycols <- c("Single amino acid" = "#1b9e77",
           "Single amino acid + PTM" = "#7fc97f",
           "Dipeptide" = "#7570b3",
           "Dipeptide + PTM" = "#beaed4",
           "Tripeptide" = "#d95f02",
           "Tripeptide + PTM" = "#fdc086",
           "Metabolite" = "#fb9a99",
           "None" = "#cccccc")

source("getMGF.R")
## gives the data table "mgf"
source("makeAAdb.R")
## gives data table "peptides"
##   (concat of "singleAA", "twoAA", "threeAA", no mods)
## and "amines" for 'metabolites'
## also:
## "elemdb", "aadb" for masses of elements and residual masses of AAs
## "singleAA", "twoAA", "threeAA" for masses+TMT for these things
## "modded", "modded2", "modded3", "peptidesAndMod" modified AAs
##   and peptides
source("getIdentifications.R")
## gives data table "identifications"
colnames(identifications)
## ScanNum = scan number
## RT = retention time
## MZ = mass over charge
## (only singly charged molecules in this data set)
## intensity = intensity of the precursor
## tmt1 - tmt6 = boolean indicating whether that reporter ion found
## tag = whether the complete tag is found
## nReporter = number of reporter ions found
## singleAA = string with the name of the amino acid found
##   in the spectrum. If multiple were found, they are separated
##   with a comma. If none, this is NA.
## similar for twoAA threeAA singleAA1mod twoAA1mod threeAA1mod
##   and metabolites
## singleAAN twoAAN threeAAN singleAA1PTMN twoAA1PTMN
##   threeAA1PTMN metabolitesN = number of respective identifications

## how many molecules of each type are in the theoretical database
table(peptides$type)

## make lists of identifications
## for each type, each spectrum has a vector of strings, which are the names of the identifications (NA if spectrum is not identified as something of that type)
sAA <- strsplit(identifications[nReporter > 2, singleAA], ",", fixed=TRUE)
sAAm <- strsplit(identifications[nReporter > 2, singleAA1mod], ",", fixed=TRUE)
dp <- strsplit(identifications[nReporter > 2, twoAA], ",", fixed=TRUE)
dpm <- strsplit(identifications[nReporter > 2, twoAA1mod], ",", fixed=TRUE)
tp <- strsplit(identifications[nReporter > 2, threeAA], ",", fixed=TRUE)
tpm <- strsplit(identifications[nReporter > 2, threeAA1mod], ",", fixed=TRUE)
met <- strsplit(identifications[nReporter > 2, metabolites], ",", fixed=TRUE)

nrids <- function (x) {
    if (any(is.na(x))) 0
    else length(x)
}

## how many identifications per spectrum
## (sum of all identifications of all types)
idn <- sapply(seq_along(sAA), function (i) {
    sum(nrids(sAA[[i]]), nrids(sAAm[[i]]),
        nrids(dp[[i]]), nrids(dpm[[i]]),
        nrids(tp[[i]]), nrids(tpm[[i]]),
        nrids(met[[i]]))
})

## make csv file with all spectra with TMT (more than 2 reporter ions)
write.table(identifications[nReporter > 2, ], file = "TMTspectra.csv",
            row.names=FALSE, sep=";")

## How many things have we identified?
identifications[, .N] # 10479 spectra total
identifications[nReporter > 2, .N] # 7689 spectra >= 3 reporter ions
identifications[nReporter > 2, .N] / identifications[, .N] # 73% >= 3 reporter ions

identifications[nReporter > 2 &
                sum(c("singleAAN":"threeAA1PTMN")) > 0, .N]
identifications[nReporter > 2 &
                (singleAAN > 0 | singleAA1PTMN > 0 |
                 twoAAN > 0 | twoAA1PTMN > 0 |
                 threeAAN > 0 | threeAA1PTMN > 0), .N] # 2220 identified not counting metabolites
identifications[nReporter > 2 &
                (singleAAN > 0 | singleAA1PTMN > 0 |
                 twoAAN > 0 | twoAA1PTMN > 0 |
                 threeAAN > 0 | threeAA1PTMN > 0 |
                 metabolitesN > 0), .N] # 2265 identified
identifications[nReporter > 2 &
                (singleAAN > 0 | singleAA1PTMN > 0 |
                 twoAAN > 0 | twoAA1PTMN > 0 |
                 threeAAN > 0 | threeAA1PTMN > 0 |
                 metabolitesN > 0), .N] /
    identifications[nReporter > 2, .N]  # 29% of TMT containing spectra identified
identifications[nReporter > 2 &
                (singleAAN > 0 | singleAA1PTMN > 0 |
                 twoAAN > 0 | twoAA1PTMN > 0 |
                 threeAAN > 0 | threeAA1PTMN > 0 |
                 metabolitesN > 0), .N] /
    identifications[, .N]  # 22% of all spectra identified

identifications[nReporter > 2 &
                (singleAAN == 0 & singleAA1PTMN == 0 &
                 twoAAN == 0 & twoAA1PTMN == 0 &
                 threeAAN == 0 & threeAA1PTMN == 0 &
                 metabolitesN == 0), .N] # 5424 with TMT not identified

sum(idn > 0) # 2265 spectra with at least one identification
sum(idn == 1) # 1384 spectra with exactly one identification
sum(idn == 1) / length(idn) # 18 % of spectra with TMT uniquely identified
sum(idn == 1) / sum(idn > 0) # 61 % of identified spectra uniquely identified
sum(idn > 1) # 881 spectra with multiple identifications


## which amino acids/peptides are identified
unique(na.omit(unlist(sAA)))   #  17
unique(na.omit(unlist(sAAm)))  #   2
unique(na.omit(unlist(dp)))    # 103
unique(na.omit(unlist(dpm)))   #  52
unique(na.omit(unlist(tp)))    # 440
unique(na.omit(unlist(tpm)))   # 355
unique(na.omit(unlist(met)))   #  21

## nr molecules identified: 990
sum(length(unique(na.omit(unlist(sAA)))),
    length(unique(na.omit(unlist(sAAm)))),
    length(unique(na.omit(unlist(dp)))),
    length(unique(na.omit(unlist(dpm)))),
    length(unique(na.omit(unlist(tp)))),
    length(unique(na.omit(unlist(tpm)))),
    length(unique(na.omit(unlist(met)))))

## nr of molecules (bottom nr) associated with nr of spectra (top nr)
table(table(c(
    na.omit(unlist(sAA)),
    na.omit(unlist(sAAm)),
    na.omit(unlist(dp)),
    na.omit(unlist(dpm)),
    na.omit(unlist(tp)),
    na.omit(unlist(tpm)),
    na.omit(unlist(met)))))

## pct of molecules (bottom nr) asociated with (top nr) nr of spectra
## 35% of molecules associated with 1 spectrum
## 0.1% of molecules (1) associated with 300 spectra
table(table(c(
    na.omit(unlist(sAA)),
    na.omit(unlist(sAAm)),
    na.omit(unlist(dp)),
    na.omit(unlist(dpm)),
    na.omit(unlist(tp)),
    na.omit(unlist(tpm)),
    na.omit(unlist(met))))) / sum (table(table(c(
    na.omit(unlist(sAA)),
    na.omit(unlist(sAAm)),
    na.omit(unlist(dp)),
    na.omit(unlist(dpm)),
    na.omit(unlist(tp)),
    na.omit(unlist(tpm)),
    na.omit(unlist(met)))))
)

## number of scans associated with each molecule
sort(
    table(c(na.omit(unlist(sAA)),
            na.omit(unlist(sAAm)),
            na.omit(unlist(dp)),
            na.omit(unlist(dpm)),
            na.omit(unlist(tp)),
            na.omit(unlist(tpm)),
            na.omit(unlist(met)))))
## pct of scans associated with each molecule
## nr of scans for each molecule, divided by total nr of identified scans
sort(
    table(c(na.omit(unlist(sAA)),
            na.omit(unlist(sAAm)),
            na.omit(unlist(dp)),
            na.omit(unlist(dpm)),
            na.omit(unlist(tp)),
            na.omit(unlist(tpm)),
            na.omit(unlist(met)))) / sum(idn > 0) )


## pct top 4 scans
sum(table(c(na.omit(unlist(sAA))))[c("Tryptophan", "Tyrosine",
                                     "Phenylalanine", "Leucine")])/
    sum (idn > 0)

## which amino acids/peptides uniquely identify at least one spectrum
unique(na.omit(unlist(sAA[which(idn == 1)])))  # 14
unique(na.omit(unlist(sAAm[which(idn == 1)]))) #  1
unique(na.omit(unlist(dp[which(idn == 1)])))   # 52
unique(na.omit(unlist(dpm[which(idn == 1)])))  # 15
unique(na.omit(unlist(tp[which(idn == 1)])))   # 92
unique(na.omit(unlist(tpm[which(idn == 1)])))  # 56
unique(na.omit(unlist(met[which(idn == 1)])))  # 11

## total number of molecules that uniquely identify a spectrum 241
sum(length(unique(na.omit(unlist(sAA[which(idn == 1)])))),
    length(unique(na.omit(unlist(sAAm[which(idn == 1)])))),
    length(unique(na.omit(unlist(dp[which(idn == 1)])))),
    length(unique(na.omit(unlist(dpm[which(idn == 1)])))),
    length(unique(na.omit(unlist(tp[which(idn == 1)])))),
    length(unique(na.omit(unlist(tpm[which(idn == 1)])))),
    length(unique(na.omit(unlist(met[which(idn == 1)])))))


## how many spectra are uniquely identified by each ...
sum(table(na.omit(unlist(sAA[which(idn == 1)]))))  # 761
sum(table(na.omit(unlist(sAAm[which(idn == 1)])))) #   4
sum(table(na.omit(unlist(dp[which(idn == 1)]))))   # 200
sum(table(na.omit(unlist(dpm[which(idn == 1)]))))  #  32
sum(table(na.omit(unlist(tp[which(idn == 1)]))))   # 221
sum(table(na.omit(unlist(tpm[which(idn == 1)]))))  # 122
sum(table(na.omit(unlist(met[which(idn == 1)]))))  #  44

## number of molecules (bottom number in table) that uniquely identify this number (top number in table) of spectra
table(table(c(
    na.omit(unlist(sAA[which(idn == 1)])),
    na.omit(unlist(sAAm[which(idn == 1)])),
    na.omit(unlist(dp[which(idn == 1)])),
    na.omit(unlist(dpm[which(idn == 1)])),
    na.omit(unlist(tp[which(idn == 1)])),
    na.omit(unlist(tpm[which(idn == 1)])),
    na.omit(unlist(met[which(idn == 1)])))))


## number of types of identification per spectrum
idtn <- sapply(seq_along(sAA), function (i) {
    sum(sapply(c(nrids(sAA[[i]]), nrids(sAAm[[i]]),
                 nrids(dp[[i]]), nrids(dpm[[i]]),
                 nrids(tp[[i]]), nrids(tpm[[i]]),
                 nrids(met[[i]])), `>`, 0))
})

