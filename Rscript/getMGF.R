library(data.table)
library(Rcpp)
sourceCpp("readMGFcpp.cpp")


mgf <- as.data.table(readMGF("../20171024_EndoCSF_TMT_Rest_Charge1.mgf", "tmt6"))

dim(mgf)
mgf

## nr ms2 scans
nrow(mgf) # 10 479
## with intact tag
sum(mgf$tag) # 7411
## with given reporter ion
mgf[tmt1 == TRUE, .N] # 7705
mgf[tmt2 == TRUE, .N] # 7584
mgf[tmt3 == TRUE, .N] # 7580
mgf[tmt4 == TRUE, .N] # 7473
mgf[tmt5 == TRUE, .N] # 7472
mgf[tmt6 == TRUE, .N] # 7460
## with n reporters
table(rowSums(mgf[, .(tmt1, tmt2, tmt3, tmt4, tmt5, tmt6)]))
##    0    1    2    3    4    5    6
## 2146  399  245  230  258  543 6658

## with all reporters and intact tag (6556)
sum(rowSums(mgf[, .(tmt1, tmt2, tmt3, tmt4, tmt5, tmt6, tag)]) == 7)

## 33 have all tmt reporters, but no intact tag
sum(rowSums(mgf[, .(tmt1, tmt2, tmt3, tmt4, tmt5, tmt6)]) == 0 & mgf[, tag])
## 102 have no reporter ions, but do have intact tag
sum(rowSums(mgf[, .(tmt1, tmt2, tmt3, tmt4, tmt5, tmt6)]) == 6 & !mgf[, tag])

## 7689 spectra with at least three reporter ions
sum(rowSums(mgf[, .(tmt1, tmt2, tmt3, tmt4, tmt5, tmt6)]) > 2)
