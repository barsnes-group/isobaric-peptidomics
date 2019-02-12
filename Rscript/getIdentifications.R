#' Function to match masses
#' @param x One line of the processed mgf file, should be data frame/table as the mass to match is taken from the MZ column
#' @param data Data frame/table containing the masses to match against (the theoretical masses of the molecules). Masses are taken from the column monomass. Tolerance is taken from the column singleTol (this is divived by 2, as we allow tolerance on both sides). Names of the molecules are taken from the names column.
#' @return Returns a string with the names of the molecules that match the mass in x. If multiple molecules match, the names are separated by commas. If none match, it returns NA_character_
matchMass <- function(x, data) {
    a <- data$name[ abs(x$MZ - data$monomass) < (data$singleTol/2) ]
    if (length(a) == 0) NA_character_ else paste(a, collapse=",")
}
#' Function to count the number of molecules in a string made by matchMass
#' @param s the string (in our case, for example identifications[, singleAA])
#' @return number of molecules in the string.
countFun <- function (s) {
    ifelse(is.na(s), 0, str_count(s,fixed(","))+1)
}

identifications <- copy(mgf)
identifications[ , nReporter := rowSums(.SD[, paste0("tmt", 1:6)])]
identifications[ , singleAA :=
                       matchMass(.SD,
                                 data = peptides[type == "single", ]),
                by = ScanNum]
identifications[ , twoAA :=
                       matchMass(.SD,
                                 data = peptides[type == "dipep", ]),
                by = ScanNum]
identifications[ , threeAA :=
                       matchMass(.SD,
                                 data = peptides[type == "tripep", ]),
                by = ScanNum]
identifications[ , singleAA1mod :=
                       matchMass(.SD,
                                 data = peptides[type == "single+1", ]),
                by = ScanNum]
identifications[ , twoAA1mod :=
                       matchMass(.SD,
                                 data = peptides[type == "dipep+1", ]),
                by = ScanNum]
identifications[ , threeAA1mod :=
                       matchMass(.SD,
                                 data = peptides[type == "tripep+1", ]),
                by = ScanNum]
## identifications[ , twoAA2mod :=
##                        matchMass(.SD,
##                                  data = peptides[type == "dipep+2", ]),
##                 by = ScanNum]
## identifications[ , threeAA2mod :=
##                        matchMass(.SD,
##                                  data = peptides[type == "tripep+2", ]),
##                 by = ScanNum]
## identifications[ , threeAA3mod :=
##                        matchMass(.SD,
##                                  data = peptides[type == "tripep+3", ]),
##                 by = ScanNum]
identifications[ , metabolites :=
                       matchMass(.SD,
                                 data = peptides[type == "metabolite", ]),
                by = ScanNum]

identifications[ ,`:=` (singleAAN = countFun(singleAA),
                        twoAAN = countFun(twoAA),
                        threeAAN = countFun(threeAA),
                        singleAA1PTMN = countFun(singleAA1mod),
                        twoAA1PTMN = countFun(twoAA1mod),
                        threeAA1PTMN = countFun(threeAA1mod),
                        ## twoAA2PTMN = countFun(twoAA2mod),
                        ## threeAA2PTMN = countFun(threeAA2mod),
                        ## threeAA3PTMN = countFun(threeAA3mod),
                        metabolitesN = countFun(metabolites))]

