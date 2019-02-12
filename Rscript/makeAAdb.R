library("XML")
library("data.table")
############################################################
## Get modification names and masses from unimod
############################################################
unimod <- xmlParse("../unimod.xml")
## xmlSApply(xmlRoot(unimod), xmlName)
## get the different parts of the unimod xml file
modifications <- xmlRoot(unimod)[[2]]
aminoacids <- xmlRoot(unimod)[[3]]
elements <- xmlRoot(unimod)[[1]]

## make a data.table from elements (as in the grand table of elements)
elemdb <-
    rbindlist(
        lapply(seq_len(xmlSize(elements)), function (i) {
            data.table(
                title = xmlAttrs(elements[[i]])["title"],
                name = xmlAttrs(elements[[i]])["full_name"],
                monomass = as.numeric(xmlAttrs(elements[[i]])["mono_mass"])
            )
        })
    )

## make a data.table from the 20 amino acids
## leave out Isoleucine, as for this purpose same as Leucine
## monomass here is residual mass
aadb <-
    rbindlist(
        lapply(seq_len(xmlSize(aminoacids)), function (i) {
            if (nchar(xmlAttrs(aminoacids[[i]])["title"]) == 1 &&
                xmlAttrs(aminoacids[[i]])["title"] != "U" &&
                xmlAttrs(aminoacids[[i]])["title"] != "I" &&
                xmlAttrs(aminoacids[[i]])["title"] != "-") {
                data.table(
                    title = xmlAttrs(aminoacids[[i]])["title"],
                    abbrv = xmlAttrs(aminoacids[[i]])["three_letter"],
                    name = xmlAttrs(aminoacids[[i]])["full_name"],
                    residmass =
                        sum(xmlSApply(aminoacids[[i]], function (x) {
                            m <- xmlAttrs(x)
                            elemdb[title == m["symbol"],
                                   monomass] * as.integer(m["number"])
                        }))
                )
            }
        })
    )

## add carbamidomethylation of C and oxidation of M as normal amino acids
carbaModi <- which(sapply(seq_len(xmlSize(modifications)), function (i) {
    xmlAttrs(modifications[[i]])["title"] == "Carbamidomethyl"
}))
carbaModMass <- as.numeric(xmlGetAttr(
    modifications[[carbaModi]][["delta"]], "mono_mass"))
aadb <-
    rbindlist(list(aadb,
                   list(title = "MO",
                        abbrv = "MetO",
                        name = "Methionine+O",
                        residmass = aadb[name == "Methionine", residmass] +
                            elemdb[name == "Oxygen", monomass]),
                   list(title = "CCAM",
                        abbrv = "CysCAM",
                        name = "Cysteine+CAM",
                        residmass = aadb[name == "Cysteine", residmass] +
                            carbaModMass)))

## terminals <- sum(ends[, monomass])
water <- elemdb[title == "O", monomass] + 2*elemdb[title == "H", monomass]
## charge is one hydrogen minus one electron
charge <- elemdb[title == "H", monomass]-elemdb[title == "e", monomass]
tmti <- which(sapply(seq_len(xmlSize(modifications)), function (i) {
    startsWith(xmlAttrs(modifications[[i]])["title"], "TMT6")
}))
tmt <- as.numeric(xmlGetAttr(
    modifications[[tmti]][["delta"]], "mono_mass"))

## single amino acids for matching
## just name, monomass, and elements
## monomass == the mass as we would find in the mass spec!
##   so: residual mass + water + charge + tmt
singleAA <- copy(aadb)
singleAA[, monomass := residmass + water + charge + tmt]
singleAA[, c("title", "abbrv", "residmass"):=NULL]

## dipeptides
## excluding reverses (Alanine+Arginine is included, Arginine+Alanine not)
## allowing doubles, 231 combinations (r+n-1)!/(r!(n-1)!) with r=2, n=21
twoAA <- rbindlist(
    lapply(seq_len(nrow(aadb)), function (i) {
        rbindlist(
            lapply(seq(i, nrow(aadb)), function (j) {
                list(name = paste(aadb[i, name], aadb[j, name], sep="::"),
                     monomass = aadb[i, residmass] + aadb[j, residmass] +
                         water + charge + tmt)
            })
        )
    })
)

## tripeptides
## allowing doubles and triplets,
## 1771 combinations (r+n-1)!/(r!(n-1)!) with r=3, n=21
threeAA <- rbindlist(
    lapply(seq_len(nrow(aadb)), function (i) {
        rbindlist( lapply(seq(i, nrow(aadb)), function (j) {
            rbindlist( lapply(seq(j, nrow(aadb)), function (k) {
                list(name = paste(aadb[i, name],
                                  aadb[j, name],
                                  aadb[k, name],
                                  sep="::"),
                     monomass = aadb[i, residmass] + aadb[j, residmass] +
                         aadb[k, residmass] + water + charge + tmt)
            })
            )
        })
        )
    })
)

## put one two three amino acids in one data table
## 21 single AA (20, minus Isoleucine, plus Meth+O and Cyst+CAM)
## 231 dipeptides
## 1771 tripeptides
## sum = 2023
peptides <-
    rbindlist(list(single = singleAA,
                   dipep = twoAA,
                   tripep = threeAA),
              idcol = "type")


### MODIFICATIONS
## make data table of all modifications
moddt <-
    rbindlist(
        lapply(seq_len(xmlSize(modifications)), function (i) {
            modName <- xmlAttrs(modifications[[i]])["full_name"]
            modTitle <- xmlAttrs(modifications[[i]])["title"]
            AAs <- which(
                xmlSApply(modifications[[i]], xmlName) == "specificity")
            specificities <-
                sapply(AAs, function (x) {
                    xmlAttrs(modifications[[i]][[x]])[
                        c("hidden", "site", "position", "classification")]
                })
            modMass <- xmlAttrs(modifications[[i]][["delta"]])["mono_mass"]
            data.table(
                name = modName,
                title = modTitle,
                common = ifelse(specificities["hidden", ] == "1", FALSE, TRUE),
                site = specificities["site", ],
                position = specificities["position", ],
                classification = specificities["classification", ],
                mass = as.numeric(modMass))
        })
    )

## moddt[common == TRUE,] # 80 common mods
## moddt[common == TRUE & classification != "Isotopic label", ]
## moddt[common == TRUE &
##       classification != "Isotopic label" &
##       site != "N-term", ]
## moddt[common == TRUE &
##       site != "N-term" &
##       position != "Protein C-term" &
##       classification == "Post-translational", ] # 11 common & PTM & !N-term & !Protein C-term
## moddt[site != "N-term" &
##       position != "Protein C-term" &
##       classification == "Post-translational", ] # 225 PTM & !N-term & !Protein C-term

## moddt[common == FALSE, .N] # 2365 uncommon mods
## moddt[common == FALSE &
##       site != "N-term" &
##       position != "Protein C-term" & position != "Protein N-term" &
##       position != "Any N-term" &
##       classification == "Post-translational", ] # 200 uncommon & PTM & !N-term & !Protein C-term
## unique(moddt[common == FALSE &
##       site != "N-term" &
##       position != "Protein C-term" & position != "Protein N-term" &
##       position != "Any N-term" &
##       classification == "Post-translational", site]) #

## common ptms, not on N-term or protein C-term
##  1:   Biotinylation      Biotin    K
##  2: Phosphorylation     Phospho    Y
##  3: Phosphorylation     Phospho    T
##  4: Phosphorylation     Phospho    S
##  5:     Methylation      Methyl    E
##  6:     Methylation      Methyl    D
##  7:   O-Sulfonation       Sulfo    S
##  8:   O-Sulfonation       Sulfo    T
##  9:   O-Sulfonation       Sulfo    Y
## 10:       dihydroxy Dioxidation    M
## 11:   Crotonylation    Crotonyl    K
cptmdt <- copy(moddt[common == TRUE &
                     site != "N-term" &
                     position != "Protein C-term" &
                     classification == "Post-translational", ])

## LOOKING FOR ALL COMMON, NOT N-TERM MODIFICATIONS (except isotopic labels)
## cptmdt <- copy(moddt[common == TRUE &
##                      classification != "Isotopic label" &
##                      site != "N-term", ])

## site == one amino acid, position == Anywhere so can just add ptm to
## each available peptide.
## no reverse of the peptide, so specific AA, but C-term can be added to all
makeOneModName <- function (peptide, mod, aa) {
    sub(paste0(aa, "(?!\\+)"), paste0(aa, "+", mod),
        peptide, perl=TRUE)
}

## 1 modification
## 2783 possiblilities for one modification (11 common, site specifc, anywhere)
modded <-
    rbindlist(
        lapply(seq_len(cptmdt[, .N]), function (i) {
            amac <- aadb[title == cptmdt[i, site], name]
            temp <- peptides[grepl(paste0(amac, "(?!\\+)"),
                                   peptides$name, perl=TRUE),]
            data.table(type = paste0(temp$type, "+1"),
                       name = makeOneModName(temp$name,
                                             cptmdt[i, title],
                                             amac),
                       monomass = temp$monomass +
                           cptmdt[i, mass])
        }))
## table(modded$type) # 11 singleAA, 231 dipep, 2541 tripep = 2783

## DO THIS INSTEAD WHEN
## LOOKING FOR ALL COMMON, NOT N-TERM MODIFICATIONS
## modded <-
##     rbindlist(
##         lapply(seq_len(cptmdt[site != "C-term", .N]), function (i) {
##             amac <- aadb[title == cptmdt[site != "C-term", ][i, site], name]
##             temp <- peptides[grepl(paste0(amac, "(?!\\+)"),
##                                    peptides$name, perl=TRUE),]
##             data.table(type = paste0(temp$type, "+1"),
##                        name = makeOneModName(temp$name,
##                                              cptmdt[site != "C-term", ][i, title],
##                                              amac),
##                        monomass = temp$monomass +
##                            cptmdt[site != "C-term", ][i, mass])
##         }))
## moddedC <-
##     rbindlist(
##         lapply(seq_len(cptmdt[site == "C-term", .N]), function (i) {
##             peps <- copy(peptides[type %in% c("single", "dipep", "tripep") &
##                                   !grepl("+", peptides$name, fixed=TRUE), ])
##             data.table(type = paste0(peps$type, "+1"),
##                        name = paste0(peps$name, "+",
##                                      cptmdt[site == "C-term", ][i, title]),
##                        monomass = peps$monomass +
##                            cptmdt[site == "C-term", ][i, mass])
##         }))


## 2 modifications
## 2662 possiblilities for two modifications
## only 1540 unique names, though...
## modded2 <-
##     rbindlist(
##         lapply(seq_len(nrow(cptmdt)), function (i) {
##             amac <- aadb[title == cptmdt[i, site], name]
##             temp <- modded[grepl(paste0(amac, "(?!\\+)"),
##                                  modded$name, perl=TRUE),]
##             data.table(type = sub("+1", "+2", temp$type, fixed=TRUE),
##                        name = makeOneModName(temp$name,
##                                              cptmdt[i, title],
##                                              amac),
##                        monomass = temp$monomass + cptmdt[i, mass])
##         }))
## table(modded2$type) # 121 dipep, 2541 tripep = 2662

## 2 modifications
## 1331 possiblilities for two modifications
## only 338 unique names, though...
## modded3 <-
##     rbindlist(
##         lapply(seq_len(nrow(cptmdt)), function (i) {
##             amac <- aadb[title == cptmdt[i, site], name]
##             temp <- modded2[grepl(paste0(amac, "(?!\\+)"),
##                                   modded2$name, perl=TRUE),]
##             data.table(type = sub("+2", "+3", temp$type, fixed=TRUE),
##                        name = makeOneModName(temp$name,
##                                              cptmdt[i, title],
##                                              amac),
##                        monomass = temp$monomass + cptmdt[i, mass])
##         }))
## ## table(modded3$type) # 1331 tripep


## put one two three amino acids and 1 2 3 ptms in one data table
## 21 single AA (20, minus Isoleucine, plus Meth+O and Cyst+CAM)
## 231 dipeptides
## 1771 tripeptides
## 2783 one ptm
## 2662 two ptm
## 1331 thee ptm
## sum = 8799
## for multiple modifications
## peptides <-
##     rbindlist(list(peptides, modded, modded2, modded3))
## for all common mods
## peptides <-
##     rbindlist(list(peptides, modded, moddedC))
peptides <-
    rbindlist(list(peptides, modded))


hmdb.csf <- xmlParse("../csf_metabolites.xml")
## xmlSApply(xmlRoot(hmdb.csf), xmlName)
metabolites <- xmlRoot(hmdb.csf)

subclasses <- xmlSApply(metabolites, function (m)
    xmlValue(m[["taxonomy"]][["sub_class"]]))
## table(subclasses)
##                               Amines
##                                   10
## Amino acids, peptides, and analogues
##                                   75
aminesi <- which(subclasses %in% c("Amines", "Amino acids, peptides, and analogues"))
sapply(aminesi, function (i) xmlValue(metabolites[[i]][["name"]]))
## these have a water more than the AA from unimod!!!
amines <-
    data.table(
        accession =
            unlist(sapply(aminesi, function (i)
                xmlValue(metabolites[[i]][["accession"]]))),
        name =
            unlist(sapply(aminesi, function (i)
                xmlValue(metabolites[[i]][["name"]]))),
        chemicalFormula =
            unlist(sapply(aminesi, function (i)
                xmlValue(metabolites[[i]][["chemical_formula"]]))),
        monoisotopicMolecularWeight =
            as.numeric(unlist(sapply(aminesi, function (i)
                xmlValue(metabolites[[i]][["monisotopic_molecular_weight"]]))))
    )

## singleAA$name
## amines that are amino acids
## sum(amines$name %in% singleAA$name)
## amines$name[which(amines$name %in% singleAA$name)]
## amines[name == "Glycine", ]
## singleAA[name == "Glycine", ]
## amines that are a specific version of amino acids
## sum(amines$name %in% paste0("L-", singleAA$name))
## amines$name[which(amines$name %in% paste0("L-", singleAA$name))]
## amines[amines$name %in% paste0("L-", singleAA$name), .(name, monomass)]
## singleAA[, .(name, monomass)]
## sum(amines$name %in% paste0("D-", singleAA$name))
## amines$name[which(amines$name %in% paste0("D-", singleAA$name))]
## amines[amines$name %in% paste0("D-", singleAA$name), .(name, monomass)]
## singleAA[, .(name, monomass)]
## remove amines that are in single amino acid db above
amines <- amines[! name %in% c(singleAA$name,
                               paste0("L-", singleAA$name),
                               paste0("D-", singleAA$name),
                               "L-Isoleucine",
                               "L-Alloisoleucine")]
amines[name == "N6,N6,N6-Trimethyl-L-lysine", name := "N6-N6-N6-Trimethyl-L-lysine"]

amines[, monomass := monoisotopicMolecularWeight + charge + tmt]
amines[, type := "metabolite"]

peptides <-
    rbindlist(list(peptides, amines[, .(type, name, monomass)]))
## tolerance of 10ppm
peptides[, singleTol := monomass*1e-5]
