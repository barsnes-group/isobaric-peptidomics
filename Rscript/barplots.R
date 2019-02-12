source("main.R")

##############################
## Bar charts for molecules and spectra

## data table for the numbers of identified spectra
specDF <- data.table(
    Spectra = c(identifications[singleAAN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(sAA[which(idn == 1)])))),
                identifications[singleAA1PTMN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(sAAm[which(idn == 1)])))),
                identifications[twoAAN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(dp[which(idn == 1)])))),
                identifications[twoAA1PTMN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(dpm[which(idn == 1)])))),
                identifications[threeAAN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(tp[which(idn == 1)])))),
                identifications[threeAA1PTMN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(tpm[which(idn == 1)])))),
                identifications[metabolitesN > 0 & nReporter > 2, .N] -
                sum(table(na.omit(unlist(met[which(idn == 1)])))),
                sum(table(na.omit(unlist(sAA[which(idn == 1)])))),
                sum(table(na.omit(unlist(sAAm[which(idn == 1)])))),
                sum(table(na.omit(unlist(dp[which(idn == 1)])))),
                sum(table(na.omit(unlist(dpm[which(idn == 1)])))),
                sum(table(na.omit(unlist(tp[which(idn == 1)])))),
                sum(table(na.omit(unlist(tpm[which(idn == 1)])))),
                sum(table(na.omit(unlist(met[which(idn == 1)]))))),
    Unique = rep(c("Not unique", "Unique"), each = 7),
    Type = factor(c("Single amino acid",
                    "Single amino acid + PTM",
                    "Dipeptide",
                    "Dipeptide + PTM",
                    "Tripeptide",
                    "Tripeptide + PTM",
                    "Metabolite",
                    "Single amino acid",
                    "Single amino acid + PTM",
                    "Dipeptide",
                    "Dipeptide + PTM",
                    "Tripeptide",
                    "Tripeptide + PTM",
                    "Metabolite"),
                  levels = rev(c("Single amino acid",
                             "Single amino acid + PTM",
                             "Dipeptide",
                             "Dipeptide + PTM",
                             "Tripeptide",
                             "Tripeptide + PTM",
                             "Metabolite"))))


spectraPlot <-
    ggplot(specDF, aes(x = Type, weight = Spectra, fill = Type)) +
    ggtitle("Spectra") +
    geom_bar(aes(alpha = Unique)) + theme_classic() +
    scale_alpha_manual(values=c(0.5, 1)) +
    theme(axis.text.y = element_text(hjust = 0.5),
          text = element_text(size = 14)) +
    coord_flip() + scale_y_continuous("Frequency", limits = c(0, 950)) +
    scale_x_discrete("") +
    scale_fill_manual(values = mycols)

emptySpectraPlot <-
    ggplot(specDF, aes(x = Type, weight = Spectra, fill = Type)) +
    ggtitle("Spectra") +
    geom_bar(aes(alpha = Unique)) + theme_classic() +
    scale_alpha_manual(values=c(0.5, 1)) +
    theme(legend.position = "none",
          text = element_text(size = 14)) +
    coord_flip() + scale_y_continuous("Frequency", limits = c(0, 1125)) +
    scale_x_discrete("", labels = NULL) +
    scale_fill_manual(values = mycols) +
    geom_text(aes(y = rep(specDF[, .(sum(Spectra)), by = Type]$V1 + 25, 2),
                  label = rep(paste(specDF[Unique == "Not unique",
                                     Spectra,
                                     by = Type]$Spectra,
                                specDF[Unique == "Unique",
                                     Spectra,
                                     by = Type]$Spectra,
                                sep = " / "), 2)
                  ), hjust = 0)



## data table for the numbers of identified molecules
idDF <- data.table(
    Identifications =
        c(length(unique(na.omit(unlist(sAA)))) -
          length(unique(na.omit(unlist(sAA[which(idn == 1)])))),
          length(unique(na.omit(unlist(sAAm)))) -
          length(unique(na.omit(unlist(sAAm[which(idn == 1)])))),
          length(unique(na.omit(unlist(dp)))) -
          length(unique(na.omit(unlist(dp[which(idn == 1)])))),
          length(unique(na.omit(unlist(dpm)))) -
          length(unique(na.omit(unlist(dpm[which(idn == 1)])))),
          length(unique(na.omit(unlist(tp)))) -
          length(unique(na.omit(unlist(tp[which(idn == 1)])))),
          length(unique(na.omit(unlist(tpm)))) -
          length(unique(na.omit(unlist(tpm[which(idn == 1)])))),
          length(unique(na.omit(unlist(met)))) -
          length(unique(na.omit(unlist(met[which(idn == 1)])))),
          length(unique(na.omit(unlist(sAA[which(idn == 1)])))),
          length(unique(na.omit(unlist(sAAm[which(idn == 1)])))),
          length(unique(na.omit(unlist(dp[which(idn == 1)])))),
          length(unique(na.omit(unlist(dpm[which(idn == 1)])))),
          length(unique(na.omit(unlist(tp[which(idn == 1)])))),
          length(unique(na.omit(unlist(tpm[which(idn == 1)])))),
          length(unique(na.omit(unlist(met[which(idn == 1)]))))),
    Unique = rep(c("Not unique", "Unique"), each = 7),
    Type = factor(c("Single amino acid",
                    "Single amino acid + PTM",
                    "Dipeptide",
                    "Dipeptide + PTM",
                    "Tripeptide",
                    "Tripeptide + PTM",
                    "Metabolite",
                    "Single amino acid",
                    "Single amino acid + PTM",
                    "Dipeptide",
                    "Dipeptide + PTM",
                    "Tripeptide",
                    "Tripeptide + PTM",
                    "Metabolite"),
                  levels = rev(c("Single amino acid",
                             "Single amino acid + PTM",
                             "Dipeptide",
                             "Dipeptide + PTM",
                             "Tripeptide",
                             "Tripeptide + PTM",
                             "Metabolite"))))

emptyIdPlot <-
    ggplot(idDF, aes(x = Type, weight = Identifications, fill = Type)) +
    ggtitle("Molecules") +
    geom_bar(aes(alpha = Unique)) + theme_classic() +
    scale_alpha_manual(values=c(0.5, 1)) +
    theme(plot.title = element_text(hjust=1),
          legend.position = "none",
          text = element_text(size = 14)) +
    coord_flip() + scale_y_reverse("Frequency", limits = c(1125, 0)) +
    scale_x_discrete("", position = "top", labels = NULL) +
    scale_fill_manual(values = mycols) +
    geom_text(aes(y = rep(idDF[, .(sum(Identifications)), by = Type]$V1 + 25, 2),
                  label = rep(paste(idDF[Unique == "Not unique",
                                     Identifications,
                                     by = Type]$Identifications,
                                idDF[Unique == "Unique",
                                     Identifications,
                                     by = Type]$Identifications,
                                sep = " / "), 2)
                  ), hjust = 1)


pdf("tablePlot.pdf", width=14, height=7)
grid.arrange(idPlot, spectraPlot, layout_matrix = matrix(c(1,1,2,2,2), nrow=1))
dev.off()

pdf("spectraPlot.pdf", width=4, height=4)
spectraPlot
dev.off()

pdf("emptySpectraPlot.pdf", width=4, height=4)
emptySpectraPlot
dev.off()

pdf("emptyIdPlot.pdf", width=4, height=4)
emptyIdPlot
dev.off()

jpeg("tablePlot.jpeg", width=14, height=14, units="in", res=600)
grid.arrange(spectraPlot, idPlot)
dev.off()

