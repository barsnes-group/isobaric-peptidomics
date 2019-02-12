## Be sure to run main.R first!

##############################
## Scatterplot


getIDtype <- function (data) {
    data <- data[, .(singleAAN, twoAAN, threeAAN,
                     singleAA1PTMN, twoAA1PTMN, threeAA1PTMN,
                     metabolitesN)]
    sapply(seq_len(nrow(data)), function (i) {
        if (sum(data[i, ]) == 0) {
            "None"
        } else if (data[i, singleAAN] > 0) {
            "Single amino acid"
        } else if (data[i, twoAAN] > 0) {
            "Dipeptide"
        } else if (data[i, threeAAN] > 0) {
            "Tripeptide"
        } else if (data[i, singleAA1PTMN] > 0) {
            "Single amino acid + PTM"
        } else if (data[i, twoAA1PTMN] > 0) {
            "Dipeptide + PTM"
        } else if (data[i, threeAA1PTMN] > 0) {
            "Tripeptide + PTM"
        } else {
            "Metabolite"
        }
    })
}

plotdt <-
    data.table(Intensity = identifications[nReporter > 2, intensity],
               mz = identifications[nReporter > 2, MZ],
               idType = getIDtype(identifications[nReporter > 2, ]))
plotdt$Identifications <-
    factor(plotdt$idType,
           levels = c("Single amino acid",
                      "Single amino acid + PTM",
                      "Dipeptide",
                      "Dipeptide + PTM",
                      "Tripeptide",
                      "Tripeptide + PTM",
                      "Metabolite", "None"))


scatterPlot <-
    ggplot(plotdt, aes(x = mz, y = Intensity,
                       col = Identifications)) +
    geom_point() +
    scale_colour_manual(values = mycols) +
    theme_classic() +
    labs(x = "m/z", y = "Intensity") +
    geom_point(data = plotdt[idType != "None", ],
               aes(x = mz, y = Intensity,
                   col = Identifications)) +
    theme(legend.position = c(.80, .70),
          plot.margin = margin(10, 10, 10, 10))

ggsave(filename="scatterplot.pdf", width=6, height=4)
ggsave(filename="scatterplot.jpeg", width=6, height=4)
