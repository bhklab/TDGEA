# Volcano plots for RMA expression
GSE <- "GSE18552";
png(
    paste("GPL570_volcano_RMA_", GSE, ".png", sep = ""),
    width = 3000,
    height = 1800,
    res = 300
);
plot(
    x = regr_values_rma[[GSE]][,"coeff"],
    y = -log10(regr_values_rma[[GSE]][,"q"]),
    # xlim = c(-xlims_rma[[GSE]], xlims_rma[[GSE]]),
    pch = ".",
    xlab = "Regression coefficient",
    ylab = "-log10(q) (FDR adjusted)"
);
abline(h = 2, lty = 2);
dev.off()