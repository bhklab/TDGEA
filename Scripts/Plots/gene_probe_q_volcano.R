# Volcano plots for RMA expression
xlims_rma <- list(
    "GSE3744" = 2e-6,
    "GSE12276" = 1e-7,
    "GSE18864" = 2e-7,
    "GSE19615" = 5e-8,
    "GSE20711" = 4e-6
)
xlims_log2mas5 <- list(
    "GSE3744" = 3e-6,
    "GSE12276" = 1e-7,
    "GSE18864" = 3e-7,
    "GSE19615" = 6e-8,
    "GSE20711" = 4e-6
);
for (GSE in dbs) {
    png(
        paste("GPL570_volcano_RMA_", GSE, ".png", sep = ""),
        width = 3000,
        height = 1800,
        res = 300
    );
    plot(
        x = regr_values_rma[[GSE]][,"coeff"],
        y = -log10(regr_values_rma[[GSE]][,"q"]),
        pch = ".",
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()

    png(
        paste("GPL570_volcano_log2MAS5_", GSE, ".png", sep = ""),
        width = 3000,
        height = 1800,
        res = 300
    );
    plot(
        x = regr_log2_mas5[[GSE]][,"coeff"],
        y = -log10(regr_log2_mas5[[GSE]][,"q"]),
        pch = ".",
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()
}
