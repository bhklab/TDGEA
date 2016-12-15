# Volcano plots for RMA expression
xlims_rma <- list(
    "GSE6800" = 1.5e-7,
    "GSE6803" = 1.7e-7,
    "GSE7161" = 1.5e-3,
    "GSE7327" = 1.2e-6,
    "GSE8139" = 2e-7,
    "GSE8140" = 1e-6,
    "GSE8141" = 2.5e-7,
    "GSE8565" = 6e-7
)

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
        # xlim = c(-xlims_rma[[GSE]], xlims_rma[[GSE]]),
        pch = ".",
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()
}
