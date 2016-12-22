### Volcano plots for regressions ###############
q_data <- data.frame();
for (GSE in dbs) {
    q_data <- rbind.data.frame(
        q_data,
        cbind.data.frame(
            "X"       = regr_values_rma[[GSE]][,"coeff"],
            "Y"       = -log10(regr_values_rma[[GSE]][,"q"]),
            "Dataset" = as.factor(GSE)
        )
    );
}

xlims <- sapply(
    dbs,
    function (GSE) {
        return(max(abs(regr_values_rma[[GSE]])));
    }
);

# Volcano plot for entire collection
volcano_multi(
    input_data = q_data,
    filename = "../reports/figures/Regression Volcano/rma_regr_metagx_breast.png",
    title = "MetaGx Breast microarray probe expression",
    xlab = "Linear regression slope",
    ylab = "-log10(q) (FDR adjusted)",
    xlim = max(xlims)
);

# Volcano plots for each dataset in collection
for (GSE in dbs) {
    volcano_single(
        input_data = q_data,
        GSE = GSE,
        filename = paste0("../reports/figures/Regression Volcano/rma_regr_metagx_breast_", GSE, ".png"),
        title = GSE,
        xlab = "Linear regression slope",
        ylab = "-log10(q) (FDR adjusted)",
        xlim = xlims[GSE],
        res = 600
    );
}