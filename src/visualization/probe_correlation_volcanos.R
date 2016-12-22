### Volcano plots for correlation tests #########
q_data <- data.frame();
for (GSE in dbs) {
    q_data <- rbind.data.frame(
        q_data,
        cbind.data.frame(
            "X"       = corr_values_rma[[GSE]][,"rho"],
            "Y"       = -log10(corr_values_rma[[GSE]][,"q"]),
            "Dataset" = as.factor(GSE)
        )
    );
}

xlims <- sapply(
    dbs,
    function (GSE) {
        return(max(abs(corr_values_rma[[GSE]])));
    }
);

# Volcano plot for entire collection
volcano_multi(
    input_data = q_data,
    filename = "../reports/figures/Correlation Volcano/rma_corr_metagx_breast.png",
    title = "MetaGx Breast microarray correlation",
    xlab = "Spearman correlation",
    ylab = "-log10(q) (FDR adjusted)",
    xlim = max(xlims)
);

# Volcano plots for each dataset in collection
for (GSE in dbs) {
    volcano_single(
        input_data = q_data,
        GSE = GSE,
        filename = paste0("../reports/figures/Correlation Volcano/rma_corr_metagx_breast_", GSE, ".png"),
        title = GSE,
        xlab = "Spearman correlation",
        ylab = "-log10(q) (FDR adjusted)",
        xlim = xlims[GSE],
        res = 600
    );
}