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
xlims <- list(
    "GSE9891" =  1e-7,
    "GSE18520" = 2.5e-7,
    "GSE26193" = 4e-8,
    "GSE44104" = 6e-7
);


# Volcano plot for entire collection
volcano_multi(
    input_data = q_data,
    filename = "rma_q_metagx_ovarian.png",
    title = "MetaGx Ovarian microarray probe expression",
    xlab = "Linear regression slope",
    ylab = "-log10(q) (FDR adjusted)",
    xlim = 1.25e-7
);

# Volcano plots for each dataset in collection
for (GSE in dbs) {
    volcano_single(
        input_data = q_data,
        GSE = GSE,
        filename = paste0("rma_q_metagx_ovarian_", GSE, ".png"),
        title = GSE,
        xlab = "Linear regression slope",
        ylab = "-log10(q) (FDR adjusted)",
        xlim = xlims[[GSE]],
        res = 600
    );
}