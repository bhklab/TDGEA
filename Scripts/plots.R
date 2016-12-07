# Histogram of dates
length_out <- 10;
delta <- (max(cel_datetimes, na.rm = TRUE) - min(cel_datetimes, na.rm = TRUE))/length_out;
ix <- seq(min(cel_datetimes, na.rm = TRUE) - delta, max(cel_datetimes, na.rm = TRUE), length.out = length_out + 1);
png("GPL570_dates_histogram.png", width = 3000, height = 1800, res = 300);
hist(
    cel_datetimes,
    xaxt = "n",
    xlab = "",
    main = "Mean expression level of GPL570 chips"
);
axis(
    1,
    at = ix,
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%m-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();

# Scatterplot of mean RMA-normalized values
png("GPL570_rma_mean.png", width = 3000, height = 1800, res = 300);
plot(
    data.table(cel_table$DateTime, cel_table$RMA_Mean),
    xaxt = "n",
    xlab = "",
    main = "Mean expression level of GPL570 chips"
);
axis(
    1,
    at = ix,
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%m-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();

# Scatterplot of median RMA-normalized values
png("GPL570_rma_median.png", width = 3000, height = 1800, res = 300);
plot(
    data.table(cel_table$DateTime, cel_table$RMA_Mean),
    xaxt = "n",
    xlab = "",
    main = "Median expression level of GPL570 chips"
);
axis(
    1,
    at = ix,
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%m-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();

# Scatterplots of individual gene expression over time
for (gene in gene_names) {
    x <- gene_exprs[[gene]]$DateTime;
    y <- gene_exprs[[gene]]$Expression;
    delta <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/length_out;
    ix <- seq(min(x, na.rm = TRUE) - delta, max(x, na.rm = TRUE), length.out = length_out + 1);

    png(paste("../Figures/Expression/expression_", gene, ".png", sep = ""), width = 3000, height = 1800, res = 300);
    plot(
        x = x,
        y = y,
        xaxt = "n",
        xlab = "Date",
        ylab = "Normalized Gene Expression",
        ylim = c(0, 15),
        main = paste("Gene expression over time for", gene)
    );
    axis(
        1,
        at = ix,
        labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%m-%d"),
        las = 2,
        cex.axis = 0.75
    );
    dev.off();
}

# Volcano plots for RMA expression
xlims <- list(
    "GSE3744" = 2e-6,
    "GSE12276" = 1e-7,
    "GSE18864" = 2e-7,
    "GSE19615" = 5e-8,
    "GSE20711" = 4e-6
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
        pch = ".",
        xlim = c(-xlims[[GSE]], xlims[[GSE]]),
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()
}
