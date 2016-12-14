# Scatterplots of individual gene expression over time
length_out <- 10;
delta <- (max(cel_datetimes, na.rm = TRUE) - min(cel_datetimes, na.rm = TRUE))/length_out;
ix <- seq(min(cel_datetimes, na.rm = TRUE) - delta, max(cel_datetimes, na.rm = TRUE), length.out = length_out + 1);
for (gene in gene_names) {
    x <- gene_exprs[[gene]]$DateTime;
    y <- gene_exprs[[gene]]$Expression;
    delta <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/length_out;
    ix <- seq(min(x, na.rm = TRUE) - delta, max(x, na.rm = TRUE), length.out = length_out + 1);

    png(paste("../../Figures/Expression/expression_", gene, ".png", sep = ""), width = 3000, height = 1800, res = 300);
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