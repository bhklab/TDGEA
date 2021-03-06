# Separate normalization
cel_rma <- matrix(, ncol = 10, nrow = 54675);
for (i in seq(1, 10)) {
    data_rma <- affy::justRMA(filenames = cel_files[i]);
    cel_rma[,i] <- exprs(data_rma);
}
colnames(cel_rma) <- cel_table$Sample[1:10];

png("RMA_separate.png", width = 3000, height = 1800, res = 300);
boxplot(
    cel_rma,
    main = "Sample expression, normalized separately",
    las = 2,
    cex.axis = 0.75
);
dev.off();

# Grouped normalization
cel_rma <- affy::justRMA(filenames = cel_files[1:10]);
data_plots <- data.table(exprs(cel_rma));
names(data_plots) <- cel_table$Sample[1:10];
png("RMA_together.png", width = 3000, height = 1800, res = 300);
boxplot(
    data_plots,
    main = "Sample expression, normalized together",
    las = 2,
    cex.axis = 0.75
);
dev.off();