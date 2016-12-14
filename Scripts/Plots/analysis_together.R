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