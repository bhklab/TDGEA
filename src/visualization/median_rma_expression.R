# Scatterplot of median RMA-normalized values
median_data <- data.frame();
for (GSE in dbs) {
    med_col <- apply(exprs(affy_dbs_rma[[GSE]]), 2, median);
    median_data <- rbind.data.frame(
        median_data,
        cbind.data.frame(
            "Date" = cel_datetimes[[GSE]],
            "Expression" = med_col,
            "Dataset" = GSE
        )
    );
}
png("GPL570_rma_median.png", width = 3000, height = 1800, res = 300);
p1 <- (
    ggplot(
        data = median_data,
        mapping = aes(x = as.POSIXct(Date, origin = "1970-01-01"), y = Expression, colour = Dataset))
    + geom_point(alpha = 0.3)
    + ggtitle("Median expression level")
    + xlab("Date")
    + ylab("Median microarray expression") 
    + ylim(c(4,6))
    + scale_x_datetime(date_labels = "%Y-%m-%d", date_breaks = "4 month")
    + theme(axis.text.x = element_text(angle = 90))
);
print(p1);
dev.off();