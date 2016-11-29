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
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%M-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();

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
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%M-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();

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
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%M-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();
