# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    require(grid);

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist);

    numPlots = length(plots);

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(
            seq(1, cols * ceiling(numPlots/cols)),
            ncol = cols,
            nrow = ceiling(numPlots/cols)
        );
    }

    if (numPlots == 1) {
        print(plots[[1]]);
    } else {
        # Set up the page
        grid.newpage();
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))));

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE));

            print(
                plots[[i]],
                vp = viewport(
                    layout.pos.row = matchidx$row,
                    layout.pos.col = matchidx$col
                )
            );
        }
    }
}

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
xlims_rma <- list(
    "GSE3744" = 2e-6,
    "GSE12276" = 1e-7,
    "GSE18864" = 2e-7,
    "GSE19615" = 5e-8,
    "GSE20711" = 4e-6
)
xlims_mas5 <- list(
    "GSE3744" = 8e-3,
    "GSE12276" = 3e-4,
    "GSE18864" = 7e-3,
    "GSE19615" = 5e-4,
    "GSE20711" = 5e-2
)
xlims_log2mas5 <- list(
    "GSE3744" = 3e-6,
    "GSE12276" = 1e-7,
    "GSE18864" = 3e-7,
    "GSE19615" = 6e-8,
    "GSE20711" = 4e-6
);
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
        xlim = c(-xlims_rma[[GSE]], xlims_rma[[GSE]]),
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()

    png(
        paste("GPL570_volcano_MAS5_", GSE, ".png", sep = ""),
        width = 3000,
        height = 1800,
        res = 300
    );
    plot(
        x = regr_values_mas5[[GSE]][,"coeff"],
        y = -log10(regr_values_mas5[[GSE]][,"q"]),
        pch = ".",
        xlim = c(-xlims_mas5[[GSE]], xlims_mas5[[GSE]]),
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()

    png(
        paste("GPL570_volcano_log2MAS5_", GSE, ".png", sep = ""),
        width = 3000,
        height = 1800,
        res = 300
    );
    plot(
        x = regr_log2_mas5[[GSE]][,"coeff"],
        y = -log10(regr_log2_mas5[[GSE]][,"q"]),
        pch = ".",
        xlim = c(-xlims_log2mas5[[GSE]], xlims_log2mas5[[GSE]]),
        xlab = "Regression coefficient",
        ylab = "-log10(q) (FDR adjusted)"
    );
    abline(h = 2, lty = 2);
    dev.off()
}
