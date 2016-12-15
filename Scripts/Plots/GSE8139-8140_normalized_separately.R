dates <- c(cel_datetimes[["GSE8139"]], cel_datetimes[["GSE8140"]]);
exprs_data <- cbind(exprs(affy_dbs_rma[["GSE8139"]]), exprs(affy_dbs_rma[["GSE8140"]])); 
test_regr <- t(apply(
    exprs_data,
    1,
    function (a) {
        probe_lm <- lm(a ~ dates);
        gene_summary <- summary.lm(probe_lm);
        return(c(gene_summary$coefficients[2,4], gene_summary$coefficients[2,1]));
    }
));
test_regr <- cbind(
    test_regr,
    p.adjust(test_regr[,1], method = correction_method)
);
colnames(test_regr) <- c("p", "coeff", "q");

plot_data_separate <- data.frame(
    "coeff" = test_regr[,"coeff"],
    "q" = -log10(test_regr[,"q"])
);

png("GSE8139-8140.png", width = 3000, height = 1800, res = 300);
p1 <- (
    ggplot(
        data = plot_data_separate,
        mapping = aes(x = coeff, y = q))
    + geom_point(alpha = 0.2)
    + geom_hline(yintercept = 2, linetype = 2)
    + ggtitle("Significant probe expression for combined GSE8139-GSE8140 regression")
    + xlab("Linear regression slope")
    + ylab("-log10(q) (FDR adjusted)")
);
print(p1);
dev.off();

plot_data_together <- data.frame(
    "coeff" = regr_values_rma[["GSE8141"]][,"coeff"],
    "q" = -log10(regr_values_rma[["GSE8141"]][,"q"])
);

png("GSE8141.png", width = 3000, height = 1800, res = 300);
p1 <- (
    ggplot(
        data = plot_data_together,
        mapping = aes(x = coeff, y = q))
    + geom_point(alpha = 0.2)
    + geom_hline(yintercept = 2, linetype = 2)
    + ggtitle("Significant probe expression")
    + xlab("Linear regression slope")
    + ylab("-log10(q) (FDR adjusted)")
);
print(p1);
dev.off();