### Spearman correlation tests ################################################
correlation_method <- "spearman";
correction_method <- "bonferroni";

print("Correlation Analysis");

# p-values of correlations between expression and time
corr_values_rma <- sapply(
    dbs,
    function (GSE) {
        print(GSE);
        # p-values and rho estimates
        corr_columns <- t(apply(
            exprs(affy_dbs_rma[[GSE]]),
            1,
            function (a) {
                # correlation test on probe expression
                test_summary <- cor.test(
                    a,
                    scaled_datetimes[[GSE]],
                    method = correlation_method
                );

                return(c(test_summary$p.value, test_summary$estimate));
            }
        ));

        # q-values after FDR correction
        corr_columns <- cbind(
            corr_columns,
            p.adjust(corr_columns[,1], method = correction_method)
        );
        colnames(corr_columns) <- c("p", "rho", "q");

        return(corr_columns);
    },
    simplify = FALSE
);