### Regression analysis #########################
correction_method <- "bonferroni";

print("Regression Analysis");

# p-values of linear regression having no slope
regr_values_rma <- sapply(
    dbs,
    function (GSE) {
        print(GSE);
        # apply regression test to each row of matrix containing expression for a probe at the time it was taken
        regr_columns <- t(apply(
            exprs(affy_dbs_rma[[GSE]]),
            1,
            function (a) {
                probe_lm <- lm(a ~ scaled_datetimes[[GSE]]);
                gene_summary <- summary.lm(probe_lm);
                return(c(gene_summary$coefficients[2,4], gene_summary$coefficients[2,1]));
            }
        ));
        # q-values after FDR correction
        regr_columns <- cbind(
            regr_columns,
            p.adjust(regr_columns[,1], method = correction_method)
        );
        colnames(regr_columns) <- c("p", "coeff", "q");

        return(regr_columns);
    },
    simplify = FALSE
);