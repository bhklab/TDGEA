### Regression analysis #########################
print("Regression Analysis");
# p-values of linear regression having no slope
regr_values_rma <- list();
for (GSE in dbs) {
    print(GSE);
    # apply regression test to each row of matrix containing expression for a probe at the time it was taken
    regr_values_rma[[GSE]] <- t(apply(
        exprs(affy_dbs_rma[[GSE]]),
        1,
        function (a) {
            probe_lm <- lm(a ~ cel_datetimes[[GSE]]);
            gene_summary <- summary.lm(probe_lm);
            return(c(gene_summary$coefficients[2,4], gene_summary$coefficients[2,1]));
        }
    ));
    # q-values after FDR correction
    regr_values_rma[[GSE]] <- cbind(
        regr_values_rma[[GSE]],
        p.adjust(regr_values_rma[[GSE]][,1], method = correction_method)
    );
    colnames(regr_values_rma[[GSE]]) <- c("p", "coeff", "q");
}

### Map significant genes #######################
hgu_probe_map <- select(hgu133plus2.db, gene_names, c("SYMBOL", "ENTREZID", "GENENAME"));
# count the number of times a particular gene probe is ranked significant in each DB
sig_genes_rma <- sapply(
    dbs,
    function (GSE) {
        return(sum(sapply(
            gene_names,
            function (gene) {
                return(regr_values_rma[[GSE]][gene, "q"] < significance);
            }
        )));
    }
);