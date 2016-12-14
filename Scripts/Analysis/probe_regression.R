### Regression analysis #########################
print("Regression Analysis");
# p-values of linear regression having no slope
regr_values_rma <- list();
regr_values_mas5 <- list();
regr_log2_mas5 <- list();
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

    # apply regression test to each row of matrix containing expression for a probe at the time it was taken
    regr_log2_mas5[[GSE]] <- t(apply(
        exprs(affy_dbs_mas5[[GSE]]),
        1,
        function (a) {
            probe_lm <- lm(log2(a) ~ cel_datetimes[[GSE]]); # log2(expression data), making it more comparable to RMA 
            gene_summary <- summary.lm(probe_lm);
            return(c(gene_summary$coefficients[2,4], gene_summary$coefficients[2,1]));
        }
    ));
    # q-values after FDR correction
    regr_log2_mas5[[GSE]] <- cbind(
        regr_log2_mas5[[GSE]],
        p.adjust(regr_log2_mas5[[GSE]][,1],
        method = correction_method)
    );
    colnames(regr_log2_mas5[[GSE]]) <- c("p", "coeff", "q");
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
sig_genes_log2_mas5 <- sapply(
    dbs,
    function (GSE) {
        return(sum(sapply(
            gene_names,
            function (gene) {
                return(regr_log2_mas5[[GSE]][gene, "q"] < significance);
            }
        )));
    }
);

# get info for probes that are consistently significant
probe_info <- t(sapply(
    # names of probes that are significant in 4 of the 5 datasets
    names(which(sig_genes_rma == 4)),
    function (probe_id) {
        return(hgu_probe_map[which(hgu_probe_map[,"PROBEID"] == probe_id),-1]);
    }
));