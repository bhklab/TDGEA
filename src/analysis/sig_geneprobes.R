### Map significant genes #######################
significance  <- 0.05;
# can do this since all dbs have the same gene probes (all on the same platform)
gene_names <- geneNames(affy_dbs_rma[[dbs[1]]]);

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