### SVA/ComBat correction #######################
pheno_rma <- sapply(
    dbs,
    function (GSE) {
        val <- cbind(
            pData(affy_dbs_rma[[GSE]]),
            format(as.POSIXct(cel_datetimes[[GSE]], origin = "1970-01-01"), format = "%Y-%m-%d"),
            format(as.POSIXct(cel_datetimes[[GSE]], origin = "1970-01-01"), format = "%H:%M:%S")
        );
        colnames(val) <- c("Sample", "Date", "Time", "Batch");
        return(val);
    },
    simplify = FALSE
);
model_rma_null <- sapply(
    dbs,
    function (GSE) {
        return(model.matrix(~1, data=pheno_rma[[GSE]]));
    },
    simplify = FALSE
);
combat_rma_edata <- sapply(
    dbs,
    function (GSE) {
        return(ComBat(
            dat = exprs(affy_dbs_rma[[GSE]]),
            batch = pheno_rma[[GSE]]$Date,
            mod = model_rma_null[[GSE]],
            par.prior = TRUE,
            prior.plots = FALSE
        ));
    },
    simplify = FALSE
);
combat_rma_pvalues <- sapply(
    dbs,
    function (GSE) {
        return(f.pvalue(combat_rma_edata[[GSE]], model_rma_full[[GSE]], model_rma_null[[GSE]]));
    },
    simplify = FALSE
);
combat_rma_qvalues <- sapply(
    dbs,
    function (GSE) {
        return(p.adjust(combat_rma_pvalues[[GSE]], method = correction_method));
    },
    simplify = FALSE
);

model_rma_full <- sapply(
    dbs,
    function (GSE) {
        return(model.matrix(~as.factor(Date), data=pheno_rma[[GSE]]));
    },
    simplify = FALSE
);
model_rma_nsv <- sapply(
    dbs,
    function (GSE) {
        return(num.sv(exprs(affy_dbs_rma[[GSE]]), model_rma_full[[GSE]], method = "leek"));
    },
    simplify = FALSE
);
model_rma_sva <- sapply(
    dbs,
    function (GSE) {
        return(sva(
            exprs(affy_dbs_rma[[GSE]]),
            model_rma_full[[GSE]],
            model_rma_null[[GSE]],
            n.sv = model_rma_nsv[[GSE]]
        ));
    },
    simplify = FALSE
);