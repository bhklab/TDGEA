#source("https://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("affxparser");
#biocLite("affyio");
#biocLite("hgu133plus2.db")
#biocLite("sva")
#biocLite("rentrez");
#install.packages("jsonlite");

### Environment ###############################################################
library("affy");
library("affyio");
library("affyPLM");
library("data.table");
library("hgu133plus2.db");
library("jsonlite");
library("RCurl");
library("rentrez");
library("sva");
library("XML");

significance      <- 0.01;
correction_method <- "fdr";

### Functions #################################################################
# list_dirs
# Description:
#   List basenames of directories in path
# Inputs:
#   path:           Path to list directories in (default: ".")
#   pattern:        Name matching pattern (default: NULL)
#   all.dirs:       List all directories or not (default: FALSE)
#   full.names:     List full path names or just basenames (default: FALSE)
#   ignore.case:    Ignore case or not for matching (default: FALSE)
# Outputs:
#   List of strings containing directory names
list_dirs <- function(
    path        = ".",
    pattern     = NULL,
    all.dirs    = FALSE,
    full.names  = FALSE,
    ignore.case = FALSE
    ) {
    # use full.names=TRUE to pass to file.info
    all <- list.files(path, pattern, all.dirs, full.names=TRUE, recursive=FALSE, ignore.case);
    dirs <- all[file.info(all)$isdir];

    # determine whether to return full names or just dir names
    if (isTRUE(full.names)) {
        return(dirs)
    } else {
        return(basename(dirs))
    }
}

# get_cel_datetime
# Description:
#   Get CEL file's date and hour information from the header (adapated from code by Matthew McCall)
# Inputs:
#   filename:   CEL file path
# Outputs:
#   datetime:   POSIXct object
get_cel_datetime <- function(filename) {
    h <- affyio::read.celfile.header(filename, info="full");
    if  (length(h$ScanDate) > 0) {
        h$ScanDate <- gsub(pattern = "T", replacement = " ", x = h$ScanDate);
        datetime <- as.POSIXct(h$ScanDate, format = "%m/%d/%y %H:%M:%S");
        return(datetime);
    } else {
        return(NA);
    }
}

# get_cel_chip
# Description:
#   Get CEL chip type (adapated from code by Matthew McCall)
# Inputs:
#   filename:    CEL file path
# Outputs:
#   chip type
get_cel_chip <- function(filename) {
	h <- affyio::read.celfile.header(filename, info="full");
	return(as.character(h$cdfName));
}

# date_to_batch
# Description:
#   Given a list of POSIXct values, map dates together and form surrogate batches
# Inputs:
#   dates:    Input vector of dates
# Outputs:
#   Vector of batch numbers, 1-N, where N is the unique number of dates in the input
date_to_batch <- function(dates) {
    dates_formatted <- format(as.POSIXct(dates, origin = "1970-01-01"), format = "%Y-%m-%d");
    dates_unique <- unique(dates_formatted);
    
    # map formatted dates against the unique dates in the list, and number them
    return(sapply(
        dates_formatted,
        function (date) {
            return(match(date, dates_unique))
        }
    ));
}


### Main ######################################################################
### Query datasets ##############################
# search parameters (these can be changed for other queries)
gds_query_params <- c(
    "GSE[ETYP]",
    "GPL570[ACCN]",
    "\"breast\"[ALL]",
    "\"homo sapiens\"[ORGN]",
    "(MCF7[SRC] OR MCF-7[SRC])"
);
# initial data query URL
gds_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=";
gds_params <- "&retmax=5000&usehistory=y&retmode=json";
# detailed information URL from previous query
gds_url <- URLencode(paste0(
    gds_url,
    do.call("paste", as.list(c(gds_query_params, sep = "+AND+"))),
    gds_params
));
# extract query information from resulting XML 
gds_query_results <- fromJSON(getURL(gds_url));

# perform detailed query
gds_detailed_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds";
gds_detailed_url <- URLencode(paste0(
    gds_detailed_url,
    "&query_key=",
    gds_query_results$esearchresult$querykey,
    "&WebEnv=",
    gds_query_results$esearchresult$webenv,
    gds_params
));
gds_detailed_results <- fromJSON(getURL(gds_detailed_url));

# extract accession names (GSENNNN) from query results
dbs <- sapply(
    gds_detailed_results$result$uids,
    function (x) {
        return(gds_detailed_results$result[[x]]$accession)
    },
    USE.NAMES = FALSE
);

# download datasets
gds_download_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series";
# removing DBs that don't have CEL files (or multiple CEL files) available, or don't match cell line, or other criteria
dbs <- dbs[c(-48, -62, -66, -69, -82, -108, -110, -112, -113)];

# downloading proper datasets
for (GSE in dbs) {
    # need substring to form full FTP URL
    series_name <- substr(GSE, 1, nchar(GSE)-3);
    # combine strings to make FTP URL
    full_url <- paste(
        gds_download_url,
        paste0(series_name, "nnn"),
        GSE,
        "suppl",
        paste0(GSE, "_RAW.tar"),
        sep = "/"
    );
    # download data file
    download.file(
        url = full_url,
        destfile = paste0("~/Documents/TDGEA/CEL/", GSE, "_RAW.tar")
    );
}

### Preprocess data ###########################################################
print("Loading CEL files");
dbs <- c(
    "GSE6800",
    "GSE6803",
    "GSE7161",
    "GSE7327",
    "GSE8139",
    "GSE8140",
    "GSE8141",
    "GSE8565"
);
cel_files <- list();
affy_dbs <- list();
affy_dbs_rma <- list();
affy_dbs_mas5 <- list();

# Load, process, and normalize CEL files based on GSE
for (GSE in dbs) {
    print(GSE);
    cel_files[[GSE]] <- list.files(
        path = paste0("../../CEL/", GSE, "_RAW"),
        pattern = "CEL",
        full.names = TRUE,
        recursive = TRUE
    );

    # Normalize by GSE
    affy_dbs[[GSE]] <- affy::ReadAffy(filenames = cel_files[[GSE]]);    # read in to AffyBatch object
    affy_dbs_rma[[GSE]] <- affy::rma(affy_dbs[[GSE]]);                  # RMA normalize
    affy_dbs_mas5[[GSE]] <- affy::mas5(affy_dbs[[GSE]]);                # MAS5 normalize
}

# build structure for recalling datetimes based on sample names
# these calls will be made often below, so this is saving processing
gene_names <- geneNames(affy_dbs[[dbs[1]]]); # can do this since all dbs have the same genes
cel_datetimes <- list();
# get datetimes for CEL files
for (GSE in dbs) {
    cel_datetimes[[GSE]] <- sapply(cel_files[[GSE]], get_cel_datetime);
}

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