#source("https://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("affxparser");
#biocLite("affyio");
#biocLite("hgu133plus2.db")

### Environment ###############################################################
library("affy");
library("affyio");
library("affyPLM");
library("data.table");
library("hgu133plus2.db");

# From DOI: 10.1186/1471-2105-6-214
sigma <- list(
    "dCHIP"    = 0.469,
    "MAS5"     = 0.475,
    "RMA"      = 0.467,
    "GCRMA-EB" = 0.459
);

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
list_dirs <- function(path=".", pattern=NULL, all.dirs=FALSE, full.names=FALSE, ignore.case=FALSE) {
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


### Main ######################################################################
# Organize data for data table
batches <- c("GSE19615", "GSE18864", "GSE3744", "GSE12276", "GSE20711");
cel_files <- list();
affy_batches <- list();
affy_batches_rma <- list();
affy_batches_mas5 <- list();

# Load, process, and normalize CEL files based on batch
for (batch in batches) {
    cel_files[[batch]] <- list.files(
        path = paste("../CEL/", batch, "_RAW", sep = ""),
        pattern = "CEL",
        full.names = TRUE,
        recursive = TRUE
    );

    # Normalize by batch
    affy_batches[[batch]] <- affy::ReadAffy(filenames = cel_files[[batch]]);    # read in to AffyBatch object
    affy_batches_rma[[batch]] <- affy::rma(affy_batches[[batch]]);              # RMA normalize
    affy_batches_mas5[[batch]] <- affy::mas5(affy_batches[[batch]]);            # MAS5 normalize
}

# build structure for recalling datetimes based on sample names
# these calls will be made often below, so this is saving processing
gene_names <- geneNames(affy_batches[[batches[1]]]); # can do this since all batches have the same genes
cel_datetimes <- list();
# get datetimes for CEL files
for (batch in batches) {
    cel_datetimes[[batch]] <- sapply(cel_files[[batch]], get_cel_datetime);
}

# Organize expression of each gene probe by DateTime
gene_exprs_rma <- list();
gene_exprs_mas5 <- list();
for (batch in batches) {
    for (gene in gene_names) {
        # RMA normalized gene expression
        gene_exprs_rma[[gene]] <- rbind(
            gene_exprs_rma[[gene]],
            data.table(
                "DateTime" = cel_datetimes[[batch]],
                "Expression" = exprs(affy_batches_rma[[batch]])[gene,]
            )
        );

        # MAS5 normalized gene expression
        gene_exprs_mas5[[gene]] <- rbind(
            gene_exprs_mas5[[gene]],
            data.table(
                "DateTime" = cel_datetimes[[batch]],
                "Expression" = exprs(affy_batches_mas5[[batch]])[gene,]
            )
        );
    }
}

# p-values of linear regression having no slope
p_values <- data.table(
    "RMA" = lapply(
        gene_names,
        function(gene) {
            gene_lm <- lm(gene_exprs_rma[[gene]]$Expression ~ gene_exprs_rma[[gene]]$DateTime);
            gene_summary <- summary.lm(gene_lm);
            return(pf(
                gene_summary$fstatistic[1],
                gene_summary$fstatistic[2],
                gene_summary$fstatistic[3],
                lower.tail = FALSE)
            );
        }
    ),
    "MAS5" = lapply(
        gene_names,
        function(gene) {
            gene_lm <- lm(gene_exprs_mas5[[gene]]$Expression ~ gene_exprs_mas5[[gene]]$DateTime);
            gene_summary <- summary.lm(gene_lm);
            return(pf(
                gene_summary$fstatistic[1],
                gene_summary$fstatistic[2],
                gene_summary$fstatistic[3],
                lower.tail = FALSE)
            );
        }
    )
);

# q-values after FDR correction
q_values <- data.table(
    "RMA" = p.adjust(p_values$RMA, method = "fdr"),
    "MAS5" = p.adjust(p_values$MAS5, method = "fdr")
);
# significant gene probes after correction
sig_gene_names <- data.table(
    "RMA" = gene_names[which(q_values$RMA < 0.01)],
    "MAS5" = gene_names[which(q_values$MAS5 < 0.01)]
);

# cel_samples <- sapply(
#     cel_files,
#     function(x) {
#         return(sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(x)));
#     },
#     USE.NAMES = FALSE
# );
# cel_datetimes <- sapply(cel_files, get_cel_datetime, USE.NAMES = FALSE);

# cel_table <- data.table(
#     "Sample" = cel_samples,
#     "Path" = cel_files,
#     "DateTime" = cel_datetimes
# );

