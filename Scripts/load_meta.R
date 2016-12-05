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
dbs <- c("GSE19615", "GSE18864", "GSE3744", "GSE12276", "GSE20711");
cel_files <- list();
affy_dbs <- list();
affy_dbs_rma <- list();
affy_dbs_mas5 <- list();

# Load, process, and normalize CEL files based on GSE
for (GSE in dbs) {
    cel_files[[GSE]] <- list.files(
        path = paste("../CEL/", GSE, "_RAW", sep = ""),
        pattern = "CEL",
        full.names = TRUE,
        recursive = TRUE
    );

    # Normalize by GSE
    affy_dbs[[GSE]] <- affy::ReadAffy(filenames = cel_files[[GSE]]);    # read in to AffyBatch object
    affy_dbs_rma[[GSE]] <- affy::rma(affy_dbs[[GSE]]);              # RMA normalize
    affy_dbs_mas5[[GSE]] <- affy::mas5(affy_dbs[[GSE]]);            # MAS5 normalize
}

# build structure for recalling datetimes based on sample names
# these calls will be made often below, so this is saving processing
gene_names <- geneNames(affy_dbs[[dbs[1]]]); # can do this since all dbs have the same genes
cel_datetimes <- list();
# get datetimes for CEL files
for (GSE in dbs) {
    cel_datetimes[[GSE]] <- sapply(cel_files[[GSE]], get_cel_datetime);
}

# p-values of linear regression having no slope
p_values_rma <- list();
p_values_mas5 <- list();
q_values_rma <- list();
q_values_mas5 <- list();
for (GSE in dbs) {
    print(GSE);
    # apply regression test to each row of matrix containing expression for a probe at the time it was taken
    p_values_rma[[GSE]] <- apply(
        exprs(affy_dbs_rma[[GSE]]),
        1,
        function (a) {
            probe_lm <- lm(a ~ cel_datetimes[[GSE]]);
            gene_summary <- summary.lm(probe_lm);
            return(pf(
                gene_summary$fstatistic[1],
                gene_summary$fstatistic[2],
                gene_summary$fstatistic[3],
                lower.tail = FALSE)
            );
        }
    );
    p_values_mas5[[GSE]] <- apply(
        exprs(affy_dbs_mas5[[GSE]]),
        1,
        function (a) {
            probe_lm <- lm(a ~ cel_datetimes[[GSE]]);
            gene_summary <- summary.lm(probe_lm);
            return(pf(
                gene_summary$fstatistic[1],
                gene_summary$fstatistic[2],
                gene_summary$fstatistic[3],
                lower.tail = FALSE)
            );
        }
    );

    # q-values after FDR correction
    q_values_rma[[GSE]] <- p.adjust(p_values_rma[[GSE]], method = "fdr");
    q_values_mas5[[GSE]] <- p.adjust(p_values_mas5[[GSE]], method = "fdr");
}

# mapping gene probes to other formats
hgu_probe_map <- select(hgu133plus2.db, gene_names, c("SYMBOL","ENTREZID"));

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

