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

# query_dbs
# Description:
#   Query GEO Datasets
# Inputs:
#   query_params:    Array of string queries
# Outputs:
#   dbs:    List of series resulting from the query
query_dbs <- function(query_params) {
    # initial data query URL
    gds_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=";
    gds_params <- "&retmax=5000&usehistory=y&retmode=json";
    # detailed information URL from previous query
    gds_url <- URLencode(paste0(
        gds_url,
        do.call("paste", as.list(c(query_params, sep = "+AND+"))),
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
}

# download_dbs
# Description:
#   Download data series from GEO Datasets
# Inputs:
#   series:         Array of series identifiers
#   download_dir:   Directory to download files to
download_dbs <- function(series, download_dir = getwd()) {
    gds_download_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series";

    # downloading proper datasets
    for (GSE in series) {
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
            destfile = paste0(download_dir, GSE, "_RAW.tar")
        );
    }
    
}

# rescale
# Description:
#   Linearly rescale array of values to within a new range
# Inputs:
#   x:              Values to be rescaled
#   lower:          Minimum value of new scale
#   upper:          Maximum value of new scale
# Outputs:
#   scaled_values:  Rescaled values
rescale <- function(x, lower = 0, upper = 1) {
    if (max(x) == min(x)) {
        stop("Cannot be rescaled, max and min values are equal.");
    }
    m <- (upper - lower)/(max(x) - min(x));
    return(m*(x-min(x)) + lower);
}


### Main ######################################################################
### Query datasets ##############################
# search parameters (these can be changed for other queries)
dbs <- c(
    "GSE9891",
    "GSE18520",
    "GSE26193",
    "GSE44104"
); 


### Find CEL files and DateTimes ################
print("Loading CEL files");
# list of CEL files
cel_files <- sapply(
    dbs,
    function (GSE) {
        return(list.files(
            path = paste0("../data/raw/", GSE, "_RAW"),
            pattern = "CEL",
            full.names = TRUE,
            recursive = TRUE
        ));
    },
    simplify = FALSE
);

# recall DateTimes based on sample names
# can do this since all dbs have the same gene probes (all on the same platform)
gene_names <- geneNames(affy_dbs[[dbs[1]]]);
# get datetimes for CEL files
cel_datetimes <- sapply(
    dbs,
    function (GSE) {
        return(sapply(cel_files[[GSE]], get_cel_datetime));
    },
    simplify = FALSE
);


### Filter microarrays with no DateTime info ####
# get CEL files with no times
no_times <- sapply(
    dbs,
    function (GSE) {
        return(-which(is.na(cel_datetimes[[GSE]])));
    },
    simplify = FALSE
);

# remove from list of CEL files
cel_files <- sapply(
    dbs,
    function (GSE) {
        return(cel_files[[GSE]][no_times[[GSE]]]);
    },
    simplify = FALSE
);

# remove from DateTimes
cel_datetimes <- sapply(
    dbs,
    function (GSE) {
        return(cel_datetimes[[GSE]][no_times[[GSE]]]);
    },
    simplify = FALSE
);

# scale DateTimes to be in [0,1]
scaled_datetimes <- sapply (
    dbs,
    function (GSE) {
        return(rescale(cel_datetimes[[GSE]]));
    },
    simplify = FALSE
);

### Normalize data ##############################
# RMA-normalized data from CEL files
affy_dbs_rma <- sapply(
    dbs,
    function (GSE) {
        print(GSE);
        return(affy::justRMA(filenames = cel_files[[GSE]]));
    },
    simplify = FALSE
);
