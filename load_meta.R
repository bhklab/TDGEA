#source("https://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("affxparser");
#biocLite("affyio");

### Environment ###############################################################
library("affy");
library("affyio");
library("data.table");

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
    require(affyio);
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
	require(affyio);
	h <- affyio::read.celfile.header(filename, info="full");
	return(as.character(h$cdfName));
}


### Main ######################################################################
# Organize data for data table
cel_files <- list.files(
    pattern = "CEL",
    full.names = TRUE,
    recursive = TRUE
);
cel_samples <- sapply(
    cel_files,
    function(x) {
        return(sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(x)));
    },
    USE.NAMES = FALSE
);
cel_datetimes <- sapply(cel_files, get_cel_datetime, USE.NAMES = FALSE);
cel_means_rma <- rep(NA, length(cel_files));
cel_means_mas5 <- rep(NA, length(cel_files));
cel_medians_rma <- rep(NA, length(cel_files));
cel_medians_mas5 <- rep(NA, length(cel_files));

# cel_table <- data.table(
#     cel_samples,
#     cel_files,
#     cel_datetimes,
#     cel_means_rma,
#     cel_medians_rma
# );
# names(cel_table) <- c("Sample", "Path", "DateTime", "RMA_Mean", "RMA_Median");

# Histogram of dates
length_out <- 10;
delta <- (max(cel_datetimes, na.rm = TRUE) - min(cel_datetimes, na.rm = TRUE))/length_out;
ix <- seq(min(cel_datetimes, na.rm = TRUE) - delta, max(cel_datetimes, na.rm = TRUE), length.out = length_out + 1);
png("GPL570_dates_histogram.png", width = 3000, height = 1800, res = 300);
hist(
    cel_datetimes,
    xaxt = "n",
    xlab = "",
    main = ""
);
title("Dates of GPL570 chips");
axis(
    1,
    at = ix,
    labels = format(as.POSIXct(ix, origin="1970-01-01"), format= "%Y-%M-%d"),
    las = 2,
    cex.axis = 0.75
);
dev.off();


