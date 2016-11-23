#source("https://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("affxparser");

library("affxparser");
files <- c("GSE20711_RAW/GSM519722.CEL");
header_meta <- strsplit(readCelHeader(files[[1]])$datheader, "\\s+")[[1]];
date <- header_meta[[8]];
time <- header_meta[[9]];

date;
time;