# ---- this is an updated method for creating unique gene names for the features table from the CellRanger output
#incorporates argparse for convenient execution
library(argparser)

p <- arg_parser("Iterate through sequencing samples from CellRanger")

p <- add_argument(p, arg = "--folder", help = "input file directory")

p <- add_argument(p, arg = "--samples", help = "input samples ID from input directory, separate by commas")

argv <- parse_args(p)

print(paste0("The currrent directory is: ", argv$folder))

print(paste0("The samples in the input directories are: ", argv$samples))

data_dir <- argv$folder

sample.ls <- unlist(strsplit(argv$samples,split = ","))

print(paste0(sample.ls))

#loop thru each
for (file in 1:length(sample.ls)){
  setwd(paste0(data_dir, sample.ls[file], "/outs/filtered_feature_bc_matrix/"))
  print(getwd())
  #Start R in directory with features.tsv.gz
  if(file.exists("./features.tsv.gz")) {
    system("mkdir backup zipped unzipped ; cp -n *.gz backup/ ; gunzip *.gz")
  }
  features <- read.delim("./features.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
  # Check for repeated genes
  grep(" \\(1 of many\\)", features[,2], value = TRUE)
  grep(" ", features[,2], value = TRUE)
  # Remove "(1 of many)" from common gene names
  features[,2] <- gsub(" \\(1 of many\\)", "", features[,2])
  # find number of duplicates
  table(duplicated(features[,2]))
  # make repeats unique
  features[,2] <- make.unique(features[,2], sep = "*")
  if("Gene Expression" %in% features[,3]) {
    features <- features[,-3]
  }
  write.table(features, "./features.tsv",
              quote = FALSE, sep = "\t", row.names = FALSE)
  system("cp -n * unzipped/ ; gzip * ; mv *gz zipped")
  
}

#output will save the original filtered mtx in "backups"
#will create the features.tsv file with uniq gene names in "unzipped" and "zipped" files