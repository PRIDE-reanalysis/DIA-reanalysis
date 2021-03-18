#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-t", "--tric"), type="character", default=NULL, 
              help="The results file from the first OpenSwath analysis pipeline.", metavar="character"),
  make_option(c("-a", "--annot"), type="character", default=NULL, 
              help="annotation.txt file for MSstats", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="goldstandard.rda", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tric) || is.null(opt$annot) ){
  print_help(opt_parser)
  stop("Required input arguments must be supplied.", call.=FALSE)
}

#con <- file(opt$tric,"r")
#first_line <- readLines(con,n=1)
#close(con)
#print(first_line)
#con <- file(opt$annot,"r")
#first_line <- readLines(con,n=1)
#close(con)
#print(first_line)
#stop("optparse test stop.")


source("/scripts/OpenSWATHtoMSstatsFormat.R")
library(MSstats)


print ("Loading tric file...")
raw.openMS <- read.csv(opt$tric,
                       sep= "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE)

print ("Loading annotation file...")
annot_osw <- read.csv(opt$annot,
                      header = TRUE,
                      sep = "\t")

print ("Parsing annotation to MSstats format...")
input.openms <- OpenSWATHtoMSstatsFormat_new(raw.openMS,
                                         annotation = annot_osw,
                                         removeProtein_with1Feature=TRUE)

print("Run MSstats analysis...")
goldstandard.proposed <-dataProcess(input.openms,
                                    normalization=FALSE,
                                    summaryMethod="TMP",
                                    cutoffCensored="minFeature", 
                                    censoredInt="0",
                                    MBimpute=TRUE,
                                    maxQuantileforCensored=0.999)

print("Save the results")
save(goldstandard.proposed, file=opt$out)

