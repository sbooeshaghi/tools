#!/usr/bin/env Rscript
library(preseqR)
library(optparse)

option_list <- list(make_option(c("-v", "--value_counts"), type="string", default="./", help="A two-column matrix of values and value counts [default %default]"))

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser, positional_arguments = 1)
opt <- args$options
file <- args$args

if(opt$count_lines) {
      print(paste(length(readLines(file)) * opt$factor))
}
