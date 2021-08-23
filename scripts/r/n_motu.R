##################################################
## Project: bat-diet (for full-species networks)
## Script purpose: Seeing the significance of various metrics across each of my networks SITEWISE
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


if(interactive()==TRUE){
  library(here)
}else{
  library(reshape2, lib.loc = '/data/home/btw863/r_packages/')
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(LOTUS, lib.loc = '/data/home/btw863/r_packages/')
}



setwd(here())

source('scripts/r/r_network_gen.r')
# getwd()

inpath <- 'data/processed_dna_data/lulu/'
filenames <- list.files(pattern = '.csv', path = inpath)
#filenames <- filenames
filenames
filenames <- paste(inpath, filenames, sep = '')
#filenames <- filenames[grep('lulu', filenames)]

rawnets <- lapply(filenames, read.csv, header = F, stringsAsFactors = F, row.names=1)
names(rawnets) <- gsub('.*\\/', '', filenames)
names(rawnets) <- gsub('_.+', '', names(rawnets))
netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))

names(netlists) <- names(rawnets)



n_MOTU <- sapply(netlists, function(a) sapply(a, function(x) nrow(x)))
n_MOTU <- t(n_MOTU)


write.csv(n_MOTU, 'results/n_MOTU.csv')
