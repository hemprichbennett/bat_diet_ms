rm(list=ls())
#To be ran in an array
if(interactive()==TRUE){
  library('here')
  library(here)
  library(bold)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(httpcode, lib.loc = '/data/home/btw863/r_packages/')
  library(urltools, lib.loc = '/data/home/btw863/r_packages/')
  library(bold, lib.loc = '/data/home/btw863/r_packages/')
  library(dplyr, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
  library(reshape2, lib.loc = '/data/home/btw863/r_packages/')
  library(data.table, lib.loc = '/data/home/btw863/r_packages/')
}
setwd(here())

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)

mydata <- read.table("data/processed_dna_data/galaxy_r_workflow/95/reps_95.tsv", header = F, sep = "\t", dec = ".", stringsAsFactors = FALSE)
colnames(mydata) <- c('seqID','seqs')

head(mydata)

mydata <- mydata[seq(nrow(mydata)/23*args, nrow(mydata)/23*(args+1)),] #Select the appropriate subset for the array



#make sure headers are not capitalized. We need to use this command to make a named list of the sequences, 
#otherwise our results from BOLD will not give the name of the OTU they correspond to
mydata2 <- as.list(setNames(mydata$seqs, mydata$seqID))  

t <- Sys.time()
output <- bold_identify(sequences = mydata2, db = "COX1_SPECIES", response=FALSE) #This can take several hours to run
t_after <- Sys.time()

cat('generating output took ', t_after - t, '\n')

outtax40 <- lapply(output, head, n=40)
outtax1 <- lapply(output, head, n = 1)
outtax1_df <- do.call('rbind', outtax1)



outtaxframe <- do.call("rbind", lapply(outtax40, data.frame))
write.csv(outtaxframe, paste('data/taxonomic_info', args,'.csv', sep = ''))

