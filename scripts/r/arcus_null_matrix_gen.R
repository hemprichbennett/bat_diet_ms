##################################################
## Project: bat-diet (for full-species networks)
## Script purpose: Seeing the significance of various metrics across each of my networks SITEWISE
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


library(here)

# get command argument

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)



setwd(here())
getwd()

source('scripts/r/r_network_gen.r')
setwd(here())
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

for(i in 1:length(rawnets)){
  rawnets[[i]][2:nrow(rawnets[[i]]),] <- ifelse(rawnets[[i]][2:nrow(rawnets[[i]]), ] == 0, 0, 1)
}

netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))

names(netlists) <- names(rawnets)


nperm <- 1 # the number of random matrices to generate per dataset
# and null-model



margins_to_preserve <- c('both', 'columns')


  for(a in 1:length(names(netlists))){
    clust <- names(netlists)[a]
    #print(clust)
    for(s in 1:length(names(netlists[[a]]))){
      #print(s)
      site <- names(netlists[[a]][s])
      #print(site)
      
      for(mar in margins_to_preserve){
        #print(mar)
        
        if(mar == 'both'){
          out_mat <- vegan::permatswap(netlists[[a]][[s]],
                                       # the row/column/both sums to be retained. Can be 'none', 'rows',
                                       # 'columns', 'both'
                                       fixedmar = mar,
                                       mtype = 'count', #it's count data
                                       times = 1000, # number of permutations
                                       method = 'quasiswap' # using quasiswap for now
          )$perm[[1]]
        }else if(mar == 'columns'){
          # fix this
          out_mat <- vegan::permatswap(netlists[[a]][[s]],
                                       # the row/column/both sums to be retained. Can be 'none', 'rows',
                                       # 'columns', 'both'
                                       fixedmar = mar,
                                       mtype = 'count', #it's count data
                                       times = 1000, # number of permutations
                                       method = 'swsh' # using shuffle for now
          )$perm[[1]]
        }
        
        outpath <- paste0('data/output_data/null_matrices/', clust, '_', site, '_', mar, 
                          '_', args, '.csv')
        write.csv(out_mat, file = outpath)
      }
      
      
    }
  }

