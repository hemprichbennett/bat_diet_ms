#library(tidyr)
library(stringr)
library(magrittr)
library(bipartite)
library(here)

setwd(here())

start_time <- lubridate::as_datetime(Sys.time())

# get command argument

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)


# set up variables for analysis
filepaths <- list.files('data/output_data/null_matrices/', full.names = T)

head(filepaths)

# make function to select the subset of matrices to be worked on in this
# iteration of the cluster job https://stackoverflow.com/a/16275428
splitter <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# make a vector of the paths of the files to be analysed in this iteration
to_analyse <- splitter(filepaths, 1000)[args] %>%
  unlist() # make it a vector, not a list

metrics <- c('weighted NODF', 'discrepancy', 'functional complementarity')
outlist <- list()
iterator <- 1



# in a loop, for each desired file 

for(i in 1:length(to_analyse)){
  # get the path of the file we're working on
  filepath <- as.character(to_analyse[i])
  filename <- gsub('.+\\/|\\..+', '', filepath)
  for(metric_used in metrics){
    # print(i)
    # print(metric_used)
    # print(filepath)
    # print(filename)
    clustering_level <- str_split(filename, pattern = '_')[[1]][1]
    site <- str_split(filename, pattern = '_')[[1]][2]
    fixedmargins <- str_split(filename, pattern = '_')[[1]][3]
    random_iteration <- str_split(filename, pattern = '_')[[1]][4]
    
    # do calc
    
    if(metric_used %in% c('discrepancy', 'functional complementarity')){
      metric_value <- read.csv(filepath) %>%
        networklevel(index = metric_used, level = 'higher')
    }else{
      metric_value <- read.csv(filepath) %>%
        networklevel(index = metric_used)
    }
    
    # print(metric_value)
    # save to list
    outlist[[iterator]] <- c(filepath, filename, metric_used, metric_value, 
                             clustering_level, site, fixedmargins, 
                             random_iteration)
    
    iterator <- iterator + 1
  }
}

out_df <- as.data.frame(do.call(rbind, outlist))
  
colnames(out_df) <- c('filepath', 'filename', 'metric_used', 'metric_value',
                      'clustering_level', 'site', 'fixedmargins', 'random_iteration')
  

write.csv(out_df, 
          file = paste0('data/output_data/null_values/', args, '.csv'),
          row.names = F, col.names = T)

end_time <- lubridate::as_datetime(Sys.time())
time_elapsed <- end_time - start_time
print(time_elapsed)