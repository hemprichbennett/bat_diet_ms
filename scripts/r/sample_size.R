## @knitr sample_count_setup
#####Header ####
## Project: All bats
## Script purpose: Calculating the number of samples per species we have for each network
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: 


dir <- getwd()
basedir <- strsplit(dir, split ='/')[[1]][2]
#print(basedir)
if(grepl('data', basedir)){
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  
}else{
  library('here')
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(reshape2)
  library(corrplot)
  library(iNEXT)
  library(DataExplorer)
  library(knitr)
  library(kableExtra)
}
setwd(here())
source('scripts/r/r_network_gen.r')

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots


#####Format the data ####


field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'), stringsAsFactors = F)
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data$Site <- gsub('DVCA', 'Danum', field_data$Site)
field_data$Site <- gsub('DANUM', 'Danum', field_data$Site)
field_data$Site <- gsub('MALIAU', 'Maliau', field_data$Site)

all_interactions <- r_network_gen(collapse_species = F, desired_species = NULL, include_malua = T, filter_species = T, lulu = T)

desired_cols <- c('MOTU', 'DANUM', 'MALIAU', 'SAFE', 'SBE')

all_interactions <- all_interactions[,-which(!all_interactions[1,] %in% desired_cols)]


colnames(all_interactions) <- all_interactions[2,]
all_interactions <- all_interactions[-c(2),]
rownames(all_interactions) <- all_interactions[,1]
all_interactions <- all_interactions[,-1]

locations <- c()

for(i in 1:ncol(all_interactions)){
  locations[i] <- all_interactions[1,i]
  names(locations)[i] <- all_interactions[1,i]
}


row_names <- rownames(all_interactions) #Store these as an object as the apply below kills them
all_interactions <- apply(all_interactions, 2, as.numeric)
rownames(all_interactions) <- row_names


#####Make a list with a network for each site####

sites_list <- list()


for(i in 1:length(unique(names(locations)))){
  loc = unique(names(locations))[i]
  sites_list[[i]] <- all_interactions[,which(names(locations)==loc)]
  sites_list[[i]] <- sites_list[[i]][-1,]#Now that we have the site, this row is worthless
  
}
names(sites_list) <- unique(names(locations))

#names(sites_list)



#Make a vector of sample IDs and species
sampleIDs <- c(field_data$Species, field_data$Species)
names(sampleIDs) <- c(field_data$Faeces_no1, field_data$Faeces_no2)
#which(colnames(sites_list$SAFE) %in% names(sampleIDs))

#Rename all the columns to species-site
for(n in 1:length(sites_list)){
  badcols <- c()
  for(i in 1:ncol(sites_list[[n]])){
    pos <- which(names(sampleIDs)==colnames(sites_list[[n]])[i])
    pos <- pos[1]
    #print(colnames(sites_list[[n]])[i])
    if(is.na(pos)){
      badcols <- c(badcols, i)
      next()
    }
    #print(paste(sampleIDs[pos], names(sampleIDs[pos]), sep = '-'))
    colnames(sites_list[[n]])[i] <- paste(sampleIDs[pos], names(sampleIDs[pos]), sep = '-')
  }
  if(length(badcols)>0){
    sites_list[[n]] <- sites_list[[n]][,-badcols]
  }
}

names(sites_list) <- gsub('DANUM', 'Danum', names(sites_list))
names(sites_list) <- gsub('MALIAU', 'Maliau', names(sites_list))

sp_long <- melt(lapply(sites_list, function(x) table(gsub('-.+', '', colnames(x)))))
sp_wide <- cast(sp_long, formula = Var.1 ~ L1)
rownames(sp_wide) <- sp_wide[,1]
sp_wide <- sp_wide[,-1]

#Remove NAs
sp_wide[is.na(sp_wide)] <- 0

#Add shannon diversity
sp_wide <- rbind(sp_wide, apply( sp_wide, 2, function(x) vegan::diversity(x)))
rownames(sp_wide)[nrow(sp_wide)] <- 'Shannon diversity'



sp_wide %>%
  kable() %>%
  kable_styling()

## @knitr saving_sample_counts

write.csv(sp_wide, 'results/sample_counts.csv')