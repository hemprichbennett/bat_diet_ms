#####Header ####
## Project: All bats
## Script purpose: Calculating how network metrics change when we adjust the number of bats in a network
## Date: 15/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: 

if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(ggridges)
  library(reshape2)
  library(forcats)
  library(netReducer)
  library(DataExplorer)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(netReducer, lib.loc = '/data/home/btw863/r_packages/')
  library(reshape2, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
  library(labeling, lib.loc = '/data/home/btw863/r_packages/')
  library(digest, lib.loc = '/data/home/btw863/r_packages/')
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
#rownames(all_interactions) <- all_interactions[,1]
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

names(sites_list)

sites_list <- sites_list[-which(names(sites_list)=='SBE')] #Theres no point doing this with SBE as it only has 1 species

#Make a vector of sample IDs and species
sampleIDs <- c(field_data$Species, field_data$Species)
names(sampleIDs) <- c(field_data$Faeces_no1, field_data$Faeces_no2)
which(colnames(sites_list$SAFE) %in% names(sampleIDs))

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

#####Data generation ####

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)

ind <- c('functional complementarity',
         'web asymmetry',
         'Alatalo interaction evenness',
         'togetherness',
         'Fisher alpha', 'mean number of shared partners',
         'niche overlap',
         'nestedness',
         'discrepancy',
         'ISA', 'weighted nestedness', 'NODF', 'weighted NODF', 'modularity')

chosen_ind <- ind[args]



#Now do the analysis
if(args !=14){
  out_df <- netreducing(input = sites_list, input_type = 'list', n_iterations = 100, min_nodes = 40, metric_chosen = chosen_ind,
                        type_chosen = 'network', level = 'higher')
}
if(args==14){
  out_df <- netreducing(input = sites_list, input_type = 'list', n_iterations = 100, min_nodes = 40, metric_chosen = 'modularity',
                       type_chosen = 'modularity', level = 'higher')
}


out_df$netnames <- gsub('DANUM', 'Danum', out_df$netnames)
out_df$netnames <- gsub('MALIAU', 'Maliau', out_df$netnames)

write.csv(out_df, paste('results/rarifying_networks/reducing_', chosen_ind,'_100.csv', sep =''))

out_df <- out_df[out_df$included==T,]

bigtax <- dcast(out_df[,c(1,5,6,7,8)], n_used + iteration +  netnames + metricval ~ Species, fun.aggregate = length, value.var = 'Species')

longtax <- melt(bigtax, id.vars = c('netnames', 'n_used', 'metricval'))
colnames(longtax)[5] <- 'Species'
bigtax$diversity <- sapply(seq(1,nrow(bigtax)), function(x) vegan::diversity(bigtax[x,seq(5, ncol(bigtax)),]))


palette <- c("#75aa56",
             "#8259b1",
             "#be7239")

n_bats_scatter <- ggplot(bigtax, aes(x = n_used, y = metricval, colour= netnames))+ 
  geom_point(alpha=0.8)+ scale_color_manual(values=palette, name = 'Site')+
  labs(x='Number of bats sampled', y= firstup(chosen_ind))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


n_bats_scatter  
pdf(paste('plots/netreducing/rarifying_', chosen_ind, 'n_bats.pdf', sep = ''))
n_bats_scatter
dev.off()

diversity_scatter <- ggplot(bigtax, aes(x = diversity, y = metricval, colour= netnames))+ 
  geom_point(alpha=0.8)+ scale_color_manual(values=palette, name = 'Site')+
  labs(x='Shannon diversity', y= firstup(chosen_ind))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


diversity_scatter  
pdf(paste('plots/netreducing/rarifying_', chosen_ind, 'n_bats.pdf', sep = ''))
diversity_scatter
dev.off()

sp_scatter <- ggplot(longtax, aes(x = metricval, y = metricval, colour= netnames))+ 
  geom_point(alpha=0.8)+ scale_color_manual(values=palette, name = 'Site')+
  labs(x='Number of individuals', y= firstup(chosen_ind))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~ Species)

sp_scatter

pdf(paste('plots/netreducing/rarifying_', chosen_ind, 'sp_n.pdf', sep = ''))
sp_scatter
dev.off()

