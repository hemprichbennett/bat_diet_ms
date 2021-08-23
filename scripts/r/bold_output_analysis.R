##################################################
## Project: all bats
## Script purpose: analysing the output of the array job 'bold_querying.R'
## Date: 06/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(here)
library(bold)
library(ggmap)
library(magrittr)


setwd(here())

files <- list.files(pattern = 'taxonomic_info*', path = 'data/')
files <- paste('data/', files, sep = '')
#One of these files is an old attempt at running a non-array job, one was the final one and returned empty
files <- files[-c(1,17)]
 
files <- lapply(files, read.csv)

alltax <- do.call(rbind, files)

pests <- scan('data/pest_species.txt', character(), sep = ',')
#remove empty values 
pests <- pests[-which(pests=='')]

#These are potential pest MOTU consumed by bats
pestconsumption <- alltax[which(alltax$taxonomicidentification %in% pests),]


write.csv(pestconsumption, 'results/pestvals.csv')

#####Plot where the pests were collected####
midpoint <- c(lon=	117.321995, lat = 4.699939)
se_asia <- get_map(midpoint, source = 'google', maptype = 'roadmap', zoom = 4)

pests <- ggmap(se_asia) +
  geom_point(aes(x = specimen_lon, y = specimen_lat), col = 'orange', alpha = 0.2,
             data = pestconsumption)

pests


#####Now plot all matches####
all_midpoint <- c(lon=	18, lat = 4.699939)
world <- get_map(all_midpoint, source = 'google', maptype = 'hybrid', zoom = 1)

map <- ggmap(world) +
  geom_point(aes(x = specimen_lon, y = specimen_lat), col = 'orange', alpha = 0.2,
             data = alltax)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"), axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  labs(x = NULL, y= NULL)
map  
pdf('plots/worldwide_bold_matches.pdf')
map
dev.off()




#Extract the bestmatches
best <- as.character(pestconsumption[which(pestconsumption$similarity > 0.97),'X'])
names(best) <- as.character(pestconsumption[which(pestconsumption$similarity > 0.97),'taxonomicidentification'])
best <- gsub('\\..+', '', best)



source('scripts/r/r_network_gen.r')


nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T, include_malua = F)
names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))



for(i in 1:length(nets)){
  rows <- which(rownames(nets[[i]]) %in% best)
  otus <- rownames(nets[[i]][which(rownames(nets[[i]]) %in% best),])
  subnet <- nets[[i]][rows,]
  #consumption <- 
}

for(i in 2:2){
#for(i in 1:lenght(nets)){
  otus <- best[best %in% rownames(nets[[i]])] # Find the otu matches for this network
  preysp <- names(best[best %in% rownames(nets[[i]])]) #get their names
  #make a subnetwork containing only the relevant otus and bats
  subnet_1 <- nets[[i]][which(rownames(nets[[i]]) %in% otus),]
  subnet_2 <- subnet_1[,which(colSums(subnet_1)>0)]
  melt(subnet_2)
}


best_df <- data.frame(best, names(best))
colnames(best_df)[2] <- 'Pest species'


#Format the nets to merge them with the pest data
edgelist <- melt(nets)
colnames(edgelist) <- c('otu', 'Bat species', 'Number of bats consuming', 'Site')
edgelist <- edgelist[edgelist$`Number of bats consuming`>0,]
edgelist <- edgelist[edgelist$otu %in% best_df$best,]



consumptiondf <- merge(y = best_df, x = edgelist, by.y = 'best', by.x= 'otu', all.x = T, all.y = F)

consumptiondf$`Bat species` %<>% 
  gsub('Hice', 'Hipposideros cervinus', .)%<>%
  gsub('Kein', 'Kerivoula intermedia', .)%<>%
  gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
  gsub('Rhtr', 'Rhinolophus trifoliatus', .)%<>%
  gsub('Hidi', 'Hipposideros diadema', .)#%<>%

consumptiondf$otu <- NULL

consumptiondf <- consumptiondf[,c(4,1,3,2)]
consumptiondf <- consumptiondf[
  with(consumptiondf, order(Site, `Bat species`)),
  ]
write.csv(consumptiondf, 'results/pest_consumption.csv', row.names = F)





