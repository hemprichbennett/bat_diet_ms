## @knitr species_removal
#### Header ####
## Project: bat-diet
## Script purpose: Finding the effect of removing different bat species from out networks
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(ggplot2)
library(magrittr)
library(forcats)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots



source('scripts/r/r_network_gen.r')

nets <- r_network_gen(lulu=T, filter_species = T)

# ind <- c('functional complementarity',
#          'web asymmetry',
#          'Alatalo interaction evenness',
#          'togetherness',
#          'Fisher alpha', 'mean number of shared partners',
#          'niche overlap',
#          'nestedness',
#          'discrepancy',
#          'ISA', 'weighted nestedness', 'NODF', 'weighted NODF')

ind <- c('functional complementarity',
         'mean number of shared partners',
         'niche overlap',
         'discrepancy',
        'NODF')


orig <- lapply(ind, function(i) lapply(nets, function(x) bipartite::networklevel(x, index = i, level = 'higher')))


names(orig) <- ind




mods <- lapply(nets, function(x) slot(bipartite::computeModules(web = x), 'likelihood'))

mods <- melt(mods)
colnames(mods) <- c('value', 'Network')
mods$Metric <- 'modularity'
ind <- c(ind, 'modularity')



orig_melted <- melt(orig)
colnames(orig_melted) <- c('value', 'Network', 'Metric')

orig_melted <- rbind(orig_melted, mods)
#orig_melted$minus_species <- "No species removed"


sp_deleter <- function(networks, chosen_index){
  outlist <- list()
  for(i in 1:length(networks)){
    if(chosen_index=='modularity'){
      a <- sapply(seq(1,ncol(networks[[i]])), function(a) slot(bipartite::computeModules(web =networks[[i]][,-a]), 'likelihood'))
    }else{
      a <- sapply(seq(1,ncol(networks[[i]])), function(a) bipartite::networklevel(networks[[i]][,-a], index = chosen_index, level = 'higher'))
    }
    
    a <- data.frame('value'= a, 'minus_species' =colnames(networks[[i]]))
    #print(a)
    outlist[[names(networks)[i]]] <- a
  }
  #cols <- lapply(networks, function(x) colnames(x))
  return(outlist)
}


#slot(bipartite::computeModules(web = nets$SAFE), 'likelihood')

#one_met <- sp_deleter(nets, chosen_index = 'nestedness')
#one_melted <- melt(one_met)

all_mets <- lapply(ind, function(i) sp_deleter(nets, chosen_index = i))
names(all_mets) <- ind
melted_all <- melt(all_mets)
melted_all <- melted_all[,-2]
colnames(melted_all) <- c('minus_species', 'value', 'Network', 'Metric')

#melted_all$Network <- gsub('DANUM', 'Danum', melted_all$Network)
#melted_all$Network <- gsub('MALIAU', 'Maliau', melted_all$Network)

#plot_str(all_mets)


master_df <- rbind(orig_melted, melted_all)
master_df$Network <- gsub('DANUM', 'Danum', master_df$Network)
master_df$Network <- gsub('MALIAU', 'Maliau', master_df$Network)

master_df$minus_species <- as.factor(master_df$minus_species)


master_df$Metric <- firstup(master_df$Metric)

#Make a second combined df, by merging the values

orig_melted$net_and_met <- paste(orig_melted$Network, orig_melted$Metric, sep = '_')
colnames(orig_melted)[1] <- 'actual_value'
melted_all$network_and_met <- paste(melted_all$Network, melted_all$Metric, sep = '_')

merged <- merge(x= orig_melted[c(1,4)], y = melted_all, by.x = 'net_and_met', by.y = 'network_and_met')

# palette <- c("#75aa56",
#              "#8259b1",
#              "#be7239")

merged$minus_species %<>% 
  gsub('Hice', 'Hipposideros cervinus', .)%<>%
  gsub('Hidi', 'Hipposideros diadema', .)%<>%
  gsub('Hidy', 'Hipposideros dyacorum', .)%<>%
  gsub('Hiri', 'Hipposideros ridleyi', .)%<>%
  gsub('Keha', 'Kerivoula hardwickii', .)%<>%
  gsub('Kein', 'Kerivoula intermedia', .)%<>%
  gsub('Kepa', 'Kerivoula papillosa', .)%<>%
  gsub('Kepe', 'Kerivoula pellucida', .)%<>%
  gsub('Rhbo', 'Rhinolophus borneensis', .)%<>%
  gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
  gsub('Rhtr', 'Rhinolophus trifoliatus', .)

merged$Network <- gsub('DANUM', 'Danum', merged$Network)
merged$Network <- gsub('MALIAU', 'Maliau', merged$Network)


merged$Metric <- firstup(merged$Metric)


merged$Metric <- gsub('Mean number of shared partners', 'Mean number of\nshared partners', merged$Metric)
merged$Metric <- gsub('Functional complementarity', 'Functional\ncomplementarity', merged$Metric)
merged$minus_species <- as.factor(merged$minus_species)

#Calculate the effect the species is having
merged$impact <- abs(merged$actual_value - merged$value)

#Now rank the impacts
merged$rankings <- NA

for(i in 1: length(unique(merged$net_and_met))){
  n_m <- unique(merged$net_and_met)[i]
  rows <- which(merged$net_and_met== n_m)
  rows_to_write <- rows[order(merged[rows,'impact'], decreasing = T)]
  merged$rankings[rows_to_write] <- seq(1:length(rows_to_write))
}



most_inf <- ggplot(merged[merged$rankings<=5,], aes(x= Network, y = fct_rev(minus_species)))+ geom_tile(aes(fill=rankings))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~ Metric)+
  scale_fill_gradient(low = 'black', high = 'lightblue', name = 'Influence\nranking')+
  labs(x = 'Site', y = 'Species')
most_inf

pdf('plots/species/species_influence_top.pdf')
most_inf
dev.off()

italic.text <- element_text(face = "italic")

all_sp <- ggplot(merged, aes(x= Network, y = fct_rev(minus_species)))+ geom_tile(aes(fill=rankings))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~ Metric)+
  scale_fill_gradient(low = 'black', high = 'lightblue', name = 'Influence\nranking')+
  labs(x = 'Site', y = 'Species')+ theme(axis.text.y = italic.text, legend.position="bottom")
all_sp  

pdf('plots/species/species_influence_all.pdf', height = 10)
all_sp
dev.off()
