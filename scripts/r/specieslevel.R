#### Header ####
## Project: Bat-diet
## Script purpose: Calculating species-level metrics for each network
## Date: 25/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(here)
library(ggplot2)
library(magrittr)
library(reshape2)
library(forcats)
library(tidyverse)

source(here('scripts', 'r', 'r_network_gen.r'))

nets <- r_network_gen(lulu = T, filter_species = T)


#####Do the actual calculations ####
sp_mets <- lapply(nets, function(x) specieslevel(x, level = 'higher'))


####Format the data ####
sp_df <- do.call(rbind, sp_mets)
sp_df$site <- gsub('\\..+', '', rownames(sp_df))
sp_df$species <- gsub('.+\\.', '', rownames(sp_df))

sp_df$species %<>%
  gsub('Hice', 'Hipposideros cervinus', .)%<>%
  gsub('Hidi', 'Hipposideros diadema', .)%<>%
  gsub('Hidy', 'Hipposideros dyacorum', .)%<>%
  gsub('Hiri', 'Hipposideros ridleyi', .)%<>%
  gsub('Keha', 'Kerivoula hardwickii', .)%<>%
  gsub('Kein', 'Kerivoula intermedia', .)%<>%
  gsub('Kemi', 'Kerivoula minuta', .)%<>%
  gsub('Kepa', 'Kerivoula papillosa', .)%<>%
  gsub('Rhbo', 'Rhinolophus borneensis', .)%<>%
  gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
  gsub('Rhtr', 'Rhinolophus trifoliatus', .)




write.csv(sp_df, 'results/species_level_data.csv')



sp_df <- read.csv(here('results', 'species_level_data.csv'), row.names = 1)

sp_df$site %<>%
  gsub('DANUM', 'Danum', .)%<>%
  gsub('MALIAU', 'Maliau', .)


degree_summaries <- sp_df %>%
  group_by(site) %>%
  summarise(mean_cc = mean(closeness),
            mean_bc = mean(betweenness),
            sd_cc = sd(closeness),
            sd_bc = sd(betweenness))

degree_summaries

#sp_df
ggplot(sp_df, aes(x=species, y = betweenness))+ geom_point()+ facet_wrap(~ site)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_melted <- melt(sp_df, id.vars = c('species', 'site'))

all_melted$variable <- gsub('\\.', ' ', all_melted$variable)

ggplot(all_melted, aes(x=site, y = value, colour = species))+ geom_point()+ facet_wrap(~ variable, scales = 'free')+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



#Make some for centrality, degree and species strength
centrality_melted <- all_melted[which(all_melted$variable %in% c('betweenness', 
                                                                 'closeness', 'degree',
                                                                 'species strength')),]

centrality_melted$variable %<>%
  gsub('betweenness', 'Betweenness centrality', .)%<>%
  gsub('closeness', 'Closeness centrality', .)%<>%
  gsub('degree', 'Degree', .) %<>%
  gsub('species strength', 'Species strength', .)

centrality_plot_list <- list()
for(i in 1:length(unique(centrality_melted$variable))){
  met <- unique(centrality_melted$variable)[i]
  if(met =="Closeness centrality" ){
    centrality_plot_list[[met]] <- ggplot(centrality_melted[which(centrality_melted$variable==met),], aes(x=site, y = fct_rev(species)))+ 
      geom_tile(aes(fill=value))+
      scale_fill_gradient(low = "white",
                          high = "black", name = gsub(' ', '\n', met), breaks = c(0.08, 0.1, 0.12)) + 
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.text.y = element_text(face = "italic"), legend.position="bottom")+
      labs(x= 'Site', y = 'Species')
    pdf(paste('plots/species/', met, '.pdf', sep = ''), height = 9)
    print(centrality_plot_list[[met]])
    dev.off()
  }
  if(met == "Species strength"){
    centrality_plot_list[[met]] <- ggplot(centrality_melted[which(centrality_melted$variable==met),], aes(x=site, y = fct_rev(species)))+ 
      geom_tile(aes(fill=value))+
      scale_fill_gradient2(low = 'white', mid = 'lightblue', high = 'black',
                          name = gsub(' ', '\n', met), breaks =c(350, 700, 1050, 1400)) + 
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.text.y = element_text(face = "italic"), legend.position="bottom",
            axis.line = element_line(colour = "black"),
            panel.background = element_rect(fill = "darkgray",
                                            colour = "darkgray"))+
      labs(x= 'Site', y = 'Species')
    pdf(paste('plots/species/', met, '.pdf', sep = ''), height = 9)
    print(centrality_plot_list[[met]])
    dev.off()
  }
  else{
    centrality_plot_list[[met]] <- ggplot(centrality_melted[which(centrality_melted$variable==met),], aes(x=site, y = fct_rev(species)))+ 
      geom_tile(aes(fill=value))+
      scale_fill_gradient(low = "white",
                          high = "black", name = gsub(' ', '\n', met)) + 
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.text.y = element_text(face = "italic"), legend.position="bottom")+
      labs(x= 'Site', y = 'Species')
    pdf(paste('plots/species/', met, '.pdf', sep = ''), height = 9)
    print(centrality_plot_list[[met]])
    dev.off()
  }
  
}


centrality_plot_list$betweenness



# ggplot(centrality_melted, aes(x=site, y = species))+ geom_tile(aes(fill=value))+
#   facet_wrap(~ variable, scales = 'free_y')+
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))






just_hice <- sp_df[which(sp_df$species=='Hipposideros cervinus'),]

just_hice$species <- NULL
melted_hice <- melt(just_hice, id.vars = c('site'))  

ggplot(melted_hice, aes(x=site, y = value))+ geom_point()+ facet_wrap(~ variable, scales = 'free')+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



####Now calculate the degree distribution of the networks####

#Make individual-based networks

ind_nets <- r_network_gen(collapse_species = F, lulu = T, filter_species = T)
ind_nets <- ind_nets[,-1] #These are the rownames but we dont need them

ind_nets <- ind_nets[,-which(!ind_nets[1,] %in% c('SAFE', 'DANUM', 'MALIAU'))]

ind_netlist <- list()

for(i in 1:length(unique(as.character(ind_nets[1,])))){
  site <- unique(as.character(ind_nets[1,]))[i]
  print(site)
  mat <- ind_nets[3:nrow(ind_nets),which(ind_nets[1,]==site)]
  ind_netlist[[site]] <- apply(mat, 2, as.numeric)
}

degree_dist <- melt(lapply(ind_netlist, function(x) colSums(x)))
colnames(degree_dist) <- c('Degree', 'Site')

ggplot(degree_dist, aes(x=Degree))+ geom_histogram(binwidth = 1)+ facet_wrap(~Site)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

fit <- aov(Degree ~ Site, degree_dist)

plot(fit)

TukeyHSD(fit)
  
