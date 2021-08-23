#### Header ####
## Project: bat-diet
## Script purpose: using the output files from the netreducing array job to make plots
## Date: 27/07/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(MASS)
library(reshape2)
library(ggplot2)
library(here)

setwd(here())

allfiles <- list.files(path = 'results/rarifying_networks/', pattern = '.csv')
allfiles <- paste('results/rarifying_networks/', allfiles, sep = '')

#Modularity was calculated separately in a big array, and its outputs aren't correctly formatted for this script
allfiles <- allfiles[-grep('modularity', allfiles)]
# some were made with a subset, and we don't care about them
allfiles <- allfiles[-grep('smaller', allfiles)]

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots


for(i in 1: length(allfiles)){
  infilename <- allfiles[i]
  
  cat('reading', infilename, '\n')
  out_df <- read.csv(infilename)
  
  cat(infilename, 'loaded\n')
  infilename <- gsub('.+/reducing_', '', infilename)
  infilename <- gsub('_100.csv', '', infilename)
  
  out_df <- out_df[out_df$included ==T,]
  
  
  
  # categories <- unique(out_df[,c('Species','netnames','iteration', 'n_used')])
  # 
  # sp <- categories$Species[1]
  # net <- categories$netnames[1]
  # it <- categories$iteration[1]
  # n_used <- categories$n_used[1]
  # 
  # n_bats <- length(which(out_df$Species == sp & out_df$netnames == net & out_df$iteration== it & out_df$n_used == n_used))
  # 
  # 
  # which(out_df$Species == categories$Species[2] && out_df$netnames == categories$netnames[2] && out_df$iteration== categories$iteration[2] && out_df$n_used == categories$n_used[2])
  # 
  
  
  
  #possible_list <- table(out_df$Species , out_df$netnames, out_df$n_used, out_df$iteration)
  
  bigtax <- dcast(out_df[,c(2,5,6,7,8,9)], n_used + netnames + metricval ~ Species, fun.aggregate = length)
  bigtax$diversity <- sapply(seq(1,nrow(bigtax)), function(x) vegan::diversity(bigtax[x,seq(5, ncol(bigtax)),]))
  
  longtax <- melt(bigtax, id.vars = c('netnames', 'n_used', 'metricval', 'diversity'))
  colnames(longtax)[5] <- 'Species'
  
  
  # palette <- c("#75aa56",
  #              "#8259b1",
  #              "#be7239")
  
  
  diversity_scatter <- ggplot(bigtax, aes(x = diversity, y = metricval, colour= netnames))+ 
    geom_point(alpha=0.3)+ 
    scale_colour_viridis_d(name = 'Site')+
    #scale_color_manual(values=palette, name = 'Site')+
    labs(x='Shannon diversity', y = firstup(infilename))+
    theme_bw() + theme(legend.position="bottom",
                       panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       text = element_text(size = 20))
  
  
  diversity_scatter  
  pdf(paste('plots/netreducing/rarifying_', infilename, '_diversity.pdf', sep = ''), width = 5)
  print(diversity_scatter)
  dev.off()
  
  sp_scatter <- ggplot(longtax, aes(x = value, y = metricval, colour= netnames))+ 
    geom_point(alpha=0.3)+ 
    scale_colour_viridis_d(name = 'Site')+
    #scale_color_manual(values=palette, name = 'Site')+
    labs(x='Number of individuals', y = infilename)+
    theme_bw() + theme(legend.position="bottom",
                       panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       text = element_text(size = 20))+
    facet_wrap(~ Species, scales = 'free_x')
  
  sp_scatter  
  pdf(paste('plots/netreducing/rarifying_', infilename, '_sp.pdf', sep = ''), width = 5)
  print(sp_scatter)
  dev.off()
  
  
  
  individuals_scatter <- ggplot(longtax, aes(x = n_used, y = metricval, colour= netnames))+ 
    geom_point(alpha=0.3)+ 
    scale_colour_viridis_d(name = 'Site')+
    #scale_color_manual(values=palette, name = 'Site')+
    labs(x='Number of individuals', y = firstup(infilename))+
    theme_bw() + theme(legend.position="bottom",
                       panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       text = element_text(size = 20))
  individuals_scatter
  pdf(paste('plots/netreducing/rarifying_', infilename, '_individuals.pdf', sep = ''), width = 5)
  print(individuals_scatter)
  dev.off()
  
  jpeg(paste('plots/netreducing/rarifying_', infilename, '_individuals.jpeg', sep = ''), width = 1800, height =2400, res = 300)
  print(individuals_scatter)
  dev.off()
  
  tiff(paste('plots/netreducing/rarifying_', infilename, '_individuals.tiff', sep = ''), width = 1800, height =2400, res = 500)
  print(individuals_scatter)
  dev.off()
  
  tiff(paste('plots/netreducing/rarifying_', infilename, '_diversity.tiff', sep = ''), width = 1800, height =2400, res = 500)
  print(diversity_scatter)
  dev.off()
  
  jpeg(paste('plots/netreducing/rarifying_', infilename, '_diversity.jpeg', sep = ''), width = 1800, height =2400, res = 300)
  print(diversity_scatter)
  dev.off()
  
  cat(infilename, ' finished\n')
}

