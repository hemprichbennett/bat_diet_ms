  ## @knitr otu_setup
  
  #setwd()
  
  #Make a nested list, where the main list has a value for each clustering level, and each item within it is a distinct network generated at that level #####
  
  library(here)
  library(bipartite)
  library(ggplot2)
  library(LOTUS)
  
  
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
  netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))
  
  names(netlists) <- names(rawnets)
  
  netlists <- lapply(netlists, function(x) lapply(x, function(y) apply(y, 2, as.numeric)))
  
  #Specify which index(s) to use
  
  ind <- c('functional complementarity',
           'mean number of shared partners',
           'niche overlap',
           'discrepancy',
           'NODF', 'modularity')
  
  
  
  
  
  m <- metcalcs(networks= netlists, indices = ind, network_level = 'higher')
  
  write.csv(m, 'data/output_data/all_bats/sitewise_otu_conclusions.csv')
  m <- read.csv('data/output_data/all_bats/sitewise_otu_conclusions.csv')
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  } #A function to capitalise the metric names when making plots
  
  
  
  
  m$site <- unlist(lapply(m$network, function(x) strsplit(as.character(x), split = ', ')[[1]][1]))
  #m$year <- unlist(lapply(m$network, function(x) strsplit(as.character(x), split = ', ')[[1]][2]))
  m$site <- gsub('DANUM', 'Danum', m$site)
  m$site <- gsub('MALIAU', 'Maliau', m$site)
  # m$habitat_type <- as.factor(ifelse(m$site == 'SAFE', 'Logged', 'Primary'))
  # m$metric <- gsub(' ', '\n', m$metric)
  # m$metric <- gsub('\\.', '\n', m$metric)
  # m$metric <- gsub('\nHL', '', m$metric)
  # 
  # m$metric <- gsub('mean\nnumber\nof\nshared\npartners', 'mean number of\nshared partners', m$metric)
  m$metric <- gsub('\\.', ' ', m$metric)
  m$metric <- gsub(' HL', '', m$metric)
  m$metric <- firstup(m$metric)
   sitescatter <- ggplot(m , aes(x = clustering, y = value, color = site)) +
     geom_point()+
     labs(x = 'Clustering (%)', y = NULL) +
     geom_smooth(method = lm, se = T)+
     scale_x_continuous(breaks = seq(91, 98, 1))+
     scale_color_brewer(type = 'qual')+
     facet_wrap(~ firstup(gsub('\\.', ' ', metric)), scales = 'free')+
     theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"))+#strip stuff sorts the facet labels, spacing adjusts the space between facets
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))
# ## @knitr otu_scatterplot
   sitescatter
#   
## @knitr otu_lineplot
  
  line_plot(input = m, metric = 'metric', network = 'network', clustering = 'clustering', value = 'value')  
  
## @knitr plot_saving  
  # sitescatter
  # pdf('plots/MOTU_sites_combined.pdf')
  # sitescatter
  # dev.off()
  # jpeg('plots/MOTU_sites_combined.jpg', width = 7, height = 7, units = 'in', res = 300)
  # sitescatter
  # dev.off()

  
  jpeg('plots/sitewise_lineplot.jpg', width = 7, height = 7, units = 'in', res = 300)
  line_plot(input = m, metric = 'metric', network = 'network', clustering = 'clustering', value = 'value')  
  dev.off()
  
  pdf('plots/sitewise_lineplot.pdf')
  line_plot(input = m, metric = 'metric', network = 'network', clustering = 'clustering', value = 'value')  
  dev.off()
  # 
  # 
mod_plot <-    ggplot(m[m$metric=='Modularity',] , aes(x = clustering, y = value, color = network)) +
     geom_point()+
     labs(x = 'Clustering (%)', y = NULL) +
     geom_smooth(method = lm, se = T)+
     scale_x_continuous(breaks = seq(91, 98, 1))+
     scale_color_brewer(type = 'qual')+
     theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"))+#strip stuff sorts the facet labels, spacing adjusts the space between facets
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf('plots/modularity.pdf')
mod_plot
dev.off()

