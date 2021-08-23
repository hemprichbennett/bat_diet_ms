#### Header ####
## Project: bat-diet
## Script purpose: Finding the effect of removing different NUMBERS of bat species from our networks
## Date: 04/01/21
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: this is a different analysis to that which I did in my thesis, which 
## was to identify the effects of removing specific species. Here I look at the
## effect of the number of bat species present in a network
##################################################

library(tidyverse)

species_subset <- function(net_no, params_index){
  # the parameters from the list to use
  to_use <- bat_params[[params_index]]
  
  # calculate the values which come from the 'networklevel' function
  netlevel_vals <- networklevel(nets[[net_no]][,to_use],
                                index = ind, level = 'higher') %>%
    t() %>%
    data.frame()
  
  # calculate the modularity value and add it to the df
  # netlevel_vals$modularity <- slot(computeModules(web = nets[[net_no]][,to_use]), 
  #                                  'likelihood')
  # return it
  return(netlevel_vals)
}


source('scripts/r/r_network_gen.r')

nets <- r_network_gen(lulu=T, filter_species = T)

# the metrics which are asked for in networklevel (modularity has to be called 
# by a separate command in the species_subset function)
ind <- c('functional complementarity', 'weighted NODF', 'discrepancy')

min_bats <- 3

i <- 1

out_list <- list()

for(net_number in 1:length(nets)){
  
  
  net_name <- names(nets)[net_number]
  
  nbats_to_keep <- seq(min_bats, ncol(nets[[net_number]]))
  
  for(nbats in nbats_to_keep){
    
    # the below gives all the possible subsets of bat columns, 
    # where the number of bat SPECIES retained is nbats_to_keep. Its in a list, 
    # and each list element is a vector of potential column numbers
    bat_params <- combn(seq(1, ncol(nets[[net_number]])),
                        m = nbats,
                        simplify = F)
    
    
    
    
    metrics_df <- lapply(seq(1,length(bat_params)), function(x) 
      species_subset(net_no = net_number, params_index = x)) %>%
      bind_rows()
    
    metrics_df$nbats_kept <- nbats
    metrics_df$network <- net_name
    
    out_list[[i]] <- metrics_df
    
    i <- i + 1
    
  }
  
  
}

out_df <- bind_rows(out_list)  

write_csv(out_df, 'results/nspecies_removals.csv')

# reformat the df to make it long-form for plotting
results_tib <- out_df %>%
  filter(nbats_kept < 10) %>%
  pivot_longer(c('weighted.NODF', 'functional.complementarity.HL', 
                 'discrepancy.HL'), names_to = 'metric_name', 
               values_to = 'metric_value') %>%
  mutate(metric_name = gsub('weighted.NODF', 'Weighted NODF', metric_name),
         metric_name = gsub('functional.complementarity.HL', 
                            'Functional complementarity', metric_name),
         metric_name = gsub('discrepancy.HL', 'Discrepancy', metric_name),
         network = gsub('DANUM', 'Danum', network),
         network = gsub('MALIAU', 'Maliau', network)) 
  


nspecies_plot <- ggplot(results_tib, aes(x = network, y = metric_value, fill = network)) +
  geom_violin() +
  theme_classic()+
  facet_grid(metric_name ~ nbats_kept , scales = 'free_y',
             # move the facet labels from default position
             switch = 'both',) +
  scale_fill_viridis_d(name = 'Network')+ 
  theme(legend.position = 'bottom',
        # put facet labels 'outside'
        strip.placement = "outside",
        strip.background = element_rect(
          colour = 'white'),
        # get rid of the axis labels, as otherwise it will have the network name
        # repeating
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(size = 20)) +
  xlab('Number of bat species in network')

nspecies_plot  
  
ggsave('plots/nspecies_plot.jpg', nspecies_plot, dpi = 700, width = 14)  
