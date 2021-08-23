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
library(bipartite)


source('scripts/r/r_network_gen.r')

nets <- r_network_gen(lulu=T, filter_species = T)

# the metrics which are asked for in networklevel
ind <- c('functional complementarity', 'weighted NODF', 'discrepancy')

field_data <- read_csv('data/Edited_all_files_with_lat_long_VKedits.csv')

# format the field data for this analysis
field_df <- field_data %>%
  select(Species, Site) %>%
  mutate(Site = gsub('DVCA', 'DANUM', Site)) %>%
  filter(Site %in% c('SAFE', 'DANUM', 'MALIAU')) %>%
  group_by(Species, Site) %>%
  summarise(nbats = n()) %>%
  filter(nbats >2) %>%
  arrange(Site, nbats, Species) %>%
  group_by(Site) %>%
  mutate(spec_rank = rank(nbats, ties.method = 'first'))# %>%


out_list <- list()
iteration <- 1
for(net_name in names(nets)){
  
  net <- nets[[which(names(nets) == net_name)]]
  
  relative_abundance <- field_df %>% 
    ungroup() %>%
    filter(Site == net_name,
           Species %in% colnames(net))
  
  
  for(n_bats_kept in c(3:nrow(relative_abundance))){
    
    # make a vector of the species to be retained (the n most abundant ones)
    to_keep <- relative_abundance %>% 
      slice_max(spec_rank, n = n_bats_kept) %>% 
      select(Species) %>%
      pull()
    
    new_net <- net[, which(colnames(net) %in% to_keep)]
    
    fun_comp <- networklevel(new_net, index = 'functional complementarity', 
                             level = 'higher')
    WNODF <- networklevel(new_net, index = 'weighted NODF')
    
    Discrepancy <- networklevel(new_net, index = 'discrepancy', 
                                level = 'higher')
    
    out_list[[iteration]] <- tibble(n_bats_kept, net_name, fun_comp, WNODF, 
                                    Discrepancy)
    iteration <- iteration + 1
  }
}

out_df <- bind_rows(out_list)

out_df <- out_df %>% pivot_longer(cols = c('WNODF', 'fun_comp', 'Discrepancy'),
                        names_to = 'metric', values_to = 'val') %>%
  mutate(net_name = gsub('DANUM', 'Danum', net_name),
         net_name = gsub('MALIAU', 'Maliau',  net_name),
         metric = gsub('fun_comp','Functional \ncomplementarity',  metric))



out_plot <- ggplot(out_df, aes(x = n_bats_kept, y = val, colour = net_name)) +
  geom_point(aes(size = 10, alpha = 0.8)) +
  facet_wrap(metric ~ ., scales = 'free',
             strip.position = 'left', ncol = 1,) +
  scale_colour_viridis_d(name = 'Network') + 
  theme_classic()+
  scale_x_continuous(breaks = seq(min(out_df$n_bats_kept), 
                                  max(out_df$n_bats_kept), 1)) +
  theme(legend.position = 'bottom',
        strip.placement = 'outside',
        strip.background = element_rect(color = "white", size = 1),
        text=element_text(size=14))+
  xlab("Number of the most abundant bats included in network")+
  ylab(element_blank())+ 
  # remove dot size from legend
  guides(size = F, alpha = F)

out_plot
ggsave('plots/n_abundant_species.pdf', out_plot,
       height = 11, width = 6)
