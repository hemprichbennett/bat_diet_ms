#### Header ####
## Project: bat diet
## Script purpose: calculating z-scores from the random files created by 
## array_site_motu_95_confidence_intervals
## Date: 2021-01-06
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

#### Setup ####

# Prevent partial-match errors 
options(warnPartialMatchArgs = TRUE)



library(tidyverse)

# Load and format the random data -----------------------------------------

rand_filenames <- list.files('data/output_data/null_values/',
                          full.names = T)

rand_vals <- lapply(rand_filenames, read_csv) %>%
  bind_rows() 



# Get the real values -----------------------------------------------

source('scripts/r/r_network_gen.r')

inpath <- 'data/processed_dna_data/lulu/'
filenames <- list.files(pattern = '.csv', path = inpath)

filenames
filenames <- paste(inpath, filenames, sep = '')


rawnets <- lapply(filenames, read.csv, header = F, stringsAsFactors = F, row.names=1)
names(rawnets) <- gsub('.*\\/', '', filenames)
names(rawnets) <- gsub('_.+', '', names(rawnets))

for(i in 1:length(rawnets)){
  rawnets[[i]][2:nrow(rawnets[[i]]),] <- ifelse(rawnets[[i]][2:nrow(rawnets[[i]]), ] == 0, 0, 1)
}

netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))

names(netlists) <- names(rawnets)



actual_vals_list <- list()
z <- 1
for(ind in unique(rand_vals$metric_used)){
  for(i in 1:length(netlists)){
    
    vals_list <- lapply(netlists[[i]], function(x) networklevel(x, index = ind, level = 'higher'))
    
    
    
    out_df <- bind_rows(vals_list) %>% 
      rename_at(1, ~"actual" ) %>%
      mutate(site = names(vals_list), 
             metric_used = ind, 
             clustering_level = names(netlists)[i]) %>%
      mutate(clustering_level = as.numeric(clustering_level))
    
    
    actual_vals_list[[z]] <- out_df
    z <- z + 1
  }
  
}

actual_vals <- bind_rows(actual_vals_list) 


# Calculate 'z-scores' ----------------------------------------------------




summary_vals <- rand_vals %>%
  group_by(site, metric_used, clustering_level, fixedmargins) %>%
  summarise(mean_rand = mean(metric_value),
            sd_rand = sd(metric_value))


z_vals <- summary_vals %>%
  ungroup() %>%
  # combine the datasets
  left_join(actual_vals) %>%
  # calculate the 'z' score
  mutate(z = (actual - mean_rand)/sd_rand) %>%
  # make some of the variables prettier for the outputs
  mutate(site = gsub('DANUM', 'Danum', site),
         site = gsub('MALIAU', 'Maliau', site),
         metric_used = gsub('functional complementarity', 'Functional complementarity',
                       metric_used),
         metric_used = gsub('weighted NODF', 'WNODF', metric_used))

# write a csv of all of the values in their raw form
write_csv(z_vals, 'results/z_scores/all_z_scores.csv')


# then make a csv for each metric, with the values given to 3 decimal places
z_vals %>%
  select(-mean_rand, -sd_rand) %>%
  mutate(actual = round(actual, digits = 2),
         z = round(z, digits = 3),
         Treatment =ifelse(site == 'SAFE', 'Logged', 'Old growth')) %>%
  rename(Site = site, Metric = metric_used, 
         `Clustering threshold` = clustering_level, 
         `Observed Value` = actual) %>%
  group_by(Metric) %>%
  do(write_csv(., paste0('results/z_scores/', unique(.$Metric), "_zscores.csv")))


# Finally, one table to rule them all...

met_list <- list()
cols <- c('actual', 'z')
for(met in unique(z_vals$metric_used)){
  met_list[[met]] <- z_vals %>%
    filter(metric_used == met) %>%
    rename_at(cols, list( ~paste( ., met, sep = '_') ) ) %>%
    select(-metric_used, - mean_rand, - sd_rand)
}

mets_df <- left_join(met_list[[1]], met_list[[2]], 
                     met_list[[3]], by = c('site', 'clustering_level'))

net_list <- list()
cols <- c('actual', 'both', 'columns')
for(net in unique(z_vals$site)){
  net_list[[net]] <- z_vals %>%
    filter(site == net) %>%
    select(-site, - mean_rand, - sd_rand) %>%
    mutate(actual = round(actual, digits = 3),
           z = round(z, digits = 3)) %>%
    pivot_wider(names_from = fixedmargins, values_from = z) %>%
    rename_at(cols, list( ~paste( ., net, sep = '_') ) )
}



nets_df <- left_join(net_list[[1]], net_list[[2]],
                     by = c('metric_used', 'clustering_level')) %>%
            left_join(net_list[[3]], by = c('metric_used', 'clustering_level'))

write_csv(nets_df, path = 'results/z_scores/grouped_by_metric.csv')

# Plot --------------------------------------------------------------------


ggplot(rand_vals, aes(x = clustering_level, y = metric_value)) +
  geom_violin() +
  facet_grid(metric_used ~ site, scales = 'free') +
  theme_classic()

for_plot <- z_vals %>%
  pivot_longer(cols = c('actual', 'z'), names_to = 'var_type', 
               values_to = 'var_value') %>%
  mutate(metric = gsub('Functional complementarity', 
                       'Functional\ncomplementarity', metric_used))

actual_points <-   ggplot(filter(for_plot, var_type == 'actual'), 
         aes(x = clustering_level, y = var_value, colour = site)) +
  geom_point() +
  scale_colour_viridis_d()+
  scale_x_continuous(breaks = c(91:98))+
  facet_wrap(. ~ metric_used, scales = 'free_y', ncol = 1,
             # sort the facet label placement
             strip.position = 'left') +
  theme_bw() +
    theme(legend.position = 'none',
          # sort the facet label placement
          strip.placement = "outside",
          strip.text.y = element_text(size = 8),
          #axis.title.y=element_blank(),
          axis.title.x=element_blank())+
  ylab('Observed value')

actual_points

z_points <-   ggplot(filter(for_plot, var_type == 'z'), 
                          aes(x = clustering_level, y = var_value, 
                              colour = site)) +
  geom_point(size = 2) +
  scale_colour_viridis_d()+
  scale_x_continuous(breaks = c(91:98))+
  facet_wrap(. ~ metric_used, ncol = 1) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        # remove facet label
        strip.text = element_blank(),
        #axis.title.y=element_blank()
        ) +
  ylab('Z-value')

z_points

# Function for extracting the legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- ggplot(for_plot, 
                 aes(x = clustering_level, y = var_value, colour = site)) +
  geom_point() +
  scale_colour_viridis_d(name = 'Network')+
  theme_bw()+
  theme(legend.position = 'bottom')

legend <- g_legend(legend)


pdf('plots/z_plot_unedited.pdf')
grid_plot <- gridExtra::grid.arrange(actual_points, z_points,
                                     grid::textGrob('Clustering threshold'),
                                     legend,
                                     layout_matrix = rbind(c(1,2),
                                                           c(3,3),
                                                           c(4,4)),
                                     heights = unit(c(3.3, 0.3, 0.3), c('in', 'in', 'in')))

dev.off()



# Alternative z and observed plot -----------------------------------------

for_alt_z_plot <- z_vals %>%
  #select only necessary columns
  select(-mean_rand, -sd_rand) %>%
  # start making it longer for ggplot
  pivot_longer(cols = c('actual', 'z'), names_to = 'plotting_var_name', 
               values_to = 'plotting_value') %>%
  unite(column_category, fixedmargins, plotting_var_name, sep = '_') %>%
  # this has given us duplicates for the 'actual' value, as we have 'both_actual'
  # and 'columns_actual', both of which are crap and unnecessary. If we strip
  # the text before the underscore, we can then remove duplicate rows easily,
  # as the duplicate rows will then be identical
  mutate(column_category = gsub('.+_actual', 'actual', column_category)) %>%
  distinct()
  
  

ggplot(for_alt_z_plot, aes(x = column_category, y = clustering_level))+ 
  geom_tile(aes(fill = plotting_value))+
  geom_text(aes(label = round(plotting_value, 3))) +
  facet_grid(site ~ metric_used)


# Now for the ranges plot -------------------------------------------------

vals_for_range_plot <- rand_vals %>%
  # group by relevant variables before we calculate the quantiles
  group_by(metric_used, clustering_level, site, fixedmargins) %>%
  # calculate the quantiles of the random values
  summarise(min_quantile = quantile(x = metric_value, probs = 0.025),
            max_quantile = quantile(x = metric_value, probs = 0.975)) %>%
  # now add the actual observed values to the tibble
  left_join(actual_vals) %>%
  # we only want to plot the dots if they fall outside of the random ranges,
  # so we need to make a column saying if they fall within the range or not
  mutate(to_plot = ifelse(actual >= min_quantile &
                            actual <= max_quantile, 
                          F, T),
         # now we need a column where the observed value is only present if it's
         # outside the ranges
         for_plotting = ifelse(to_plot == T, actual, NA)) %>%
  # make the names of the variables nicer
  mutate(site = gsub('DANUM', 'Danum', site),
         site = gsub('MALIAU', 'Maliau', site),
         metric_used = gsub('discrepancy', 'Discrepancy', metric_used),
         metric_used = gsub('weighted NODF', 'Weighted NODF', metric_used),
         metric_used = gsub('functional complementarity', 'Functional Complementarity', metric_used),
         fixedmargins = gsub('both', 'Both margin\nsums retained', fixedmargins),
         fixedmargins = gsub('columns', 'Column\nsums retained', fixedmargins))


# palette to use

cbPalette <- c("#5a3fa8",
               "#d4c200",
               "#a4c7ff")


# pasted in from other script, needs editing!
#plot
metrics_facet <- ggplot(vals_for_range_plot, aes(y =clustering_level, colour = site))+
  geom_errorbarh(aes(xmin=min_quantile, xmax=max_quantile, colour = site),
                 height = 0.4, alpha = 0.8, show.legend = F)+
  geom_point(aes(y = clustering_level, x = for_plotting, size = 1.4))+
  facet_grid(fixedmargins ~ metric_used, scales = 'free_x'#, ncol = 2
             )+
  scale_colour_manual(values=cbPalette, name = 'Observed network\nvalue')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.placement = "outside", 
        strip.background =element_rect(fill="white", colour = 'white'),
        text = element_text(size=20),
        legend.position="bottom")+
  labs(y = 'Clustering %', x = NULL) + 
  # tell ggplot to plot the colour in the legend but not the size
  guides(colour = "legend", size = "none") +
  # make the points in the legend bigger
  guides(colour = guide_legend(override.aes = list(size=4.5)))
metrics_facet

ggsave('plots/randomized_ranges/both_nullmodels.pdf', metrics_facet, width = 14)

