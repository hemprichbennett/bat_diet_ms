#### Header ####
## Project: Bat-diet
## Script purpose: Getting inter-site-and year beta-diversity
## Date: 04/08/21
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(here)
library(tidyverse)


source(here('scripts', 'r', 'r_network_gen.r'))

nets <- r_network_gen(lulu = T, filter_species = T, include_malua = F, 
                      split_by = 'site and year')

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))
# reorder the list alphabetically, for happier plotting later
nets <- nets[order(names(nets))]

# make some basic stats from vegan
motu_degrees_list <- list()
for(n in 1: length(nets)){
  netname <- names(nets)[n]
  motu_tib <- tibble(motu_name = rownames(nets[[n]]),
                     degree = as.numeric(rowSums(nets[[n]])))
  colnames(motu_tib)[2] <- netname
  motu_degrees_list[[n]] <- motu_tib
}

motu_site_matrix <- motu_degrees_list %>%
  reduce(full_join, by = "motu_name") %>% 
  replace(is.na(.), 0) %>%
  as.matrix() %>%
  t()

colnames(motu_site_matrix) <- motu_site_matrix[1,]
motu_site_matrix <- motu_site_matrix[-1,]

class(motu_site_matrix) <- "numeric"

vegandiv <- betadiver(motu_site_matrix, "w")

# functipon to get the upper triangle of a correlation matix, taken from
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

distplot <- vegandiv %>% 
  as.matrix() %>% 
  get_upper_tri() %>%
  melt(na.rm = T) %>%
  filter(value != 0) %>%
  ggplot(aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(Var2, Var1, label = round(value,3))) +
  scale_fill_gradient(low = 'blue', high = 'red',
                      name = 'MOTU beta-diversity')+ 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
  
distplot

ggsave('plots/motu_betadiversity.pdf')


