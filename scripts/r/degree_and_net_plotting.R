#### Header ####
## Project: Bat-diet
## Script purpose: Plotting degree distributions
## Date: 04/08/21
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(here)
library(tidyverse)
library(ggridges)
library(ggraph)
library(tidygraph)

source(here('scripts', 'r', 'r_network_gen.r'))

nets <- r_network_gen(lulu = T, filter_species = T, include_malua = F)

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))
       
for(n in 1:length(nets)){
  write.csv(x = nets[[n]],
            file = paste0('data/output_data/adjacency_matrices_for_SI/',
                          names(nets)[n], '.csv'))
}

# get degrees of each site
degree_list <- list()
for(i in 1:length(nets)){
  chosen_net <- nets[[i]]
  sitename <- names(nets)[i]
  bat_degrees <- colSums(chosen_net)
  motu_degrees <- rowSums(chosen_net)
  
  bat_df <- data.frame(nodename = names(bat_degrees), degree = bat_degrees, 
                       type = 'Bats', site = sitename)
  motu_df <- data.frame(nodename = names(motu_degrees), degree = motu_degrees, 
                        type = 'MOTU', site = sitename)
  
  degree_list[[i]] <- rbind(bat_df, motu_df)
}

degree_df <- bind_rows(degree_list) %>%
  mutate(habitat_type = ifelse(site == 'SAFE', 'Logged', 'Old-growth')) %>%
  # make the factor levels alphabetical
  mutate(site = fct_relevel(site, sort))

# make a simple histogram
ggplot(degree_df, aes(x = degree)) +
  geom_histogram(bins = 10) +
  facet_wrap(site ~ type, scales = 'free', ncol = 2) +
  theme_bw()

# try a ridgeplot, they look nicer
degree_plot <- ggplot(degree_df, aes(x = degree, y = fct_rev(site), fill = habitat_type)) +
  geom_density_ridges() +
  facet_wrap(type ~ ., scales = 'free', ncol = 1) +
  theme_ridges() + #improve formatting
  theme_bw()+ 
  # stop the top of the plot being cropped
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .9)))+
  # sort the fill colours
  scale_fill_cyclical(values = c("#d0ca9f", "#85d7da"), guide = "legend", 
                      name = "Habitat type") +
  xlab('Degree')+
  ylab('Site')+
  theme(legend.position = 'bottom',
        text = element_text(size=20))

degree_plot
ggsave('plots/all_degree.pdf', degree_plot, width  = 14)




# Network plotting --------------------------------------------------------



for(n in 1:length(nets)){
  net <- nets[[n]]
  netname <- names(nets)[n]
  molten_net <- melt(net)
  molten_net_2 <- molten_net[(which(molten_net[,3] != 0)),]
  bat_igraph <- graph.data.frame(molten_net_2, directed = F)
  V(bat_igraph)$type <- V(bat_igraph)$name %in% molten_net_2[,1]
  # make the node size proportionate to degree
  VS = sqrt(degree(bat_igraph)) 

  
  tidy_bats <- as_tbl_graph(bat_igraph, directed = T) %>%
    activate(nodes) %>%
    # make a weighting for the degree
    mutate(degree = centrality_degree(), weights = VS) %>%
    activate(nodes) %>%
    left_join(data.frame(node_names = V(bat_igraph)$name, 
                         node_type = ifelse(grepl('denovo', V(bat_igraph)$name), 
                                            'MOTU', 'Bat')),
              by = c('name' = 'node_names'))
  
  
  
  
  
  tidy_bat_graph <- ggraph(tidy_bats, layout = 'kk') + 
    geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE) + 
    geom_node_point(aes(size = degree, colour = node_type)) + 
    scale_colour_manual(values = c('yellow', '#01c2cd'),
                        name = 'Node type')+
    # improve the names of the legend items
    scale_size_continuous(name = 'Degree') +
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white')  +
    # put the legend at the bottom, on two rows
    theme(legend.position =  "bottom",
          legend.box = "vertical"#,
          #text=element_text(size=16,  family="ArialMT")
          )+ 
    labs(title=netname)
  
  ggsave(filename = paste0('plots/network_plots/', netname, '.jpeg'), 
         plot = tidy_bat_graph)
}

stats_list <- list()
for(i in 1:length(nets)){
  stats_list[[i]] <- tibble(Site = names(nets)[i],
                       Connectance = networklevel(web = nets[[i]], 
                                    index = 'connectance',
                                    ),
                       Generality = networklevel(web = nets[[i]], 
                                    index = 'generality', level = 'higher'),
                       Vulnerability = networklevel(web = nets[[i]], 
                                    index = 'vulnerability', level = 'lower'))
}

basic_stats <- bind_rows(stats_list) %>%
  arrange(Site) %>%
  write_csv(x = ., path = 'results/basic_network_stats.csv')


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}




# bray-curtis dissimilarity
for(i in 1:length(nets)){
  
  print(names(nets)[i])
  
  
  colnames(nets[[i]]) <- gsub("Hice", "Hipposideros cervinus", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Hiri", "Hipposideros ridleyi", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Hidi", "Hipposideros diadema", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Hidy", "Hipposideros dyacorum", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Kein", "Kerivoula intermedia", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Keha", "Kerivoula hardwickii", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Kepa", "Kerivoula papillosa", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Rhbo", "Rhinolophus borneensis", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Rhse", "Rhinolophus sedulus", colnames(nets[[i]]))
  colnames(nets[[i]]) <- gsub("Rhtr", "Rhinolophus trifoliatus", colnames(nets[[i]]))
  
  # order the network's columns alphabetically, for happier plotting
  desired_net <- nets[[i]]
  desired_net <- desired_net[,order(colnames(desired_net))]
  
  vegandiv <- vegan::vegdist(t(desired_net), 'bray')
  
  dist_df <- vegandiv %>% 
    as.matrix() %>% 
    get_upper_tri() %>%
    melt(na.rm = T) %>%
    filter(value != 0)
    
  
  
  distplot <- ggplot(dist_df, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(Var2, Var1, label = round(value,3))) +
    scale_fill_gradient(low = 'blue', high = 'red',
                        name = 'Bray-curtis dissimilarity',
                        breaks = c(0.8, 0.9, 1),
                        limits = c(0.8, 1))+ 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'bottom',
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                 face = "italic"),
      axis.text.y = element_text(face = 'italic'), 
      text = element_text(size=20)
    )+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    ggtitle(names(nets)[i])
  
  ggsave(filename = paste0('plots/bray_', names(nets[i]), '.pdf'), distplot,
         width =12)
}

