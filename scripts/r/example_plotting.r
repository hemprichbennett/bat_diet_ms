#### Header ####
## Project: 
## Script purpose: 
## Date: 
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################
options(warnPartialMatchArgs = TRUE)

#
library(tidyverse)
# library(ggplotify)
# library(igraph)
# library(gridExtra)
# library(grid)

# Setup -------------------------------------------------------------------


# Load in the external data
nest <- read.csv("data/dummy_data/nested_2.csv")




# Functions ---------------------------------------------------------------

# A function which turns the input matrix into a nicely formatted edgelist,
# which ggplot will happily play with
matrix_formatting <- function(mat) {
  out <- mat %>%
    # Rename the prey species column
    rename(`Prey species` = X) %>%
    # make data long format, instead of wide
    gather(-`Prey species`, key = "Predator species", value = "Interacts") %>%
    # turn the interaction column into T/F instead of binary
    # mutate(Interacts = ifelse(Interacts == 1, T, F)) %>%
    # Reorder the factors, as ggplot is plotting them wrong
    mutate(`Prey species` = fct_rev(`Prey species`))


  return(out)
}


# A function to make a random edgelist with the exact same number of nodes and
# connectance as the input edgelist
non_fun <- function(x) {
  
  # Decide which rows should have interactions
  to_fill <- sample(
    seq(1, length(x$Interacts)),
    sum(x$Interacts)
  )
  
  output <- x %>%
    # Make all the interactions empty before adding the new ones
    mutate(Interacts = 0)
  
  # Add the new interactions
  output$`Interacts`[to_fill] <- 1
  
  return(output)
}


# Bipartite data manipulation ---------------------------------------------

# set up the object to be used for plotting modularity and nestedness

Nested <- matrix_formatting(nest)
#Modular <- matrix_formatting(mod)

set.seed(3)

`Non-nested` <- non_fun(Nested)
#`Non-modular` <- non_fun(Modular)

outlist <- list()
nets <- list(
  #Modular, `Non-modular`, 
  Nested, `Non-nested`
  )
names(nets) <- c(
  #'Modular', 'Non-modular', 
  'Nested', 'Non-nested'
  )

for (i in 1:length(nets)) {
  nets[[i]]$dataset <- names(nets)[i]
  
}

all_df <- do.call(rbind, nets) %>%
  # Turn the interaction column from a binary variable to a T/F
  mutate(Interacts = ifelse(Interacts == 1, T, F)) %>%
  # Now turn it into a factor, where T is before F
  mutate(Interacts = fct_rev(as.factor(Interacts))) %>%
  # Rename the predators and prey for the plotting
  mutate(`Predator species` = gsub('Pred', 'Predator ', `Predator species`),
         `Prey species` = gsub('Bug', 'Prey ', `Prey species`)) %>%
  # Reverse the factor orders for prey species, as ggplot plots them from bottom
  # to top, not vice versa as would be helpful here
  mutate(`Prey species` = fct_rev(`Prey species`)) %>%
  mutate(dataset = gsub('Nested', 'A) Nested', dataset),
        dataset = gsub('Non-nested', 'B) Non-nested', dataset))


# Plot the bipartite networks ---------------------------------------------

bipartite_nets <- ggplot(all_df, aes(x = `Predator species`, y = `Prey species`, fill = Interacts)) +
  # Add the tiles. 'Colour = black' adds gridlines
  geom_tile(colour = "black") +
  scale_fill_manual(values = c("black", "#f2f2f2"),
                    labels = c('True', 'False')) +
  facet_wrap(. ~ dataset, strip.position = 'bottom') +
  theme_classic() +
  theme(
    # Put the legend at the bottom
    legend.position = "bottom",
    # # Get rid of the axis labels
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = -0.1),
    # Get rid of the background colours on facet labels
    strip.background = element_rect(
      color="black", fill="white", size=1, linetype="solid"),
    panel.spacing = unit(2, "lines"),
    text = element_text(size=20)
  ) +
  scale_x_discrete(position = "top") 

bipartite_nets
ggsave('example_images/biparite_networks.png', bipartite_nets,
       width = 9.5, height = 9)
