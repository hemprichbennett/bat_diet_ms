##################################################
## Project: All bats
## Script purpose: Calculating individual-based metrics and taxonomic composition for individual bats
## Date: 12/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: This is a modified version of the all_species_individual_analysis script, here the analysis looks at the sites with years combined,
## rather than splitting the metaweb by site AND year
##################################################
if (interactive() == TRUE) {
  library("here")
  library(ggplot2)
  library(ggridges)
  library(reshape2)
  library(forcats)
  library(dplyr)
} else {
  library(here, lib.loc = "/data/home/btw863/r_packages/")
}

setwd(here())
source("scripts/r/r_network_gen.r")

field_data <- read.csv(here("data/Edited_all_files_with_lat_long_VKedits.csv"))
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ", ")
field_data$Faeces_no1 <- gsub("T", "", field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub("T", "", field_data$Faeces_no2)
field_data$Site <- gsub("DVCA", "Danum", field_data$Site)
field_data$Site <- gsub("DANUM", "Danum", field_data$Site)
field_data$Site <- gsub("MALIAU", "Maliau", field_data$Site)

all_interactions <- r_network_gen(collapse_species = F, desired_species = NULL, include_malua = F, filter_species = T, lulu = T)

desired_cols <- c("MOTU", "DANUM", "MALIAU", "SAFE")

all_interactions <- all_interactions[, -which(!all_interactions[1, ] %in% desired_cols)]


colnames(all_interactions) <- all_interactions[2, ]
all_interactions <- all_interactions[-c(2), ]
rownames(all_interactions) <- all_interactions[, 1]
all_interactions <- all_interactions[, -1]

locations <- c()

for (i in 1:ncol(all_interactions)) {
  locations[i] <- all_interactions[1, i]
  names(locations)[i] <- all_interactions[1, i]
}


row_names <- rownames(all_interactions) # Store these as an object as the apply below kills them
all_interactions <- apply(all_interactions, 2, as.numeric)
rownames(all_interactions) <- row_names

prey_data <- read.csv("data/taxonomy/order.csv")
colnames(prey_data) <- c("MOTU", "Taxa")
prey_data$MOTU <- as.character(prey_data$MOTU)




### Add the taxonomic information to the interactions for everything
taxa_mat <- matrix(nrow = 0, ncol = ncol(all_interactions))
colnames(taxa_mat) <- colnames(all_interactions)
z <- 1
for (i in 1:nrow(all_interactions)) {
  rowname <- rownames(all_interactions)[i]
  if (!rowname %in% prey_data$MOTU) {
    next()
  }
  tax <- as.character(prey_data[which(prey_data$MOTU == rowname), "Taxa"])
  if (is.null(nrow(taxa_mat))) { # If its the first iteration there won't be any rownames yet, so the next if statement will fail
    taxa_mat <- rbind(taxa_mat, as.numeric(all_interactions[i, ]))
    rownames(taxa_mat)[z] <- tax
    z <- z + 1
  }
  if (tax %in% rownames(taxa_mat)) {
    to_merge <- which(rownames(taxa_mat) == tax)
    taxa_mat[to_merge, ] <- taxa_mat[to_merge, ] + as.numeric(all_interactions[i, ])
  } else {
    taxa_mat <- rbind(taxa_mat, as.numeric(all_interactions[i, ]))
    rownames(taxa_mat)[z] <- tax
    z <- z + 1
  }
}

##### Make a list with a network for each site####

sites_list <- list()





for (i in 1:length(unique(names(locations)))) {
  loc <- unique(names(locations))[i]
  sites_list[[i]] <- all_interactions[, which(names(locations) == loc)]
}
names(sites_list) <- unique(names(locations))

names(sites_list)

##### Do some ecology ####
all_ecology <- matrix(nrow = 0, ncol = 1 + ncol(specieslevel(matrix(sample(0:1, size = 100, replace = T), nrow = 10, ncol = 10), level = "higher")))

for (i in 1:length(sites_list)) {
  starttime <- Sys.time()
  if (length(which(duplicated(colnames(sites_list[[i]])))) > 0) {
    sites_list[[i]] <- sites_list[[i]][, -which(duplicated(colnames(sites_list[[i]])))]
  }
  splevel <- specieslevel(sites_list[[i]], level = "higher")
  endtime <- Sys.time()
  cat(names(sites_list)[i], "took", endtime - starttime, "\n")
  all_ecology <- rbind(all_ecology, cbind(rep(names(sites_list)[i], nrow(splevel)), splevel))
}
all_ecology <- cbind(rownames(all_ecology), all_ecology)
# write.csv(all_ecology, 'all_ecology_1.csv')
# all_ecology <- read.csv('all_ecology_1.csv')
##### Add some taxonomic information
write.csv(taxa_mat, "shiny/order_composition/order_mat.csv")

colnames(all_ecology)[c(1:2)] <- c("Sample", "Site")
t_taxa <- t(taxa_mat)
t_taxa <- cbind(t_taxa, rownames(t_taxa))
t_taxa <- as.data.frame(t_taxa)
colnames(t_taxa)[ncol(t_taxa)] <- "Sample_no"

all_ecology <- merge(all_ecology, t_taxa, by.x = "Sample", by.y = "Sample_no")

all_ecology <- merge(x = all_ecology, y = field_data, by.x = "Sample", by.y = "Faeces_no1")

colnames(all_ecology)[which(colnames(all_ecology) == "Site.y")] <- "Site"


##### Work out the nestedness of each bat, then add it to the df#####



write.csv(all_ecology, "data/output_data/all_bats/sitewise_all_individual_info.csv")
save.image("data/output_data/sitewise_all_individual_info.RDS")
##### Local work ####
load("data/output_data/sitewise_all_individual_info.RDS")


# Remove SBE (aka Malua) as we don't want to use it for the manuscript

all_ecology <- all_ecology[-which(all_ecology$Site == "MALUA"), ]

# Look at the occurence of taxa in each species

tax_df <- all_ecology[, c(1, 2, 66, seq(23, 47))]
tax_df <- melt(tax_df, id.vars = c("Sample", "Species", "Site.x"))
colnames(tax_df)[c(3:5)] <- c("Site", "Order", "Present/absent")
# tax_df$`Present/absent` <- as.integer(tax_df$`Present/absent`)
tax_df$`Present/absent` <- ifelse(tax_df$`Present/absent` == 0, 0, 1)

prop_present <- sapply(unique(tax_df[, c("Species", "Site", "Order")]), function(x) as.character(x))
prop <- c()
nbats <- c()
# for(i in 1: 1){
for (i in 1:nrow(prop_present)) {
  tem <- tax_df[which(tax_df$Species == as.character(prop_present[i, 1]) & tax_df$Site == as.character(prop_present[i, 2])
  & tax_df$Order == as.character(prop_present[i, 3])), ]
  prop <- c(prop, sum(tem$`Present/absent`) / nrow(tem)) # The number of bats that consumed the order, divided by total bats
  nbats <- c(nbats, nrow(tem))
}


prop_present <- cbind(prop_present, prop, nbats)
prop_present <- as.data.frame(prop_present)
prop_present$prop <- as.numeric(as.character(prop_present$prop))
prop_present$nbats <- as.integer(as.character(prop_present$nbats))

prop_present$Site <- gsub("DANUM", "Danum", prop_present$Site)
prop_present$Site <- gsub("MALIAU", "Maliau", prop_present$Site)
prop_present$Species <- gsub("Hice", "Hipposideros\ncervinus", prop_present$Species)
prop_present$Species <- gsub("Hiri", "Hipposideros\nridleyi", prop_present$Species)
prop_present$Species <- gsub("Hidi", "Hipposideros\ndiadema", prop_present$Species)
prop_present$Species <- gsub("Hidy", "Hipposideros\ndyacorum", prop_present$Species)
prop_present$Species <- gsub("Kein", "Kerivoula\nintermedia", prop_present$Species)
prop_present$Species <- gsub("Keha", "Kerivoula\nhardwickii", prop_present$Species)
prop_present$Species <- gsub("Kepa", "Kerivoula\npapillosa", prop_present$Species)
prop_present$Species <- gsub("Rhbo", "Rhinolophus\nborneensis", prop_present$Species)
prop_present$Species <- gsub("Rhse", "Rhinolophus\nsedulus", prop_present$Species)
prop_present$Species <- gsub("Rhtr", "Rhinolophus\ntrifoliatus", prop_present$Species)



tiles <- ggplot(data = prop_present[which(prop_present$nbats > 5), ], aes(y = fct_rev(Order), x = Site)) + geom_tile(aes(fill = prop), colour = "white") +
  scale_fill_gradient2(
    low = "white", mid = 'blue',
    high = "black", name = "Proportion of bats consuming",
    breaks = c(0, 0.25, 0.5, 0.75, 1), midpoint = 0.5
  ) +
  labs(x = "Site", y = "Prey order") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray")) +
  theme(strip.text.x = element_text(size = 12)) +
  facet_wrap(~Species) +
  theme(
    legend.position = "bottom", legend.justification = "left", strip.background = element_rect(fill = "white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"),
    strip.text = element_text(face = "italic")
  ) # strip stuff sorts the facet labels, spacing adjusts the space between facets

tiles

pdf("plots/sitewise_proportion_of_bats_containing.pdf", height = 12)
tiles
dev.off()

jpeg("plots/sitewise_proportion_of_bats_containing.jpg", units = "in", width = 7.5, height = 12, res = 300)
tiles
dev.off()


degree_df <- all_ecology %>%
  select(Site, degree, Species)

# make some basic summary numbers
sp_summaries <- degree_df %>%
  group_by(Site, Species) %>%
  summarise(mean_degree  = mean(degree),
            sd_degree = sd(degree))

colnames(degree_df) <- gsub("all_ecology\\.", "", colnames(degree_df))
degree_df$genus <- rep(NA, nrow(degree_df))
degree_df$genus[grepl("Hi", degree_df$Species)] <- "Hipposideros"
degree_df$genus[grepl("Rh", degree_df$Species)] <- "Rhinolophus"
degree_df$genus[grepl("Ke", degree_df$Species)] <- "Kerivoula"
degree_df$hab_type <- rep(NA, nrow(degree_df))
degree_df$hab_type[grepl("Danum", degree_df$Site)] <- "Old growth"
degree_df$hab_type[grepl("Maliau", degree_df$Site)] <- "Old growth"
degree_df$hab_type[grepl("SAFE", degree_df$Site)] <- "Logged"
degree_df$hab_type[grepl("SBE", degree_df$Site)] <- "Logged, replanted"
degree_df$degree <- as.integer(degree_df$degree)


library(plm)
library(car) # Companion to applied regression
library(magrittr)

degree_df$Species %<>%
  gsub("Hice", "Hipposideros cervinus", .) %<>%
  gsub("Hidi", "Hipposideros diadema", .) %<>%
  gsub("Hidy", "Hipposideros dyacorum", .) %<>%
  gsub("Hiri", "Hipposideros ridleyi", .) %<>%
  gsub("Keha", "Kerivoula hardwickii", .) %<>%
  gsub("Kein", "Kerivoula intermedia", .) %<>%
  gsub("Kepa", "Kerivoula papillosa", .) %<>%
  gsub("Kepe", "Kerivoula pellucida", .) %<>%
  gsub("Rhbo", "Rhinolophus borneensis", .) %<>%
  gsub("Rhse", "Rhinolophus sedulus", .) %<>%
  gsub("Rhtr", "Rhinolophus trifoliatus", .)


# Make a table of samples used, for thesis

degree_df$SiteAndYear <- paste0(degree_df$Site, ", ", degree_df$Year)
sample_table <- table(degree_df$Species, degree_df$SiteAndYear)

write.csv(sample_table, "results/sample_table.csv")

# Run two fixed effects models, then compare them
fixed_site <- lm(degree ~ factor(Site) + factor(Species) - 1, data = degree_df)



fixed_hab <- lm(degree ~ factor(hab_type) + factor(Species) - 1, data = degree_df)



anova(fixed_site, fixed_hab)

# Run a fixed effect using site AND habitat type, then AIC it
both <- lm(degree ~ factor(hab_type) + factor(Site) + factor(Species) - 1, data = degree_df)
library(MASS)
step_mod <- MASS::stepAIC(both, 
  direction = "backward")

# step_mod shows us that the best-fitting model drops no terms. However we want to drop either
# hab_type or Site. Dropping Site reduces AIC the least, therefore we drop Site
# as dropping hab_type reduces the AIC of our model most.

summary(fixed_hab)


library(broom); library(dplyr)
model_df <- left_join(tidy(both, quick = TRUE),
                tidy(both, quick = FALSE),
                by = c("term", "estimate")) #Save the model output as a df
# code from https://stackoverflow.com/questions/47891753/exporting-lm-summary-output-to-dataframe-including-na

model_df$term <- gsub('.+)', '', model_df$term)

write.csv(model_df, 'results/degree_model_coefficients.csv')

sp_ridge <- ggplot(degree_df, aes(y = fct_rev(Site), x = degree, fill = fct_rev(hab_type))) +
  geom_density_ridges(scale = 0.85, panel_scaling = F,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x)) + # The scale determines the space between the rows
  theme_ridges() + # This changes the theme to make it more aesthetically pleasing
  scale_fill_cyclical(values = c("#d0ca9f", "#85d7da"), guide = "legend", name = "Habitat type") +
  scale_x_continuous(expand = c(0.01, 0)) + # Make the space between the labels and plot smaller
  scale_y_discrete(expand = c(0.01, 0)) + # Make it so the top series actually fits in the plot
  ylab(NULL) + xlab("Number of prey MOTUs") +
  facet_wrap(~Species, ncol = 5) + # free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(
    strip.background = element_rect(fill = "white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"), # strip stuff sorts the facet labels, spacing adjusts the space between facets
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 12)
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "italic"),
    legend.position = "bottom"
  )+ # This makes the facet titles italic
  theme(axis.text.y = element_text(vjust=-2.5)) # Make the y axis labels a little higher

sp_ridge

pdf("plots/degree_ridges.pdf", width = 12)
sp_ridge
dev.off()

jpeg("plots/degree_ridges.jpg", units = "in", width = 12, height = 7, res = 500)
sp_ridge
dev.off()

mean_degree <- degree_df %>% 
  group_by(Species, Site) %>% 
  summarise(m = mean(degree))

tall <- ggplot(degree_df, aes(y = fct_rev(Site), x = degree, fill = fct_rev(hab_type))) +
  geom_density_ridges(scale = 0.85) + # The scale determines the space between the rows
  theme_ridges() + # This changes the theme to make it more aesthetically pleasing
  scale_fill_cyclical(values = c("#d0ca9f", "#85d7da"), guide = "legend", name = "Habitat type") +
  scale_x_continuous(expand = c(0.01, 0)) + # Make the space between the labels and plot smaller
  scale_y_discrete(expand = c(0.01, 0)) + # Make it so the top series actually fits in the plot
  ylab(NULL) + xlab("Number of prey OTUs") +
  facet_wrap(~Species, ncol = 2) + # free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(
    strip.background = element_rect(fill = "white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"), # strip stuff sorts the facet labels, spacing adjusts the space between facets
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 12)
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "italic"),
    legend.position = "bottom"
  ) # This makes the facet titles italic)

tall

pdf("plots/degree_ridges_tall.pdf", height = 12)
tall
dev.off()

# For presentation

smallersp_ridge <- ggplot(degree_df[-which(!degree_df$Species %in% c(c("Kerivoula intermedia", "Kerivoula papillosa", "Rhinolophus borneensis", "Rhinolophus sedulus", "Rhinolophus trifoliatus"))), ], aes(y = fct_rev(Site), x = degree, fill = fct_rev(hab_type))) +
  geom_density_ridges(scale = 0.85) + # The scale determines the space between the rows
  theme_ridges() + # This changes the theme to make it more aesthetically pleasing
  scale_fill_cyclical(values = c("#d0ca9f", "#85d7da"), guide = "legend", name = "Habitat type") +
  scale_x_continuous(expand = c(0.01, 0)) + # Make the space between the labels and plot smaller
  scale_y_discrete(expand = c(0.01, 0)) + # Make it so the top series actually fits in the plot
  ylab("Number of bats, Site") + xlab("Number of prey OTUs") +
  facet_wrap(~Species, ncol = 5) + # free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(
    strip.background = element_rect(fill = "white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"), # strip stuff sorts the facet labels, spacing adjusts the space between facets
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 12)
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "italic"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(colour = , size = 15),
    legend.title = element_text(size = 15, face = "bold")
  )

smallersp_ridge

pdf("plots/presentation_degree_ridges.pdf", width = 13)
smallersp_ridge
dev.off()



#### for full corrplot ####

library(corrplot)

taxa_for_cor <- t(taxa_mat)
cormat <- round(cor(taxa_for_cor), 2)
res1 <- cor.mtest(taxa_mat, conf.level = .95)

taxa_corr <- corrplot(cormat,
  method = "circle", p.mat = res1$p, sig.level = .05, type = "upper", order = "AOE",
  tl.col = "black", tl.srt = 45, insig = "blank",
  bg = "black"
)

#### manipulate the degree for corplot ####


n_matches <- 5
for_bigmat <- t(taxa_mat[, which(colSums(taxa_mat) > n_matches)])
for_bigmat <- for_bigmat[, -which(colSums(for_bigmat) <= 20)]


big_cor <- for_bigmat
bigcormat <- round(cor(big_cor), 2)
resbig <- cor.mtest(for_bigmat)



pdf("plots/order_correlations.pdf")

corrplot(bigcormat,
  method = "circle", p.mat = res1$p, sig.level = .05, type = "upper", order = "AOE",
  tl.col = "black", tl.srt = 45, insig = "blank", col = c("black", "white"),
  bg = "lightblue",
  cl.pos = "b"
)


dev.off()
