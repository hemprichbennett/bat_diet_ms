#### Header ####
## Project: bat-diet
## Script purpose: visualising the family-level diets of each bat species
## Date: 20/07/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(reshape2)
  library(corrplot)
  library(DataExplorer)
  library(dplyr)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')


#####Data import and cleaning #####
field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data_2 <- field_data

field_data$Faeces_no2 <- NULL

field_data_2$Faeces_no1 <- field_data_2$Faeces_no2
field_data_2$Faeces_no2 <- NULL

field_data <- rbind(field_data, field_data_2) #We do this so the merge can work for both columns later


field_data$Species <- gsub('Hebl', 'Hesperoptenus blanfordi', field_data$Species)
field_data$Species <- gsub('Hice', 'Hipposideros cervinus', field_data$Species)
field_data$Species <- gsub('Hidi', 'Hipposideros diadema', field_data$Species)
field_data$Species <- gsub('Hidy', 'Hipposideros dyacorum', field_data$Species)
field_data$Species <- gsub('Hiri', 'Hipposideros ridleyi', field_data$Species)
field_data$Species <- gsub('Keha', 'Kerivoula hardwickii', field_data$Species)
field_data$Species <- gsub('Kein', 'Kerivoula intermedia', field_data$Species)
field_data$Species <- gsub('Kepa', 'Kerivoula papillosa', field_data$Species)
field_data$Species <- gsub('Rhbo', 'Rhinolophus borneensis', field_data$Species)
field_data$Species <- gsub('Rhse', 'Rhinolophus sedulus', field_data$Species)
field_data$Species <- gsub('Rhtr', 'Rhinolophus trifoliatus', field_data$Species)
field_data$Species <- gsub('Nytr', 'Nycteris tragata', field_data$Species)
field_data$Species <- gsub('Muro', 'Murina rozendaali', field_data$Species)
field_data$Species <- gsub('Kepe', 'Kerivoula pellucida', field_data$Species)

#Import and sort the interaction matrix
all_interactions <- read.csv('data/processed_dna_data/lulu/95/lulu_95.csv', header = F, stringsAsFactors = F)#
all_interactions[1,1] <- 'MOTU'
all_interactions[1,] <- gsub('X', '', all_interactions[1,])
all_interactions[2:nrow(all_interactions),2:ncol(all_interactions)] <- ifelse(all_interactions[2:nrow(all_interactions), 2:ncol(all_interactions)]==0,0,1)





#This finds all samples with 'GC' in the name and gives them a useful name
gc <- grep('GC',all_interactions[1,])
for(i in 1:length(gc)){
  #print(str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7])
  temp <- str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7]
  all_interactions[1,gc[i]] <- str_split(temp, pattern='_')[[1]][1]
}

#This finds all samples without 'GC' in the name and gives them a useful name
non_gc <- seq(1, ncol(all_interactions))[-gc]
str_split(all_interactions[1,non_gc[2]], pattern = '\\.')[[1]][1] #Finish this
for(i in 1:length(non_gc)){
  all_interactions[1,non_gc[i]] <- str_split(all_interactions[1,non_gc[i]], pattern = '\\.')[[1]][1]
}

badcols <- c('1774','4437', '2070', '2275', '4260', '4531', "1004", "1007", "1107", "1134", "1165", "1180",  "198",  "209",  "210","387",  "426",  "459",  "497",  "536",  "541",  "567",  "591",
             "689","796",  "806",  "822",  "841",  "843",  "899",  "910",  "918",  "986","996", "3712", "4341", "4361",'1774','4437', '2070', '2275', '4260', '4531', '841', '843')#Sadly these columns match two different samples, so must be removed for now until checked against the field data

all_interactions <- all_interactions[,-which(all_interactions[1,] %in% badcols)]

for_rownames <- all_interactions[,1]
all_interactions <- all_interactions[,-1]

all_interactions <- apply(all_interactions, 2, as.numeric)
rownames(all_interactions) <- for_rownames #This is a roundabout way of renaming the rownames, but if we do it in the normal manner before the apply, the apply removes the names

#Import the prey data

prey_data <- read.csv('data/taxonomy/family.csv')
colnames(prey_data) <- c('MOTU', 'Taxa')
prey_data$MOTU <- as.character(prey_data$MOTU)




###Add the taxonomic information to the interactions for everything
taxa_mat <- all_interactions[1,]
z <- 2
for(i in 1: nrow(all_interactions)){
  rowname = rownames(all_interactions)[i]
  if(!rowname %in% prey_data$MOTU){
    next()
  }
  tax = as.character(prey_data[which(prey_data$MOTU == rowname),'Taxa'])
  if(is.null(nrow(taxa_mat))){ #If its the first iteration there won't be any rownames yet, so the next if statement will fail
    taxa_mat <- rbind(taxa_mat, all_interactions[i,])
    rownames(taxa_mat)[z] <- tax
    z <- z+1
  }
  if(tax %in% rownames(taxa_mat)){
    to_merge = which(rownames(taxa_mat)==tax)
    taxa_mat[to_merge,] <- taxa_mat[to_merge,]+ all_interactions[i,]
  }else{
    taxa_mat <- rbind(taxa_mat, all_interactions[i,])
    rownames(taxa_mat)[z] <- tax
    z <- z+1
  }
}
#The colnames currently have the sample numbers, we want them to be a useable value in the data frame
colnames(taxa_mat) <- taxa_mat[1,]
taxa_mat <- taxa_mat[-1,]

t_taxa_mat <- data.frame(t(taxa_mat))

t_taxa_mat$Sample <- colnames(taxa_mat)

t_taxa_mat <- merge(x=t_taxa_mat, y =field_data, by.x = 'Sample', by.y = 'Faeces_no1')



tax_df <- t_taxa_mat[,c(seq(2,146),164)]
tax_df <- melt(tax_df, id.vars = c('Species'))
colnames(tax_df)[c(2,3)] <- c('Family', 'Present/absent')
tax_df$`Present/absent` <- as.integer(as.character(tax_df$`Present/absent`))
tax_df$`Present/absent` <- ifelse(tax_df$`Present/absent`== 0, 0, 1)

prop_present <- sapply(unique(tax_df[,c('Species', 'Family')]),  function(x) as.character(x))
prop <- c()
nbats <- c()
#for(i in 1: 1){
for(i in 1: nrow(prop_present)){
  tem <- tax_df[which(tax_df$Species== as.character(prop_present[i,1])
                      & tax_df$Family==as.character(prop_present[i,2])),]
  prop <- c(prop, sum(tem$`Present/absent`)/nrow(tem)) #The number of bats that consumed the Family, divided by total bats
  nbats <- c(nbats, nrow(tem))
}


prop_present <- cbind(prop_present, prop, nbats)
prop_present <- as.data.frame(prop_present)
prop_present$prop <- as.numeric(as.character(prop_present$prop))
prop_present$nbats <- as.integer(as.character(prop_present$nbats))

prop_present$for_x <- paste(prop_present$Species, ' (', prop_present$nbats, ' samples)', sep = '')
prop_present$for_x <- gsub('\\(1 samples\\)', '\\(1 sample\\)', prop_present$for_x)

prop_present$order <- rep(NA, nrow(prop_present)) #Make a column of the order which each family belongs to
prop_present$order[prop_present$Family %in% c('Sparassidae', 'Salticidae', 'Araneidae', 'Theridiidae', 'Clubionidae', 'Linyphiidae', 'Pholcidae', 
                                              'Uloboridae')] <- 'Araneae (Spiders)'
prop_present$order[prop_present$Family %in% c('Ectobiidae', 'Blaberidae', 'Termitidae', 'Blattidae', 'Rhinotermitidae', 'Cryptocercidae')] <- 'Blattodea (Termites\nand cockroaches)'
prop_present$order[prop_present$Family %in% c('Chrysomelidae', 'Cerambycidae', 'Curculionidae', 'Elateridae', 'Carabidae', 'Ptilodactylidae', 'Mordellidae', 'Eucnemidae',
                                              'Latridiidae', 'Dermestidae', 'Mycetophagidae', 'Dytiscidae', 'Nitidulidae', 'Cleridae', 'Leiodidae', 'Byrrhidae',
                                              'Throscidae', 'Scarabaeidae', 'Psephenidae', 'Tenebrionidae', 'Zopheridae', 'Cantharidae')] <- 'Coleoptera (Beetles)'
prop_present$order[prop_present$Family %in% c('Culicidae', 'Chironomidae', 'Tachinidae', 'Sciaridae', 'Psychodidae', 'Phoridae', 'Ceratopogonidae', 
                                              'Muscidae', 'Mycetophilidae', 'Chloropidae', 'Calliphoridae', 'Tabanidae', 'Stratiomyidae', 'Cecidomyiidae', 'Milichiidae',
                                              'Syrphidae', 'Asilidae', 'Ephydridae', 'Dolichopodidae', 'Pediciidae', 'Drosophilidae', 'Keroplatidae', 'Tipulidae',
                                              'Platypezidae')] <- 'Diptera (Flies)'
prop_present$order[prop_present$Family %in% c('Entomobryidae', 'Isotomidae')] <- 'Entomobryomorpha (Springtails)'
prop_present$order[prop_present$Family %in% c('Baetidae', 'Heptageniidae', 'Ephemerellidae')] <- 'Ephemeroptera (Mayflies)'
prop_present$order[prop_present$Family %in% c('Euphausiidae')] <- 'Euphausiacea (krill)'
prop_present$order[prop_present$Family %in% c('Armadillidiidae')] <- 'Isopoda (Woodlice)'
prop_present$order[prop_present$Family %in% c('Cicadidae', 'Cicadellidae', 'Pentatomidae', 'Derbidae', 'Aphididae', 'Miridae', 'Fulgoridae', 'Flatidae', 'Cydnidae',
                                              'Pemphigidae', 'Kinnaridae', 'Rhyparochromidae', 'Dictyopharidae', 'Delphacidae', 'Nogodinidae', 'Cixiidae', 'Hormaphididae',
                                              'Lygaeidae')] <- 'Hemiptera (True bugs)'
prop_present$order[prop_present$Family %in% c('Ichneumonidae', 'Mutillidae', 'Perilampidae', 'Formicidae', 'Braconidae', 'Dryinidae', 'Agaonidae', 'Vespidae',
                                              'Tenthredinidae')] <- 'Hymenoptera (Ants, wasps, etc)'
prop_present$order[prop_present$Family %in% c('Geometridae', 'Lymantriidae', 'Tortricidae', 'Noctuidae', 'Lycaenidae', 'Sphingidae', 'Crambidae', 'Tineidae', 'Erebidae',
                                              'Nolidae', 'Adelidae', 'Pterophoridae', 'Limacodidae', 'Pyralidae', 'Blastobasidae', 'Hesperiidae', 'Lecithoceridae')] <- 'Lepidoptera (Butterflies and moths)'
prop_present$order[prop_present$Family %in% c('Mantidae', 'Tarachodidae')] <- 'Mantodea (Mantises)'
prop_present$order[prop_present$Family %in% c('Laelapidae')] <- 'Mesostigmata (Mites)'
prop_present$order[prop_present$Family %in% c('Chrysopidae', 'Hemerobiidae', 'Coniopterygidae')] <- 'Neuroptera (Net-winged insects)'
prop_present$order[prop_present$Family %in% c('Brachychthoniidae')] <- 'Oribatida (Moss mites)'
prop_present$order[prop_present$Family %in% c('Acrididae', 'Gryllidae', 'Tettigoniidae', 'Acrididae', 'Chorotypidae')] <- 'Orthoptera (Grasshoppers)'
prop_present$order[prop_present$Family %in% c('Diapheromeridae', 'Phasmatidae')] <- 'Phasmatodea (Stick insects)'
prop_present$order[prop_present$Family %in% c('Paradoxosomatidae')] <- 'Polydesmida (Millipedes)'
prop_present$order[prop_present$Family %in% c('Caeciliusidae', 'Archipsocidae', 'Stenopsocidae', 'Lepidopsocidae', 'Psocidae', 'Epipsocidae')] <- 'Psocoptera (Barkflies)'
prop_present$order[prop_present$Family %in% c('Eremaeidae', 'Brachychthoniida', 'Scheloribatidae', 'Ceratozetidae', 'Oppiidae')] <- 'Sarcoptiformes (Mites)'
prop_present$order[prop_present$Family %in% c('Philopotamidae', 'Psychomyiidae', 'Leptoceridae')] <- 'Trichoptera (Caddisflies)'
prop_present$order[prop_present$Family %in% c('Phlaeothripidae')] <- 'Thysanoptera (Thrips)'
prop_present$order[prop_present$Family %in% c('Tarsonemidae', 'Eupodidae', 'Bdellidae', 'Cheyletidae', 'Ereynetidae', 'Teneriffiidae', 'Cunaxidae')] <- 'Trombidiformes (Mites)'



prop_present$order <-  gsub(' .+', '', prop_present$order)

write.csv(prop_present, 'results/family_proportions.csv')

longtiles <- ggplot(data = prop_present[which(prop_present$prop!=0),], 
                    aes(y = fct_rev(Family), x = for_x)) + 
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.5, 
                       name ='Proportion of\nsamples containing', limits = c(0,1)) +
  labs(x ="Bat species", y = 'Prey')+
  theme_linedraw()+
  theme(strip.background =element_rect(fill="black", size = 3))+
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(colour = 'white', size = 8, angle = 180), 
        strip.placement = "outside",
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray",
                                        size = 0.5, linetype = NULL))+
  facet_grid(order ~ ., scales = 'free', space="free", switch = "y")

longtiles

pdf('plots/familyplots/tile_long_family.pdf', height = 14)
longtiles
dev.off()

# Now make a series of these family plots over multiple figures, it'll be
# easier to read than squidging it all into a single one


# Orders A-D --------------------------------------------------------------

prop_1 <- prop_present %>%
  filter(order %in% unique(prop_present$order)[grep('^[A-D].+', unique(prop_present$order))]) %>%
  filter(prop != 0)

long_1 <- ggplot(data = prop_1, 
                 aes(y = fct_rev(Family), x = for_x)) + 
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.5, 
                       name ='Proportion of\nsamples containing', limits = c(0,1)) +
  labs(x ="Bat species", y = 'Prey')+
  theme_linedraw()+
  theme(strip.background =element_rect(fill="black", size = 3))+
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(colour = 'white', size = 8, angle = 180), 
        strip.placement = "outside",
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray",
                                        size = 0.5, linetype = NULL))+
  facet_grid(order ~ ., scales = 'free', space="free", switch = "y")

long_1
ggsave('plots/familyplots/long_1.pdf', long_1, height = 14)


# Orders E-L --------------------------------------------------------------

prop_2 <- prop_present %>%
  filter(order %in% unique(prop_present$order)[grep('^[E-L].+', unique(prop_present$order))]) %>%
  filter(prop != 0)

long_2 <- ggplot(data = prop_2, 
                 aes(y = fct_rev(Family), x = for_x)) + 
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.5, 
                       name ='Proportion of\nsamples containing', limits = c(0,1)) +
  labs(x ="Bat species", y = 'Prey')+
  theme_linedraw()+
  theme(strip.background =element_rect(fill="black", size = 3))+
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(colour = 'white', size = 8, angle = 180), 
        strip.placement = "outside",
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray",
                                        size = 0.5, linetype = NULL))+
  facet_grid(order ~ ., scales = 'free', space="free", switch = "y")

long_2
ggsave('plots/familyplots/long_2.pdf', long_2, height =14)


# Orders M-Z --------------------------------------------------------------

prop_3 <- prop_present %>%
  filter(order %in% unique(prop_present$order)[grep('^[M-Z].+', unique(prop_present$order))]) %>%
  filter(prop != 0)

long_3 <- ggplot(data = prop_3, 
                 aes(y = fct_rev(Family), x = for_x)) + 
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.5, 
                       name ='Proportion of\nsamples containing', limits = c(0,1)) +
  labs(x ="Bat species", y = 'Prey')+
  theme_linedraw()+
  theme(strip.background =element_rect(fill="black", size = 3))+
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(colour = 'white', size = 8, angle = 180), 
        strip.placement = "outside",
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray",
                                        size = 0.5, linetype = NULL))+
  facet_grid(order ~ ., scales = 'free', space="free", switch = "y")

long_3
ggsave('plots/familyplots/long_3.pdf', long_3, height = 14)



widetiles <- ggplot(data = prop_present[which(prop_present$prop!=0),], aes(x = Family, y = fct_rev(for_x))) + 
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.5, name ='Proportion of\nsamples\ncontaining', limits = c(0,1)) +
  labs(y ="Bat species", x = 'Prey')+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(colour = 'white', size = 8, angle = 90), strip.placement = "outside",
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray",
                                        size = 0.5, linetype = NULL),
        strip.background =element_rect(fill="black", size = 3))+
  facet_grid(. ~order, scales = 'free', space="free", switch = 'both')#+
  
widetiles

pdf('plots/familyplots/tile_wide_family.pdf', width = 13)
widetiles
dev.off()


balloons <- ggplot(data = prop_present[which(prop_present$prop!=0),], aes(x = Family, y =fct_rev(for_x))) + geom_point(aes(size=prop))+ 
  labs(y ="Bat species", x = 'Prey', size = 'Proportion of\nsamples\ncontaining')+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(colour = 'white', size = 8, angle = 90), strip.placement = "outside",
        strip.background =element_rect(fill="black", size = 3))+
  facet_grid(. ~order, scales = 'free', space="free", switch = 'both')#+
balloons


pdf('plots/familyplots/balloon_wide_family.pdf', width = 16)
balloons
dev.off()

