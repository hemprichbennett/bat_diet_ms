#### Header ####
## Project: Bat diet
## Script purpose: Calculating the (s)hit rate of all bat species
## Date: 22/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################
library(here)
library(magrittr)
#Load in data
field_df <- read.csv(here('data', 'Edited_all_files_with_lat_long_VKedits.csv'))

#Remove pregnant or lactating individuals, as they weren't kept
field_df <- field_df[-which(field_df$Reproductive_condition %in% c('PR', 'P', 'L', 'LA')),]

#get rid of the sites which were never analysed
field_df <- field_df[-which(field_df$Block %in% c('D', 'DV88', 'DV89', 'E',  'F')),]

#Strip it to the required columns
field_df <- field_df[,c('Species', "Faeces_no1", "Faeces_no2")]
field_df$Species <- gsub('Hici', 'Hice', field_df$Species)
field_df$Species <- gsub('Mely', 'Mesp', field_df$Species)
field_df$Species <- gsub('Rbho', 'Rhbo', field_df$Species)
field_df$Species <- gsub('Nyja', 'Nytr', field_df$Species)


#Remove the zero vals
has_samples <- field_df[-which(field_df$Faeces_no1==0 & field_df$Faeces_no2==0),]
has_samples <- has_samples[-which(has_samples$Faeces_no1=='' & has_samples$Faeces_no2==''),]

sampletable <- table(has_samples$Species)

all_species <- table(field_df$Species)

#Make dataframe, calculate rate
out_df <- as.data.frame(cbind(all_species[which(names(all_species) %in% names(sampletable))], sampletable))
out_df$Rate <-out_df[,2]/out_df[,1]
colnames(out_df) <- c('Number captured', 'Number of samples obtained', 'Rate')

#Kill rownames that are either poor identifications or typos
badrows <- c('0', 'Em--', 'Hi--', 'Ke--', 'Hici', 'Mi--', 'Mu--', 'My--', 'Pi--', 'Pip?', 'Rbho', 'Rh--', 'Unknown')
out_df <- out_df[-which(rownames(out_df) %in% badrows),]

# Replace NAs with zeros
out_df[which(is.na(out_df$Rate)), 'Rate'] <- 0

rownames(out_df) %<>%
  gsub('Bama', 'Balionycteris maculata', .)%<>%
  gsub('Emal', 'Emballonura alecto', .)%<>%
  gsub('Emmo', 'Emballonura monticola', .)%<>%
  gsub('Haha', 'Harpiocephalus harpia', .)%<>%
  gsub('Hebl', 'Hesperoptenus blanfordi', .)%<>%
  gsub('Hiat', 'Hipposideros ater', .)%<>%
  gsub('Hibi', 'Hipposideros bicolor', .)%<>%
  gsub('Hice', 'Hipposideros cervinus', .)%<>%
  gsub('Hidi', 'Hipposideros diadema', .)%<>%
  gsub('Hidy', 'Hipposideros dyacorum', .)%<>%
  gsub('Higa', 'Hipposideros galeritus', .)%<>%
  gsub('Hiri', 'Hipposideros ridleyi', .)%<>%
  gsub('Hisa', 'Hipposideros sabanus', .)%<>%
  gsub('Keha', 'Kerivoula hardwickii', .)%<>%
  gsub('Kein', 'Kerivoula intermedia', .)%<>%
  gsub('Kele', 'Kerivoula lenis', .)%<>%
  gsub('Kemi', 'Kerivoula minuta', .)%<>%
  gsub('Kepa', 'Kerivoula papillosa', .)%<>%
  gsub('Kepe', 'Kerivoula pellucida', .)%<>%
  gsub('Kewh', 'Kerivoula whiteheadi', .)%<>%
  gsub('Mami', 'Macroglossus minimus', .)%<>%
  gsub('Mesp', 'Megaderma spasma', .)%<>%
  gsub('Mewe', 'Megaerops wetmorei', .)%<>%
  gsub('Muae', 'Murina aenea', .)%<>%
  gsub('Mucy', 'Murina cyclotis', .)%<>%
  gsub('Mupe', 'Murina peninsularis', .)%<>%
  gsub('Muro', 'Murina rozendaali', .)%<>%
  gsub('Musu', 'Murina suilla', .)%<>%
  gsub('Mymu', 'Myotis muricola', .)%<>%
  gsub('Myri', 'Myotis ridleyi', .)%<>%
  gsub('Nytr', 'Nycteris tragata', .)%<>%  
  gsub('Phat', 'Phoniscus atrox', .)%<>%
  gsub('Phbr', 'Philetor brachypterus', .)%<>%
  gsub('Pija', 'Pipistrellus javanicus', .)%<>%
  gsub('Pite', 'Pipistrellus tenuis', .)%<>%
  gsub('Rhac', 'Rhinolophus acuminatus', .)%<>%
  gsub('Rhaf', 'Rhinolophus affinis', .)%<>%
  gsub('Rhbo', 'Rhinolophus borneensis', .)%<>%
  gsub('Rhcr', 'Rhinolophus creaghi', .)%<>%
  gsub('Rhlu', 'Rhinolophus luctus', .)%<>%
  gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
  gsub('Rhtr', 'Rhinolophus trifoliatus', .)

#Round it to three decimal places
out_df$Rate <- round(out_df$Rate, digits =3)

focal_species <- c('Hipposideros cervinus', 'Hipposideros diadema', 'Hipposideros dyacorum', 'Kerivoula hardwickii', 'Kerivoula intermedia', 'Kerivoula papillosa', 'Rhinolophus borneensis', 'Rhinolophus sedulus', 'Rhinolophus trifoliatus')

focal_df <- out_df[which(rownames(out_df) %in% focal_species),]


write.csv(out_df, here('results', 'all_shit_rate.csv'))
write.csv(focal_df, here('results', 'focal_sp_shit_rate.csv'))
