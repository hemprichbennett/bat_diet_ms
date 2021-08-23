# 6 traps, 5 nights, 3 sites, 2 rounds, 3 years
safe <- 6  * 5 * 3 *2 * 3

danum <- 6 * 10 * 1 * 1 * 2

maliau <- 6 * 10 * 1 * 1 * 2

safe + danum + maliau


library(tidyverse)
field_data <- read_csv('data/Edited_all_files_with_lat_long_VKedits.csv')



focal_sp <- c('Hice',
              'Hidi',
              'Hidy',
              'Hiri',
              'Keha',
              'Kein',
              'Kemi',
              'Kepa',
              'Rhbo',
              'Rhse',
              'Rhtr')


temp <- field_data %>% 
  filter(Site %in% c('SAFE', 'DVCA', "MALUA", "MALIAU", "DANUM")) %>%
  filter(!Block %in% c('D', 'E', 'F')) %>%
  unite(LatLong, Lat, Long) 


length(unique(temp$Latlong))

all_biopsies <- field_data %>%
  filter(!is.na(Biopsy_Dave),
         Biopsy_Dave != 0)

all_biopsies %>% filter(!Species %in% focal_sp)

