library(ShortRead)
library(dplyr)

paths <- c('~/Dropbox/work/projects/in_progress/Bat-diet/data/raw_dna_data/GC-EC-7330-66241242/FASTQ_Generation_2018-02-28_16_59_54Z-82611867/',
           '~/Dropbox/work/DNA_sequence_files_general usage/Dave_summer_17_sequencing/')

all_fastqs <- list.files(path = paths, 
           pattern = '.gz', recursive = T, full.names = T)



out_list <- list()
for(f_name in all_fastqs){
  out_list[[f_name]] <- length(readFastq(f_name))
}

# need to look at the blanks in here
out_df <- bind_rows(out_list)
#length(readFastq('~/Dropbox/work/projects/in_progress/Bat-diet/data/raw_dna_data/GC-EC-7330-66241242/FASTQ_Generation_2018-02-28_16_59_54Z-82611867/ds.7d6c69b4674945aaa02530101e41f678/GC-EC-7330-Dave-2017-2-TT240_S360_L001_R1_001.fastq.gz'))
better_formatted <- data.frame(file_path = colnames(out_df),
                               nreads = as.numeric(out_df[1,]))

summary_vals <- better_formatted %>%
  mutate(is_blank = grepl('BLANK', file_path),
         file_name = gsub('.+\\/', '', file_path),
         sample_id_etc = gsub('_.+', '', file_name)) %>%
  select(-file_path, -file_name) %>%
  distinct() %>%
  group_by(is_blank) %>%
  summarise(mean = mean(nreads),
            sd = sd(nreads))

