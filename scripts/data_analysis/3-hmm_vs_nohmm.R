library(tidyverse)
theme_set(theme_bw(base_size = 12))

# Auxiliary functions

source('scripts/data_analysis/stat_smooth_func.R')

# Load data

hmm_table_wide <- 
  read_tsv('data/tables/MPN_hmm_comparison_table.txt')

# Reformat

hmm_table <- 
  hmm_table_wide %>% 
  gather(key = Sample,
         value = Abundance,
         contains('_')) %>% 
  separate(Sample,
           into = c('Sample','Protocol'),
           sep = '_') %>%
  filter(!is.na(Supergroup),
         Supergroup != 'Eukaryota') %>% 
  group_by(Supergroup, Sample, Protocol) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  spread(key = Protocol,
         value = Abundance,
         fill = 0)

supergroups_by_diff <- 
  hmm_table %>% 
  group_by(Supergroup) %>% 
  summarise(HMM = sum(HMM), 
            noHMM = sum(noHMM), 
            perc = 100*HMM/noHMM) %>% 
  arrange(perc) %>% 
  .$Supergroup

hmm_table$Supergroup <- 
  factor(hmm_table$Supergroup,
        levels = supergroups_by_diff)

# Scatter plots

figureS3 <- 
  ggplot(data = hmm_table, 
         aes(x = noHMM, y = HMM)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = 'lm', 
              color = 'black', 
              size = 0.2,
              se = FALSE) +
  facet_wrap(~Supergroup, scales = 'free', ncol = 3) +
  theme(legend.position = 'none',
        strip.background = element_blank())+
  stat_smooth_func(geom="text", 
                   method="lm", 
                   hjust=0, 
                   parse=TRUE,
                   se = FALSE) +
  xlab('\nNumber of miTags (BLAST)') +
  ylab('Number of miTags (HMM)\n') +
  geom_blank(data = hmm_table %>% 
               group_by(Supergroup) %>% 
               summarise(max = max(c(noHMM,
                                     HMM),
                                   na.rm = TRUE)) %>% 
               mutate(x = ceiling(max), 
                      y = ceiling(max)) %>% 
               select(-max),
             aes(x = x,
                 y = y))

ggsave(plot = figureS3,
       filename = 'figureS3.pdf',
       height = 220,
       width = 220,
       units = 'mm')
