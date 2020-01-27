library(tidyverse)
theme_set(theme_bw(base_size = 12))
library(RColorBrewer)

# Load data

trimming_comparison_wide <- 
  read_tsv('data/tables/MPN_nano_trimming_comparison_table.txt')

# Reformat

trimming_comparison <- 
  trimming_comparison_wide %>% 
  gather(key = Sample,
         value = Abundance,
         contains('_')) %>%
  separate(Sample, 
           into = c('Sample', 'Length'), 
           sep = '_')
  
# Add tax levels

trimming_comparison_taxlev <- 
  trimming_comparison %>% 
  mutate(Length = factor(Length, levels = c('151bp', '101bp')),
         tax_level = case_when(is.na(Supergroup) ~ 'Ambiguous',
                               is.na(Group) ~ 'Supergroup',
                               is.na(OTU) ~ 'Group',
                               TRUE ~ 'OTU'),
         tax_level = factor(tax_level, levels = c('OTU','Group','Supergroup','Ambiguous'))) %>%
  filter(!is.na(Supergroup)) %>% 
  group_by(Supergroup,
           tax_level,
           Length) %>% 
  summarise(Abundance = sum(Abundance)) 
  
overall_taxlev_perc <- 
  trimming_comparison_taxlev %>% 
  group_by(Length, tax_level) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  mutate(perc = 100*Abundance/sum(Abundance))

order_sg <- 
  trimming_comparison_taxlev %>% 
  group_by(Supergroup) %>% 
  summarise(total = sum(Abundance)) %>% 
  arrange(desc(total)) %>% 
  .$Supergroup

trimming_comparison_taxlev$Supergroup <- 
  factor(trimming_comparison_taxlev$Supergroup,
         levels = order_sg)

# Comparison

figureS2 <- 
  ggplot(data = trimming_comparison_taxlev,
         aes(x = tax_level, y = Abundance)) +
    geom_col(aes(fill = Supergroup)) +
    geom_label(data = overall_taxlev_perc,
               aes(x = tax_level,
                   y = Abundance + 6000,
                   label = paste0(round(perc,1),'%'),
                   color = tax_level)) +
    facet_wrap(~Length) +
    scale_fill_brewer(palette = 'Paired') +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 14, hjust = 0)) +
    xlab('Taxonomic level') +
    ylab('Number of mTags') +
    guides(color = FALSE) +
    scale_color_grey(start = 0, 
                     end = 0.5) +
    scale_x_discrete(breaks = c('OTU','Group','Supergroup'),
                     labels = c(expression(OTU[97 * phantom("+")]),
                                'Group',
                                'Supergroup'))

ggsave(plot = figureS2,
       filename =  'figureS2.pdf',
       height = 100,
       width = 200,
       units = 'mm')  

