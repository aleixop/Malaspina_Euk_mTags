library(tidyverse)
library(ggpubr)
theme_set(theme_bw(base_size = 12))

# Load data

mtags <- 'data/tables/MPN_VP_metaG_OTUtable.txt'
mtags_table_wide <- read_tsv(mtags)

# Reformat to long format

mtags_table <- 
  mtags_table_wide %>% 
  gather(key = Sample,
         value = Abundance,
         contains('MP')) 

# Taxonomic levels

## overall taxonomic levels 

tax_levels_all <- 
  mtags_table %>% 
  mutate(tax_level = case_when(is.na(Supergroup) ~ 'Ambiguous',
                               is.na(Group) ~ 'Supergroup',
                               is.na(OTU) ~ 'Group',
                               TRUE ~ 'OTU'),
         tax_level = factor(tax_level, levels = c('OTU','Group','Supergroup','Ambiguous'))) %>% 
  group_by(tax_level) %>% 
  summarise(number_mtags = sum(Abundance)) %>% 
  mutate(perc = round(100*number_mtags/sum(number_mtags),digits = 2))

## taxonomic levels by supergroup 

tax_levels_sg <- 
  mtags_table %>% 
  mutate(tax_level = case_when(is.na(Supergroup) ~ 'Ambiguous',
                               is.na(Group) ~ 'Supergroup',
                               is.na(OTU) ~ 'Group',
                               TRUE ~ 'OTU'),
         tax_level = factor(tax_level, levels = c('OTU','Group','Supergroup','Ambiguous'))) %>% 
  group_by(Supergroup, tax_level) %>%
  filter(!is.na(Supergroup),
         Supergroup != 'Eukaryota') %>% 
  summarise(number_mtags = sum(Abundance)) %>% 
  mutate(perc = round(100*number_mtags/sum(number_mtags),digits = 1))

tax_levels_sg$Supergroup <- 
  factor(tax_levels_sg$Supergroup, 
         levels = tax_levels_sg %>% arrange(tax_level, perc) %>% .$Supergroup %>% unique())

## figures

figure1B <- 
  ggplot(data = tax_levels_all,
         aes(x = tax_level,
             y = number_mtags/1000)) + # put to 1000-scale
  geom_col(aes(fill = tax_level), color = 'black', size = 0.15) +
  geom_text(aes(label = paste0(perc,'%'), 
                hjust = ifelse(number_mtags < 50000, -0.1, 1.1), 
                color = ifelse(perc < 10, 'black', 'white'))) +
  scale_fill_manual(name = NULL,
                  breaks = c('OTU','Group','Supergroup'),
                  labels = c(expression(OTU[97 * phantom("+")]),
                             'Group',
                             'Supergroup'),
                  values = c("#35B779FF","#31688EFF","#440154FF","#FDE725FF")) +
  ylab(bquote('Number of mTags (x'*10^3*')')) + 
  xlab(NULL) +
  theme(legend.position = 'top',
        panel.grid = element_blank()) + 
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,180),
                     breaks = seq(0,180,45)) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(tax_levels_all$tax_level)),
                   labels = rev(c(expression(OTU[97 * phantom("+")]),
                                  'Group',
                                  'Supergroup',
                                  'Ambiguous'))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE)

figure1C <- 
  ggplot(data = tax_levels_sg,
         aes(x = Supergroup,
             y = perc)) +
  geom_col(aes(fill = tax_level),
           color = 'black',
           size = 0.15,
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values =c("#35B779FF","#31688EFF","#440154FF","#FDE725FF")) +
  guides(fill = FALSE) +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(5.5,10,5.5,5.5), "pt")) +
  ylab('Percentage of mTags (%)') + 
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  xlab(NULL)

figure1BC <- 
  ggarrange(figure1B, 
            figure1C, 
            labels = c('B','C'), 
            common.legend = TRUE,
            nrow = 2,
            heights = c(1.2,2), 
            legend = 'bottom',
            align = 'v')

ggsave(plot = figure1BC,
       filename = 'figure1BC.pdf',
       height = 120,
       width = 120,
       units = 'mm')
