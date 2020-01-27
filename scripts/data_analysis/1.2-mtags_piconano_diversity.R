library(tidyverse)
library(ggpubr)
library(stringr)
library(vegan)
library(data.table)
library(indicspecies)
theme_set(theme_bw(base_size = 12))

# Auxiliary functions

str_pad_custom <- function(labels){
  new_labels <- stringr::str_pad(labels, 20, "right")
  return(new_labels)
}

standard_error <- function(x) sd(x)/sqrt(length(x))

# Load data

mtags <- 'data/tables/MPN_VP_metaG_OTUtable.txt'
mtags_table_wide <- read_tsv(mtags)

metadata <- 'data/metadata/MPN_VP_metadata.txt'
metadata_table <- read_tsv(metadata)

# Reformat to long format and add metadata

mtags_table_md <- 
  mtags_table_wide %>% 
  gather(key = Sample,
         value = Abundance,
         contains('MP')) %>% 
  left_join(metadata_table, by = 'Sample')

# Diversity

## Vertical distribution along the water column

mtags_table_relab_group <- 
  mtags_table_md %>% 
  filter(!is.na(Group)) %>% 
  group_by_at(setdiff(colnames(.),
                      c('OTU','Abundance'))) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  group_by(Sample) %>% 
  mutate(rel_abun = 100*Abundance/sum(Abundance)) %>% 
  mutate(Depth_2 = ifelse(Depth <= 200, 
                          Depth, 
                          ifelse(Depth > 200 & Depth <= 1000, 
                                 200 + (Depth-200)/4, 
                                 400 + (Depth-1000)/15)),
         Ocean = case_when(Ocean %in% c('SAB','IO') ~ 'Indian',
                           Ocean == 'AO' ~ 'Atlantic',
                           Ocean == 'PO' ~ 'Pacific'))

groups_to_represent <- 25

vertical_abundant_pico <- 
  mtags_table_relab_group %>%
  filter(!grepl('InSed', Group),
         Fraction == 'Pico') %>% 
  group_by(Group) %>% 
  summarise(total = sum(Abundance)) %>% 
  arrange(desc(total)) %>% 
  .$Group %>% 
  head(n = groups_to_represent)

mtags_table_relab_group_pico <- 
  mtags_table_relab_group %>% 
  filter(Fraction == 'Pico',
         Group %in% vertical_abundant_pico)

mtags_table_relab_group_pico$Group <- 
  factor(mtags_table_relab_group_pico$Group,
         levels = vertical_abundant_pico)

vertical_abundant_nano <- 
  mtags_table_relab_group %>%
  filter(!grepl('InSed', Group),
         Fraction == 'Nano') %>% 
  group_by(Group) %>% 
  summarise(total = sum(Abundance)) %>% 
  arrange(desc(total)) %>% 
  .$Group %>% 
  head(n = groups_to_represent)

mtags_table_relab_group_nano <- 
  mtags_table_relab_group %>% 
  filter(Fraction == 'Nano',
         Group %in% vertical_abundant_nano)

mtags_table_relab_group_nano$Group <- 
  factor(mtags_table_relab_group_nano$Group,
         levels = vertical_abundant_nano)

figureS5A <- 
  ggplot(data = mtags_table_relab_group_pico,
         aes(y = rel_abun, 
             x = Depth_2)) +
  geom_point(aes(color = Ocean), 
             size = 1.5, 
             alpha = 0.85) +
  scale_color_manual(values = c('#a36249', '#a5bb4e','#00abad')) +
  facet_wrap(~Group, 
             scales = "free_x", 
             nrow = 5) +
  geom_smooth(method = "lm", 
              formula =  y ~ poly(x, 4), 
              se = FALSE, 
              color = "goldenrod1", 
              size = 1) +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        panel.background = element_blank()) + 
  scale_x_continuous(breaks = c(0, 200, 400, 600), 
                     labels = c(0, 200, 1000, 4000), 
                     trans = 'reverse') +
  coord_flip() +
  ylab('Relative abundance (%)') + 
  xlab('Depth (m)') +
  ggtitle('Pico fraction')

figureS5B <- 
  ggplot(data = mtags_table_relab_group_nano,
         aes(y = rel_abun, 
             x = Depth_2)) +
  geom_point(aes(color = Ocean), 
             size = 1.5, 
             alpha = 0.85) +
  scale_color_manual(values = c('#a36249', '#a5bb4e','#00abad')) +
  facet_wrap(~Group, 
             scales = "free_x", 
             nrow = 5) +
  geom_smooth(method = "lm", 
              formula =  y ~ poly(x, 4), 
              se = FALSE, 
              color = "goldenrod1", 
              size = 1) +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        panel.background = element_blank()) + 
  scale_x_continuous(breaks = c(0, 200, 400, 600), 
                     labels = c(0, 200, 1000, 4000), 
                     trans = 'reverse') +
  coord_flip() +
  ylab('Relative abundance (%)') + 
  xlab('Depth (m)') +
  ggtitle('Nano fraction')

figureS5 <- 
  ggarrange(figureS5A,
            figureS5B,
            common.legend = TRUE,
            legend = 'top',
            labels = c('A','B'),
            nrow = 2)

ggsave(plot = figureS5,
       filename = 'figureS5.pdf',
       height = 350,
       width = 200,
       units = 'mm')

## Vertical distribution by layer and fraction

mtags_table_relab_group_lf <- 
  mtags_table_relab_group %>% 
  mutate(Layer_Fraction = paste(Fraction, Layer, sep = '-')) %>% # add column merging Layer and Fraction
  filter(!grepl(x = Group, pattern = 'InSed'))
  
abundant_groups_median <- 
  mtags_table_relab_group_lf %>% 
  group_by(Group, Layer_Fraction) %>% 
  summarise(mediana = median(rel_abun)) %>% 
  arrange(desc(mediana)) %>% 
  .$Group %>% 
  unique()

mtags_table_relab_group_lf$Layer_Fraction <- 
  factor(mtags_table_relab_group_lf$Layer_Fraction,
         levels = c('Pico-Photic',
                    'Pico-Aphotic',
                    'Nano-Photic',
                    'Nano-Aphotic'))

mtags_table_relab_group_lf$Fraction <- 
  factor(mtags_table_relab_group_lf$Fraction,
         levels = c('Pico','Nano'))

mtags_table_relab_group_lf$Group <- 
  factor(mtags_table_relab_group_lf$Group,
         levels = abundant_groups_median)

f_groups_boxplot <- function(vector_groups){
  plot <- 
    ggplot(data = mtags_table_relab_group_lf %>% filter(Group %in% abundant_groups_median[vector_groups]),
           aes(x = Group,
               y = rel_abun,
               fill = Layer_Fraction)) +
    geom_blank() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 10, fill = "gray90", color = NA) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.01, ymax = 0.1, fill = "gray90", color = NA) +
    geom_boxplot(outlier.size = 1,
                 outlier.alpha = 0.5,
                 size = 0.4) +
    scale_fill_manual(values = c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4"), # first 4 values from 'Paired' palette reordered
                      labels = str_pad_custom) + # function to add more whitespace between legend items
    ylab('Relative abundance (%)') +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~Fraction, 
               ncol = 1,
               strip.position = 'right') +
    theme(strip.background = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          panel.grid = element_blank()) +
    scale_y_log10(breaks = c(100,10,1,0.1,0.01),
                  labels = c(100,10,1,0.1,0.01)) +
    geom_vline(xintercept = 1.5:9.5, linetype = 'dotted') +
    geom_blank(aes(y = 100))
  return(plot)
}

p_groups_boxplot_1 <- f_groups_boxplot(1:10)
p_groups_boxplot_2 <- f_groups_boxplot(11:20)
p_groups_boxplot_3 <- f_groups_boxplot(21:30)

figure3 <- 
  ggarrange(p_groups_boxplot_1,
            p_groups_boxplot_2,
            p_groups_boxplot_3,
            nrow = 3,
            common.legend = TRUE)

ggsave(plot = figure3,
       filename = 'figure3.pdf',
       width = 170,
       height = 270,
       units = 'mm') # then manually add 'P' and 'A' legends

## Summary of vertical distribution by layer and fraction

mtags_table_median_lf <- 
  mtags_table_relab_group_lf %>% 
  group_by(Group, Layer_Fraction) %>% 
  summarise(median_rel_abun = median(rel_abun),
            mean = mean(rel_abun),
            std_err = standard_error(rel_abun)) %>% 
  ungroup() %>% 
  mutate(Group = as.character(Group)) %>% 
  filter(Group %in% abundant_groups_median[1:30]) %>% 
  group_by(Group) %>% 
  mutate(perc = 100*median_rel_abun/sum(median_rel_abun))

mtags_table_median_lf$Layer_Fraction <- 
  factor(mtags_table_median_lf$Layer_Fraction,
         levels = c('Pico-Photic',
                    'Pico-Aphotic',
                    'Nano-Photic',
                    'Nano-Aphotic'))

mtags_table_higher_median_lf <- 
  mtags_table_median_lf %>% 
  filter(perc == max(perc)) %>% 
  mutate(higher_layer = Layer_Fraction,
         higher_perc = perc) %>% 
  select(Group, higher_layer, higher_perc)

mtags_table_merge_median_lf <- 
  left_join(mtags_table_median_lf,
            mtags_table_higher_median_lf,
            by = 'Group') %>% 
  arrange(desc(higher_layer), higher_perc)

mtags_table_merge_median_lf$Group <-  
  factor(mtags_table_merge_median_lf$Group,
         levels = unique(mtags_table_merge_median_lf$Group))

figure4 <- 
  ggplot(data = mtags_table_merge_median_lf,
         aes(x = Layer_Fraction,
             y = Group)) +
  geom_point(aes(size = perc,
                 fill = higher_layer),
             shape = 21,
             color = 'black',
             stroke = 0.3) +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        plot.margin = margin(t = 0, r = 15, b = 15, l = 0, unit = "mm"),
        legend.position = c(.2,-.05),
        legend.direction = 'horizontal',
        legend.text = element_text(size = 8)) +
  scale_radius(range = c(0,6), 
               breaks = seq(0,100,25), 
               labels = paste0(seq(0,100,25),'%'),
               limits = c(0,100)) + 
  scale_x_discrete(position = "top") +
  guides(fill = FALSE, 
         size = guide_legend(title = ''))

ggsave(plot = figure4,
       filename = 'figure4.pdf',
       height = 200,
       width = 80,
       units = 'mm')

## table S4

tableS4 <- 
  mtags_table_median_lf %>% 
  select(Group, median_rel_abun, mean, std_err, Layer_Fraction) %>% 
  ungroup() %>% 
  mutate(median_rel_abun = round(median_rel_abun, 1),
         mean = round(mean, 1),
         std_err = round(std_err, 1),
         mean_stderr = paste0(mean," ± ",std_err),
         Group = factor(Group, levels = abundant_groups_median[1:30])) %>%
  select(-mean, -std_err) %>%
  unite(temp, median_rel_abun, mean_stderr, sep = '/') %>% 
  spread(key = Layer_Fraction, value = temp) %>% 
  arrange(Group)

write_tsv(tableS4,
          'tableS4.txt') # put the table pretty in Excel

## NMDS at OTU level

shared_samples_df <- 
  inner_join(metadata_table %>% select(Fraction,Sample,Station,Depth,Approach) %>% unique %>% filter(Approach == 'mTags', Fraction == 'Pico'),
             metadata_table %>% select(Fraction,Sample,Station,Depth,Approach) %>% unique %>% filter(Approach == 'mTags', Fraction == 'Nano'),
             by = c('Station','Depth'))

shared_samples <- 
  c(shared_samples_df$Sample.x, shared_samples_df$Sample.y)

mtags_table_relab_otu_shared <- 
  mtags_table_md %>% 
  filter(!is.na(OTU),
         Sample %in% shared_samples) %>% 
  group_by(Sample, OTU) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  mutate(rel_abun = 100*Abundance/sum(Abundance)) %>% 
  select(OTU, Sample, rel_abun) %>% 
  spread(key = OTU, value = rel_abun, fill = 0) %>% 
  column_to_rownames('Sample')

bc_mtags_otu <- 
  vegdist(mtags_table_relab_otu_shared, 
          method = "bray")

nmds_mtags_otu <- 
  metaMDS(bc_mtags_otu, 
          autotransform = FALSE, 
          distance = 'bray', 
          trymax = 100)

nmds_mtags_otu_table_md <- 
  as.data.frame(scores(nmds_mtags_otu)) %>% 
  rownames_to_column(var = 'Sample') %>% 
  as_tibble() %>% 
  left_join(metadata_table, by = 'Sample') %>% 
  mutate(Fraction = factor(Fraction, levels = c('Pico','Nano')),
         Layer = factor(Layer, levels = c('Photic','Aphotic')))

figureS6 <- 
  ggplot(data = nmds_mtags_otu_table_md, mapping = aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = Fraction, color = Layer), size = 3, alpha = 0.7) +
  scale_color_manual(values = c('#fca311', '#1546af'))

ggsave(plot = figureS6,
       filename = 'figureS6.pdf',
       height = 80,
       width= 115,
       units = 'mm')

# PERMANOVA

metadata_for_permanova <- 
  metadata_table %>% 
  filter(Sample %in% rownames(mtags_table_relab_otu_shared))

adonis2(formula = bc_mtags_otu ~ Fraction + Layer2 + Ocean + temp + salinity + O2 + cond,
        data = metadata_for_permanova,
        permutations = 999, 
        method = "bray")
