library(tidyverse)
library(vegan)
library(data.table)
library(ggpubr)
theme_set(theme_bw(base_size = 12))

# Auxiliary functions

str_pad_custom <- function(labels){
  new_labels <- stringr::str_pad(labels, 15, "right")
  return(new_labels)
}

# Load data

amplicons_mtags <- 'data/tables/MPN_amplicons_mtags_comparison_table.txt'
amplicons_mtags_table_wide <- read_tsv(amplicons_mtags)

metadata <- 'data/metadata/MPN_VP_metadata.txt'
metadata_table <- read_tsv(metadata)

# Reformat and add metadata

amp_mtags_table_md <- 
  amplicons_mtags_table_wide %>% 
  gather(key = Sample2,
         value = Abundance,
         3:ncol(.)) %>% 
  separate(Sample2, 
           into = c('Sample',
                    'Approach'),
           sep = '_',
           remove = FALSE) %>% 
  left_join(metadata_table, by = c('Sample','Approach')) %>% 
  mutate(Approach = factor(Approach, levels = c('mTags','ampliconV4','ampliconV9')))

# Comparison 

amp_mtags_table_relab <- 
  amp_mtags_table_md %>% 
  filter(!is.na(Group)) %>% 
  group_by(Sample2) %>%
  mutate(rel_abun_by_sample = 100*Abundance/sum(Abundance)) %>% 
  ungroup()

abundant_groups_forreg_mean <-
  amp_mtags_table_relab %>%
  group_by(Group, Approach) %>%
  summarise(mean = mean(rel_abun_by_sample)) %>% 
  filter(Approach == 'mTags',
         !grepl(x = Group, pattern = 'InSed')) %>%
  arrange(desc(mean)) %>%
  head(n=20)

## Diversity comparison

groups_median <- 
  amp_mtags_table_relab %>% 
  filter(rel_abun_by_sample > 0) %>% 
  group_by(Group, Approach) %>% 
  summarise(median = median(rel_abun_by_sample)) %>% 
  filter(Approach == 'mTags',
         Group %in% abundant_groups_forreg_mean$Group) %>% 
  arrange(desc(median))

missing_data <- # patch so missing groups in certain approaches also appear in boxplots
  tibble(Group = c('Prymnesiophyceae','Diplonemea','Discosea','Kinetoplastida'),
         Layer = c('Photic','Aphotic','Aphotic','Aphotic'),
         Approach = factor(rep('ampliconV4', 4),
                         levels = c('mTags', 'ampliconV4', 'ampliconV9')),
         rel_abun_by_sample = rep(0.01,4))

amp_mtags_table_relab_wmiss <- 
  amp_mtags_table_relab %>% 
  select(Group, Layer, Approach, Station, Depth, rel_abun_by_sample) %>% 
  bind_rows(missing_data)

### Wilcoxon paired test

wilcox_paired <- compare_means(data = amp_mtags_table_relab %>% arrange(Approach, Sample, Group), 
                               formula = rel_abun_by_sample ~ Approach, 
                               group.by = 'Group', 
                               paired = TRUE,
                               method = 'wilcox.test') %>% 
  filter(Group %in% groups_median$Group) %>% 
  mutate(Group = factor(Group,
                        levels = groups_median$Group))

### Boxplots for comparison

amp_mtags_table_relab_all <- 
  amp_mtags_table_relab_wmiss %>% 
    filter(Group %in% groups_median$Group) %>% 
    mutate(Group = factor(Group,
                          levels = groups_median$Group),
           Approach = factor(Approach, levels = c('mTags', 'ampliconV4', 'ampliconV9')))

figure2_alt <- 
  ggplot(data = amp_mtags_table_relab_all,
         aes(x = Group,
             y = rel_abun_by_sample)) +
  geom_blank() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 10, fill = "gray90", color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.01, ymax = 0.1, fill = "gray90", color = NA) +
  geom_boxplot(aes(fill = Approach),
               position = 'dodge',
               color = 'black',
               size = 0.2,
               outlier.size = 0.4,
               outlier.alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10)) +
  scale_fill_manual(name = NULL,
                    values = c('#84a955','#965da7','#bc5d41'),
                    labels = str_pad_custom) +
  geom_vline(xintercept = 1.5:19.5, linetype = 'dotted') +
  xlab(NULL) + ylab("Relative abundance (%)") +
  ylim(0,32) +
  scale_y_log10(breaks = c(100,10,1,0.1,0.01),
                labels = c(100,10,1,0.1,0.01)) +
  geom_blank(aes(y = 100)) +
  geom_blank(aes(y = 0.00516)) +
  theme(legend.position = 'top',
        axis.text.x = element_text(size = 8))


family = 'sans'
size = 2.3
y_mv4 = 51
y_v4v9 = 81
y_mv9 = 120
height = 5

wilcox_paired_df <- 
  wilcox_paired %>% 
  filter(!(Group %in% c('Discosea','Diplonemea','Kinetoplastida','Prymnesiophyceae') &
           (group1 == 'ampliconV4' | group2 == 'ampliconV4')),
         p.signif != 'ns') %>% 
  mutate(group.num = as.numeric(Group),
         p.signif = case_when(p.signif == '****' ~ '***',
                              TRUE ~ p.signif),
         segment.x = case_when(group1 == 'mTags' ~ group.num - 0.25,
                               group1 == 'ampliconV4' ~ group.num),
         segment.xend = case_when(group1 == 'mTags' & group2 == 'ampliconV4' ~ group.num,
                                  group1 == 'mTags' & group2 == 'ampliconV9' ~ group.num + 0.25,
                                  group1 == 'ampliconV4' ~ group.num + 0.25),
         segment.y = case_when(group1 == 'mTags' & group2 == 'ampliconV4' ~ y_mv4,
                               group1 == 'mTags' & group2 == 'ampliconV9' ~ y_mv9,
                               group1 == 'ampliconV4' ~ y_v4v9),
         pos.signif.x = case_when(group1 == 'mTags' & group2 == 'ampliconV4' ~ (group.num*2-0.25)/2,
                                  group1 == 'mTags' & group2 == 'ampliconV9' ~ group.num,
                                  group1 == 'ampliconV4' ~ (group.num*2+0.25)/2),
         pos.signif.y = case_when(group1 == 'mTags' & group2 == 'ampliconV4' ~ y_mv4 + height,
                                group1 == 'mTags' & group2 == 'ampliconV9' ~ y_mv9 + height,
                                group1 == 'ampliconV4' ~ y_v4v9 + height))


figure2_alt_signif <- 
  figure2_alt +
    geom_segment(data = wilcox_paired_df,
                 aes(x = segment.x,
                     xend = segment.xend,
                     y = segment.y,
                     yend = segment.y),
                 size = 0.2) +
    geom_text(data = wilcox_paired_df,
              aes(x = pos.signif.x,
                   y = pos.signif.y,
                  label = p.signif,
                  family = family),
              size = size)

ggsave(figure2_alt_signif,
       filename = 'figure2.pdf',
       width = 180,
       height = 80,
       units = 'mm')

### NMDS for the three approaches

amp_mtags_table_relab_t <- 
  amp_mtags_table_relab %>% 
  select(Sample2, Group, rel_abun_by_sample) %>% 
  spread(key = Group,
         value = rel_abun_by_sample,
         fill = 0) %>% 
  column_to_rownames('Sample2')

bc_amp_mtags <- 
  vegdist(amp_mtags_table_relab_t, 
          method = "bray")

nmds_amp_mtags <- 
  metaMDS(bc_amp_mtags, 
          autotransform = FALSE, 
          distance = 'bray', 
          trymax = 100)

nmds_amp_mtags_table_md <- 
  as.data.frame(scores(nmds_amp_mtags)) %>% 
  rownames_to_column(var = 'Sample2') %>% 
  as_tibble() %>% 
  separate(Sample2, into = c('Sample','Approach'), sep = '_') %>% 
  left_join(metadata_table, by = c('Sample','Approach')) %>% 
  mutate(Fraction = factor(Fraction, levels = c('Pico','Nano')),
         Layer = factor(Layer, levels = c('Photic','Aphotic')),
         Approach = factor(Approach, levels = c('mTags','ampliconV4','ampliconV9'))) %>% 
  unite(col = 'Station_depth', Station, Depth, sep = '_')

figureS4 <- 
  ggplot(data = nmds_amp_mtags_table_md, 
         aes(x = NMDS1, y = NMDS2)) +
  geom_line(aes(group = Station_depth, 
                color = Layer), 
            size = 0.2,
            alpha = 0.5) +
  geom_point(aes(color = Layer, 
                 shape = Approach), 
             size = 3,
             alpha = 0.6) +
  theme(legend.position = 'right') +
  scale_color_manual(values = c('#fca311', '#1546af'))

ggsave(plot = figureS4,
       filename = 'figureS4.pdf',
       height = 100,
       width = 140,
       units = 'mm')

### PERMANOVA

metadata_for_permanova_amp <- 
  metadata_table %>% 
  unite(Sample, Approach, sep = '_', col = 'Sample2', remove = FALSE) %>% 
  filter(Sample2 %in% rownames(amp_mtags_table_relab_t)) %>% 
  arrange(Sample2)

adonis2(formula = bc_amp_mtags ~ Approach + Layer2 + Ocean,
        data = metadata_for_permanova_amp,
        permutations = 999, 
        method = "bray")
