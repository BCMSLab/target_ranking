library(tidyverse)
library(reshape2)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(cowplot)
library(ggupset)
library(target)

# run analysis
source('scripts/analysis.R')

# generate figure
c1 <- diff_exp %>%
  ggplot(aes(x = fc, y = -log10(pvalue))) +
  geom_point(color = 'darkgray', alpha = .5) +
  facet_wrap(~tf, scales = 'free_x') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Fold-change (log_2)',
       y = 'P-value (-log_10)')

c2 <- diff_exp2 %>%
  ggplot(aes(x = YY1_fc, y = YY2_fc)) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_hline(yintercept = 0, lty = 2, color = 'red') +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  theme_bw() +
  lims(x = c(-3, 3), y = c(-3, 3)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'YY1 Fold-change (log_2)',
       y = 'YY2 Fold-change (log_2)')

c3 <- regulated_groups %>%
  filter(dir != 'None') %>%
  unite(group, dir, tf) %>%
  group_by(gene) %>%
  summarise(groups = list(group)) %>%
  ggplot(aes(x = groups)) +
  geom_bar() +
  scale_x_upset() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA, size = 1.5)) +
  labs(x = 'Target Sets', y = "Count of Targets")

plot_grid(c1, c2, c3,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10) %>%
  ggsave(filename = 'manuscript/fc_stats.png',
         width = 24, height = 8, units = 'cm')
