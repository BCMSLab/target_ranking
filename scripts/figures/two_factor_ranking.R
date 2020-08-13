library(tidyverse)
library(reshape2)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(cowplot)

# run analysis
source('scripts/analysis.R')

# figure
p1 <- mcols(dt1) %>%
  as_tibble() %>%
  mutate(group = cut(stat, 
                     breaks = 3,
                     labels = c('Down', 'None', 'Up'))) %>%
  ggplot(aes(x = score_rank, color = group)) +
  stat_ecdf() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential (YY1)', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))

p2 <- mcols(dt2) %>%
  as_tibble() %>%
  mutate(group = cut(stat, 
                     breaks = 3,
                     labels = c('Down', 'None', 'Up'))) %>%
  ggplot(aes(x = score_rank, color = group)) +
  stat_ecdf() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential (YY2)', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))

p3 <- mcols(dt) %>%
  as_tibble() %>%
  mutate(group = cut(stat, 
                     breaks = 3,
                     labels = c('Competitive', 'None', 'Cooperative'))) %>%
  ggplot(aes(x = score_rank, color = group)) +
  stat_ecdf() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Interaction', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))

plot_grid(p1, p2, p3,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10) %>%
  ggsave(filename = 'manuscript/two_factor_ranking.png',
         width = 24, height = 8, units = 'cm')
