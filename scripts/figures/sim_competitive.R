# load libraries
library(target)
library(tidyverse)
library(cowplot)

# load data
data("sim_peaks")
data("sim_transcripts")

# select transcripts to change and shift
ap <- associated_peaks(sim_peaks, sim_transcripts, 'tx_id')
nearby <- unique(ap$assigned_region[order(abs(ap$distance))[1:5000]])

ind <- sim_transcripts$tx_id %in% nearby

new_transcripts <- sim_transcripts
new_transcripts$stat1[ind] <- ifelse(new_transcripts$stat1[ind] > 0,
                                     new_transcripts$stat1[ind] * 3,
                                     new_transcripts$stat1[ind])

new_transcripts$stat2[ind] <- ifelse(new_transcripts$stat2[ind] > 0,
                                     new_transcripts$stat2[ind] * -3,
                                     new_transcripts$stat2[ind])

# get targets
dt1 <- direct_targets(sim_peaks, 
                      new_transcripts,
                      'tx_id',
                      'stat1')

dt2 <- direct_targets(sim_peaks, 
                      new_transcripts,
                      'tx_id',
                      'stat2')

dt <- direct_targets(sim_peaks, 
                     new_transcripts,
                     'tx_id',
                     c('stat1', 'stat2'))

# figures
tibble(stat1 = sim_transcripts$stat1,
       stat2 = sim_transcripts$stat2) %>%
  ggplot(aes(x = stat1, y = stat2)) +
  geom_point(alpha = .5, color = 'darkgray') +
  geom_hline(yintercept = 0, lty = 2, color = 'red') +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Factor One', y = 'Factor Two') -> p1

tibble(stat1 = new_transcripts$stat1,
       stat2 = new_transcripts$stat2) %>%
  ggplot(aes(x = stat1, y = stat2)) +
  geom_point(alpha = .5, color = 'darkgray') +
  geom_hline(yintercept = 0, lty = 2, color = 'red') +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Factor One', y = 'Factor Two') -> p2

tibble(distance = ap$distance,
       peak_score = ap$peak_score) %>%
  ggplot(aes(x = distance, y = peak_score)) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Peak Distance', y = 'Peak Score') -> p3

tibble(stat = dt1$stat1,
       rank = dt1$score_rank) %>%
  mutate(group = cut(stat,
                     breaks = quantile(stat, c(0, .25, .75, 1)),
                     labels = c('Down', 'None', 'Up'))) %>%
  na.omit() %>%
  ggplot(aes(x = rank, color = group)) +
  stat_ecdf() +
  theme_bw() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        legend.background = element_rect(fill = NA),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred')) -> p4

tibble(stat = dt2$stat2,
       rank = dt2$score_rank) %>%
  mutate(group = cut(stat,
                     breaks = quantile(stat, c(0, .25, .75, 1)),
                     labels = c('Down', 'None', 'Up'))) %>%
  na.omit() %>%
  ggplot(aes(x = rank, color = group)) +
  stat_ecdf() +
  theme_bw() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        legend.background = element_rect(fill = NA),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred')) -> p5

tibble(stat = dt$stat,
       rank = dt$score_rank) %>%
  mutate(group = cut(stat,
                     breaks = quantile(stat, c(0, .25, .75, 1)),
                     labels = c('Competitive', 'None', 'Cooperative'))) %>%
  na.omit() %>%
  ggplot(aes(x = rank, color = group)) +
  stat_ecdf() +
  theme_bw() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        legend.background = element_rect(fill = NA),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred')) -> p6

plot_grid(p1, p2, p3, p4, p5, p6,
          nrow = 2,
          scale = .9,
          labels = 'AUTO',
          label_size = 10,
          label_fontface = 'plain') %>%
  ggsave(filename = 'manuscript/sim_competitive.png',
         height = 16, width = 24, units = 'cm')
