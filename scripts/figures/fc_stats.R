library(tidyverse)
library(limma)
library(cowplot)

pharma <- read_rds('~/workingon/curatedAdipoArray/cleandata/pharmacological_perturbation.rds')
ind <- pharma$series_id == 'GSE26207'

eset <- pharma[, ind]

ind <- apply(exprs(eset), 1, function(x) sum(is.na(x)))

eset_d <- eset[ind == 0, 1:12]
eset_g <- eset[ind == 0, 13:24]


res <- list(PPARG = eset_g,
            PPARD = eset_d) %>%
  map(function(x) {
    x$treatment <- factor(x$treatment, levels = c('none', 'drug'))
    mod <- model.matrix(~treatment, data = pData(x))
    fit <- lmFit(x,
                 design = mod)
    fit <- eBayes(fit)
    tt <- topTable(fit,
                   number = Inf,
                   adjust.method = 'fdr',
                   genelist = featureNames(x))
    as.tibble(tt)
  }) %>%
  bind_rows(.id = 'target')

p1 <- res %>%
  ggplot(aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = .5, color = 'darkgray') +
  geom_vline(xintercept = c(-1, 1), lty = 2, color = 'red') +
  geom_hline(yintercept = -log10(.01), lty = 2, color = 'red') +
  facet_wrap(~target) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0,"null"),
        panel.border = element_rect(size = 1.5)) +
  lims(x = c(-2, 2)) +
  labs(x = 'Fold-change (log_2)',
       y = 'P-value (-log_10)')

p2 <- res %>%
  gather(type, value, logFC, t) %>%
  dplyr::select(target, ID, type, value) %>%
  spread(target, value) %>%
  ggplot(aes(x = PPARD, y = PPARG)) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  geom_hline(yintercept = 0, lty = 2, color = 'red') +
  facet_wrap(~type, scales = 'free') +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5))

plot_grid(p1, p2, 
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_size = 10,
          label_fontface = 'plain') %>%
  ggsave(filename = 'manuscript/figures/fc_stats.png',
         height = 7, width = 24, units = 'cm')
