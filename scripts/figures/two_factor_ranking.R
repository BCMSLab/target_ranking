library(tidyverse)
library(limma)
library(cowplot)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(target)

pharma <- read_rds('~/workingon/curatedAdipoArray/pharmacological_perturbation.rds')
ind <- pharma$series_id == 'GSE26207'

eset <- pharma[, ind]

ind <- apply(exprs(eset), 1, function(x) sum(is.na(x)))

eset_d <- eset[ind == 0, 1:12]
eset_g <- eset[ind == 0, 13:24]


res <- list(Pparg = eset_g,
            Ppard = eset_d) %>%
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

res2 <- res %>%
  filter(adj.P.Val < .2) %>%
  dplyr::select(SYMBOL = ID, target, t) %>%
  spread(target, t) %>%
  na.omit()

symbol_entrez <- AnnotationDbi::select(org.Mm.eg.db,
                                       res2$SYMBOL,
                                       'ENTREZID',
                                       'SYMBOL')

regions <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                     columns = c(ENTREZID = 'GENEID'),
                     upstream = 100000, downstream = 100000,
                     filter = list(gene_id = symbol_entrez$ENTREZID)) %>%
  as_tibble() %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  left_join(symbol_entrez) %>%
  inner_join(res2) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

peaks <- read_tsv('~/workingon/auto_adipo_diff/autoreg/data/peaks/GSM1199130_peaks.xls',
                  comment = '#') %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

ap <- associated_peaks(peaks, regions, 'SYMBOL')
dt1 <- direct_targets(peaks, regions, 'SYMBOL', 'Ppard')
dt2 <- direct_targets(peaks, regions, 'SYMBOL', 'Pparg')
dt3 <- direct_targets(peaks, regions, 'SYMBOL', c('Pparg', 'Ppard'))

p1 <- tibble(stat = dt1$stat,
       rank = dt1$score_rank) %>%
  mutate(group = cut(stat,
                     breaks = quantile(stat, c(0, .1, .9, 1)),
                     labels = c('Down', 'None', 'Up'))) %>%
  na.omit() %>%
  ggplot(aes(x = rank, color = group)) +
  stat_ecdf() +
  theme_bw() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential (PPARD)', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))

p2 <- tibble(stat = dt2$stat,
             rank = dt2$score_rank) %>%
  mutate(group = cut(stat,
                     breaks = quantile(stat, c(0, .1, .9, 1)),
                     labels = c('Down', 'None', 'Up'))) %>%
  na.omit() %>%
  ggplot(aes(x = rank, color = group)) +
  stat_ecdf() +
  theme_bw() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential (PPARG)', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))

p3 <- tibble(stat = dt3$stat,
             rank = dt3$score_rank) %>%
  mutate(group = cut(stat,
                     breaks = quantile(stat, c(0, .1, .9, 1)),
                     labels = c('Competitive', 'None', 'Cooperative'))) %>%
  na.omit() %>%
  ggplot(aes(x = rank, color = group)) +
  stat_ecdf() +
  theme_bw() +
  geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
  theme(panel.grid = element_blank(),
        legend.position = c(.7,.3),
        panel.border = element_rect(size = 1.5)) +
  labs(x = 'Regulatory Potential', y = 'ECDF', color = '') +
  scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))

plot_grid(p1, p2, p3,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_size = 10,
          label_fontface = 'plain') %>%
  ggsave(filename = 'manuscript/figures/two_factor_ranking.png',
         height = 8, width = 24, units = 'cm')
