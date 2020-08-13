library(tidyverse)
library(reshape2)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# run analysis
source('scripts/analysis.R')

# generate table
tsts <- map2_df(
  list('Down vs Up' = dt1,
       'Up vs Down' = dt2,
       'Cooperate vs Compete' = dt),
  list(c('negative', 'positive'),
       c('positive', 'negative'),
       c('positive', 'negative')),
  function(x, y) {
    r <- x$rank
    g <- cut(x$stat,
             breaks = 3,
             labels = c('negative', 'none', 'positive'))
    broom::tidy(
      test_predictions(r,
                       g,
                       compare = y)
    )
  }, .id = 'Test')

tsts %>%
  mutate(Factor = c('YY1', 'YY2', 'Two Factors'),
         Statistic = round(statistic, 2),
         `P Value` = format(p.value, scientific = TRUE, digits = 1)) %>%
  dplyr::select(Factor, Test, Statistic, `P Value`) %>%
  xtable::xtable(align = 'cllcc') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        math.style.exponents = TRUE,
        file = 'manuscript/tests.tex')
