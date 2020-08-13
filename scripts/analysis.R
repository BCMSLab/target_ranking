library(tidyverse)
library(reshape2)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(target)

# load peaks
peaks1 <- GRanges(import.bed('data/Oth.Utr.05.YY1.AllCell.bed'))
mcols(peaks1) <- NULL
peaks2 <- GRanges(import.bed('data/Oth.Utr.05.YY2.AllCell.bed'))

# merge and name peaks
peaks <- reduce(subsetByOverlaps(peaks1, peaks2))
peaks$name <- paste0('common_peak_', 1:length(peaks))

# load expression data
diff_exp <- map_df(c('data/DataSet_01_18.tsv',
                     'data/DataSet_01_19.tsv'), 
                   function(x) {
                     read_tsv(x, col_names = FALSE) %>%
                       dplyr::select(2,3, 7, 9) %>% #9
                       setNames(c('tf', 'gene', 'fc', 'pvalue')) %>%
                       filter(tf %in% c('YY1', 'YY2'))
                   }) %>%
  na.omit() 

diff_exp2 <- diff_exp %>%
  unite(new, fc, pvalue) %>%
  spread(tf, new) %>%
  separate(YY1, c('YY1_fc', 'YY1_pvalue'), sep = '_') %>%
  separate(YY2, c('YY2_fc', 'YY2_pvalue'), sep = '_') %>%
  mutate_at(vars(starts_with('YY')), as.numeric)
#%>% filter_at(vars(ends_with('pvalue')), function(x) x < .2)

# load genome data
symbol_entrez <- select(org.Hs.eg.db,
                        unique(diff_exp$gene),
                        'ENTREZID',
                        'SYMBOL') %>%
  setNames(c('gene', 'gene_id'))

genome <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                      filter = list(gene_id = symbol_entrez$gene_id),
                      columns = c('tx_id', 'tx_name', 'gene_id')) %>%
  promoters(upstream = 100000) %>%
  as_tibble() %>%
  filter(length(gene_id) > 1) %>%
  mutate(gene_id = as.character(gene_id))

regions <- genome %>%
  inner_join(symbol_entrez) %>%
  inner_join(diff_exp2) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# run target analysis
ap <- associated_peaks(peaks, regions, 'tx_id')

dt1 <- direct_targets(peaks1, regions, 'tx_id', 'YY1_fc')
dt2 <- direct_targets(peaks2, regions, 'tx_id', 'YY2_fc')
dt <- direct_targets(peaks, regions, 'tx_id', c('YY1_fc', 'YY2_fc'))

# figure
regulated_groups <- map(c(YY1 = dt1, YY2 = dt2), function(x) {
  df <- as_tibble(mcols(x)) %>%
    mutate(group = cut(stat, 
                       breaks = 3,
                       labels = c('Down', 'None', 'Up')))
  split(df$gene, df$group)
}) %>%
  melt() %>%
  as_tibble() %>%
  unique() %>%
  setNames(c('gene', 'dir', 'tf'))
