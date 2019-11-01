library(tidyverse)
library(limma)
library(cowplot)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(BCRANK)
library(seqLogo)
library(target)

pharma <- read_rds('~/workingon/curatedAdipoArray/cleandata/pharmacological_perturbation.rds')
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
    as_tibble(tt)
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

peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')

ap <- associated_peaks(peaks, regions, 'SYMBOL')
dt1 <- direct_targets(peaks, regions, 'SYMBOL', 'Ppard')
dt2 <- direct_targets(peaks, regions, 'SYMBOL', 'Pparg')
dt3 <- direct_targets(peaks, regions, 'SYMBOL', c('Pparg', 'Ppard'))

# group targets by directions
direction <- cut(dt3$stat, 
                 breaks = quantile(dt3$stat, c(0, .1, .9, 1)),
                 labels = c('down', 'none', 'up'))

# group peaks by their assigned targets
peak_groups <- split(dt3$SYMBOL, direction)

# reorder peaks and get top n peaks
peak_groups <- lapply(peak_groups, function(x) {
  # get peaks in x targets group
  p <- ap[ap$assigned_region %in% unique(x)]
  
  # order peaks by score
  p <- p[order(p$peak_score, decreasing = TRUE)]
  
  # get n top peaks
  n <- 500
  p[1:n]
})

# make a temporary file
tmp_fasta <- tempfile()

# extract sequences of top peaks from the hg19 genome
pseq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,
               names = peak_groups$up)

# write sequences to fasta file
writeXStringSet(pseq, tmp_fasta)

# set random see
set.seed(123)

# call bcrank with the fasta file
bcout <- bcrank(tmp_fasta, silent = TRUE)

# extract the top motif
top_motif <- toptable(bcout, 1)

# print top motif
top_motif

# plot top motif
png(filename = 'manuscript/figures/motifs_up.png',
    height = 10, width = 10, units = 'cm', res = 300)
plot(top_motif)
dev.off()

# plot sequence logo
png(filename = 'manuscript/figures/logo_up.png',
    height = 10, width = 10, units = 'cm', res = 300)
seqLogo(pwm(top_motif))
dev.off()

# extract sequences of top peaks from the hg19 genome
pseq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,
               names = peak_groups$down)

# write sequences to fasta file
writeXStringSet(pseq, tmp_fasta)

# set random see
set.seed(123)

# call bcrank with the fasta file
bcout <- bcrank(tmp_fasta, silent = TRUE)

# extract the top motif
top_motif <- toptable(bcout, 1)

# print top motif
top_motif

# plot top motif
png(filename = 'manuscript/figures/motifs_down.png',
    height = 10, width = 10, units = 'cm', res = 300)
plot(top_motif)
dev.off()

# plot sequence logo
png(filename = 'manuscript/figures/logo_down.png',
    height = 10, width = 10, units = 'cm', res = 300)
seqLogo(pwm(top_motif))
dev.off()

