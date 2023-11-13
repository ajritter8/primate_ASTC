library(DEXSeq)
library(dplyr)
library(tidyr)

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/'
jc <- read.csv(paste0(dataDir, 'all_species_junctioncounts_restricted.tsv'), sep="\t")
gff <- readGFF(paste0(dataDir, 'human_splice_lib_events_restricted.gff')) %>%
  filter(type == "exonic_part") %>%
  mutate(bin = paste0(gene_id, ":", exonic_part_number)) %>%
  separate(transcripts, c("A", "B"), "-") %>%
  mutate(event_type = substr(B, 1, 2)) %>%
  dplyr::select(9, 12, 13)

# Generate counts for comparing between conditions (to get significance of splicing events across 2 conditions).
for (sample in unique(jc$sample_name)) {
  outFile <- paste0(dataDir, sample, ".txt")
  x <- jc3 %>% filter(sample_name == sample)
  included <- gff %>% filter(exonic_part_number == "001") %>%
    left_join(., x[c(2, 10)], by=c("gene_id" = "event_id")) %>% rename(counts = avg_ijc)
  excluded <- gff %>% filter(exonic_part_number == "002") %>%
    left_join(., x[c(2, 11)], by=c("gene_id" = "event_id")) %>% rename(counts = avg_ejc)
  y <- rbind(included, excluded) %>% arrange(bin) %>%
    dplyr::select(3, 4) %>%
    mutate_at(2, ~replace_na(.,0)) %>%
    mutate(counts = floor(0.5 + counts))
  write.table(y, outFile, row.names=F, col.names=F, quote=F, sep="\t")
}

runDEXSeq <- function(counts, samples, comparison, subset) {
  testCounts <- counts[subset]
  testSamples <- samples[subset,]
  dxd <- DEXSeqDataSetFromHTSeq(
    testCounts,
    sampleData=testSamples,
    design= ~ sample + exon + condition:exon,
    flattenedfile=flattenedFile)
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
  dxr1 <- DEXSeqResults(dxd)
  gene_qval <- as.data.frame(perGeneQValue(dxr1)) %>%
    tibble::rownames_to_column("event_id") 
  names(gene_qval)[2] <- "event_qval"
  
  dxrDf <- as.data.frame(dxr1) %>%
    separate(transcripts, into = c("A", "B"), sep = ", ") %>%
    mutate(meanA = rowMeans(dplyr::select(., c(colnames(.)[16:17])), na.rm = TRUE),
           meanB = rowMeans(dplyr::select(., c(colnames(.)[18:19])), na.rm = TRUE))
  names(dxrDf)[22] <- paste0(testSamples$condition[1], "_baseMean")
  names(dxrDf)[23] <- paste0(testSamples$condition[4], "_baseMean")
  finalDf <- dxrDf %>%
    filter(grepl('included', A)) %>%
    dplyr::rename(event_id = groupID, exon_id = featureID, exon_padj = padj) %>%
    inner_join(., gene_qval, by="event_id") %>%
    dplyr::select(1, 2, 3, 22, 23, 24, 7, 10)
  
  write.csv(finalDf, paste0(dataDir, '/', comparison, '.csv'), row.names=F)
}

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/chimp'
flattenedFile <- paste0(dataDir, '/human_splice_lib_events_restricted_filtered.gff')
sampleTable <- read.csv(paste0(dataDir, '/chimp_samples.csv'))
countFiles <- list.files(dataDir, pattern=".txt$", full.names=TRUE)

comps <- c('chimp_CvMono', 'chimp_CvPOLH', 'chimp_CvPOLL', 'chimp_CvPOLM')
subsets <- list(c(1, 2, 3, 4), c(1, 2, 5, 6), c(1, 2, 7, 8), c(1, 2, 9, 10))
for (i in length(comps)) {
  runDEXSeq(countFiles, sampleTable, comps[i], unlist(subsets[i]))
}

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/human'
sampleTable <- read.csv(paste0(dataDir, '/human_samples.csv'))
countFiles <- list.files(dataDir, pattern=".txt$", full.names=TRUE)
comps <- c('human_CvMono', 'human_CvPOLH', 'human_CvPOLL', 'human_CvPOLM')
for (i in length(comps)) {
  runDEXSeq(countFiles, sampleTable, comps[i], unlist(subsets[i]))
}

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/orang'
sampleTable <- read.csv(paste0(dataDir, '/orang_samples.csv'))
countFiles <- list.files(dataDir, pattern=".txt$", full.names=TRUE)
comps <- c('orang_CvMono', 'orang_CvPOLH', 'orang_CvPOLL', 'orang_CvPOLM')
for (i in length(comps)) {
  runDEXSeq(countFiles, sampleTable, comps[i], unlist(subsets[i]))
}

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/human_chimp'
sampleTable <- read.csv(paste0(dataDir, '/samples.csv'))
countFiles <- list.files(dataDir, pattern=".txt$", full.names=TRUE)
comps <- c('human_chimp_CvC', 'human_chimp_MvM', 'human_chimp_POLHvPOLH', 
           'human_chimp_POLLvPOLL', 'human_chimp_POLMvPOLM')
subsets <- list(c(1, 2, 11, 12), c(3, 4, 13, 14), c(5, 6, 15, 16), 
                c(7, 8, 17, 18), c(9, 10, 19, 20))
for (i in length(comps)) {
  runDEXSeq(countFiles, sampleTable, comps[i], unlist(subsets[i]))
}

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/human_orang'
sampleTable <- read.csv(paste0(dataDir, '/samples.csv'))
countFiles <- list.files(dataDir, pattern=".txt$", full.names=TRUE)
comps <- c('human_orang_CvC', 'human_orang_MvM', 'human_orang_POLHvPOLH', 
           'human_orang_POLLvPOLL', 'human_orang_POLMvPOLM')
for (i in length(comps)) {
  runDEXSeq(countFiles, sampleTable, comps[i], unlist(subsets[i]))
}

dataDir <- '/public/groups/sanfordlab/people/alexritter/Data/Julia/chimp_orang'
sampleTable <- read.csv(paste0(dataDir, '/samples.csv'))
countFiles <- list.files(dataDir, pattern=".txt$", full.names=TRUE)
comps <- c('chimp_orang_CvC', 'chimp_orang_MvM', 'chimp_orang_POLHvPOLH', 
           'chimp_orang_POLLvPOLL', 'chimp_orang_POLMvPOLM')
for (i in length(comps)) {
  runDEXSeq(countFiles, sampleTable, comps[i], unlist(subsets[i]))
}