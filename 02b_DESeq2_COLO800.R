# de novo transcript analysis for COLO800 data
# 9.04.2025
setwd("/omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800")

file1 <- "/omics/odcf/project/OE0420/setdb1_human/sequencing/rna_sequencing/view-by-pid/OE0420-SETDB1-h_xxxx/untreated04/paired/merged-alignment/untreated04_OE0420-SETDB1-h_xxxx_merged.mdup.bam"
file2 <- "/omics/odcf/project/OE0420/setdb1_human/sequencing/rna_sequencing/view-by-pid/OE0420-SETDB1-h_xxxx/untreated05/paired/merged-alignment/untreated05_OE0420-SETDB1-h_xxxx_merged.mdup.bam"
file3 <- "/omics/odcf/project/OE0420/setdb1_human/sequencing/rna_sequencing/view-by-pid/OE0420-SETDB1-h_xxxx/untreated06/paired/merged-alignment/untreated06_OE0420-SETDB1-h_xxxx_merged.mdup.bam"
list_s <- c("C15", "C3", "C4", "C7", "C8", "C9")
test <- c(sapply(list_s, function(x) {gsub("xxxx", x, file1)}),
          sapply(list_s, function(x) {gsub("xxxx", x, file2)}),
          sapply(list_s, function(x) {gsub("xxxx", x, file3)}))
test2 <- paste("cp --preserve=links", test, " .")
write.csv(test2, file="bam_copying_command.sh", row.names = F, quote = F)

# 10.04.2025, 16.04.2025
setwd("/omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2")

# results from /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/stringtie_output2
# DE
library(data.table)
library(DESeq2)
library(readr)
library(tximport)

files <- list.files(path="/omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/stringtie_output2", pattern="t_data.ctab", recursive = T, full.names = T) 
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene, txOut=T)
counts <- txi$counts
colnames(counts) <- c("15_04", "3_04",  "4_04",  "7_04",  "8_04",  "9_04",
                      "15_05", "3_05",  "4_05",  "7_05",  "8_05",  "9_05",
                      "15_06", "3_06",  "4_06",  "7_06",  "8_06",  "9_06")

# deseq2
rownames_coldata <- c("15_04", "15_05", "15_06", 
                      "3_04", "3_05", "3_06",
                      "4_04", "4_05", "4_06",
                      "7_04", "7_05", "7_06",
                      "8_04", "8_05", "8_06",
                      "9_04", "9_05", "9_06")
coldata <- data.frame(sample_type=rep("WT", 18), row.names = rownames_coldata)
coldata$sample_type[c(1:3, 7:12)] <- "KO"


counts <- counts[,rownames(coldata)]
write.csv(counts, file="counts_COLO800.csv")

rowdata = data.frame(tmp[,1:10])
rownames(rowdata) <- tmp$t_name
rowdata = makeGRangesFromDataFrame(rowdata, keep.extra.columns=T)

# repeat annotation for rowdata
repeats<- read.csv("../de_novo_rna2/reference_data/hg19_repeats.csv")
repeats$genoName <- gsub("chr", "", repeats$genoName)
colnames(repeats)[6:8] <- c("chr", "start", "end")
repeats = makeGRangesFromDataFrame(repeats, keep.extra.columns=T)

#distance to closest repeaat based of TSS
dist_N <- distanceToNearest(resize(rowdata,1), repeats, ignore.strand=FALSE)
rowdata$nearest_repeat_repName <- NA
rowdata$nearest_repeat_repClass <- NA
rowdata$nearest_repeat_repFamily <- NA
rowdata$dist_nearest_repeat <- NA
rowdata[queryHits(dist_N),]$nearest_repeat_repName<- as.character(repeats[subjectHits(dist_N), ]$repName)
rowdata[queryHits(dist_N),]$nearest_repeat_repClass<- as.character(repeats[subjectHits(dist_N), ]$repClass)
rowdata[queryHits(dist_N),]$nearest_repeat_repFamily<- as.character(repeats[subjectHits(dist_N), ]$repFamily)
rowdata[queryHits(dist_N),]$dist_nearest_repeat<- as.numeric(dist_N@elementMetadata@listData$distance)


dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              rowData = rowdata,
                              colData = coldata,
                              design= ~ sample_type)

dds$sample_type <- relevel(dds$sample_type, ref = "WT")
dds <- DESeq(dds)
cols_list= c("gene_name", "nearest_repeat_repName", "nearest_repeat_repClass", "nearest_repeat_repFamily", "dist_nearest_repeat")
res1_colo800 <- results(dds, contrast=c("sample_type", "KO", "WT"), format = "GRanges", saveCols=cols_list)

res1_colo800_s <- res1_colo800[which(abs(res1_colo800$log2FoldChange)>2 & res1_colo800$padj<0.01)]

summary(as.factor(substr(names(res1_colo800),1,5)))
summary(as.factor(substr(names(res1_colo800_s),1,5)))

write.csv(res1_colo800, file="DE_KO_WT_colo800.csv", row.names = names(res1_colo800))
write.csv(res1_colo800_s, file="DE_KO_WT_subset_colo800.csv", row.names = names(res1_colo800_s))

# intersction with A375 results
int_names1 <- intersect(names(res1_colo800_s), names(res1_s))
int_names2 <- intersect(names(res1_colo800_s), int_names)
summary(as.factor(substr(int_names1,1,5)))
summary(as.factor(substr(int_names2,1,5)))
write.csv(res1_colo800_s[int_names1], file="DE_cKO_WT_subset_intersect_A375_Colo800.csv", row.names = int_names1)
write.csv(res1_colo800_s[int_names2], file="DE_cKO_WT_subset_intersect_A375.pKO.cKO_Colo800.csv", row.names = int_names2)

