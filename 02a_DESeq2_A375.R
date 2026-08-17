# 20.03.2025, 14.06.2025
setwd("/omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2")

# results from /omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/stringtie_output2
# DE
library(data.table)
library(DESeq2)
library(readr)
library(tximport)

files <- list.files(path="/omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/stringtie_output2", pattern="t_data.ctab", recursive = T, full.names = T) 
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene, txOut=T)
counts <- txi$counts
colnames(counts) <- c("12_01", "13_01",  "16_01",  "18_01",  "2_01",  "20_01", "3_01", "4_01",  "9_01",
                      "12_02", "13_02",  "16_02",  "18_02",  "2_02",  "20_02", "3_02", "4_02",  "9_02",
                      "12_03", "13_03",  "16_03",  "18_03",  "2_03",  "20_03", "3_03", "4_03",  "9_03")

# deseq2
rownames_coldata <- c("12_01", "12_02", "12_03", 
                          "13_01", "13_02", "13_03",
                          "16_01", "16_02", "16_03",
                          "18_01", "18_02", "18_03",
                          "2_01", "2_02", "2_03",
                          "20_01", "20_02", "20_03",
                          "3_01", "3_02", "3_03",
                          "4_01", "4_02", "4_03",
                          "9_01", "9_02", "9_03")
coldata <- data.frame(sample_type=rep("WT", 27), row.names = rownames_coldata)
coldata$sample_type[c(1:6, 25:27)] <- "cKO"
coldata$sample_type[c(7:12, 16:18)] <- "pKO"

counts <- counts[,rownames(coldata)]
write.csv(counts, file="counts_A375.csv")

rowdata = data.frame(tmp[,1:10])
rownames(rowdata) <- tmp$t_name
rowdata = makeGRangesFromDataFrame(rowdata, keep.extra.columns=T)

# repeat annotation for rowdata
repeats<- read.csv("reference_data/hg19_repeats.csv")
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
res1 <- results(dds, contrast=c("sample_type", "cKO", "WT"), format = "GRanges", saveCols=cols_list)
res2 <- results(dds, contrast=c("sample_type", "pKO", "WT"), format = "GRanges", saveCols=cols_list)
#res1$gene_name <- dds@rowRanges$gene_name
#res2$gene_name <- dds@rowRanges$gene_name
res1_s <- res1[which(abs(res1$log2FoldChange)>2 & res1$padj<0.01)]
res2_s <- res1[which(abs(res2$log2FoldChange)>2 & res2$padj<0.01)]
summary(as.factor(substr(names(res1),1,5)))
summary(as.factor(substr(names(res1_s),1,5)))
summary(as.factor(substr(names(res2_s),1,5)))
length(intersect(names(res1_s), names(res2_s)))
int_names <- intersect(names(res1_s), names(res2_s))
summary(as.factor(substr(intersect(names(res1_s), names(res2_s)),1,5)))
write.csv(res1, file="DE_cKO_WT.csv", row.names = names(res1))
write.csv(res2, file="DE_pKO_WT.csv", row.names = names(res2))
write.csv(res1_s, file="DE_cKO_WT_subset.csv", row.names = names(res1_s))
write.csv(res2_s, file="DE_pKO_WT_subset.csv", row.names = names(res2_s))
write.csv(res1_s[int_names], file="DE_cKO_WT_subset_intersect_w_subset_pKO.csv", row.names = int_names)
