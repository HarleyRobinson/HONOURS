library(DESeq2)
rawdata<- read.csv(file="miRNA_raw_counts_PC3_filtered.csv", header= TRUE, row.names=1)
Cav2PN4= rawdata[c(3,7, 11, 1, 5, 9)]
Cav1PN4= rawdata[c(2, 6, 10, 1, 5, 9)]
#CaveolinPN4= rawdata[c(3, 8, 13, 1, 5, 14)]
Cav3PN4= rawdata[c(4, 8, 12, 1, 5, 9)]
GFP= rawdata[c(1, 5, 9)]
rd<- read.csv(file="miRNA_raw_counts_HEK_filtered.csv", header= TRUE, row.names=1)
CaveolinPN4= rd[c(3, 8, 13)]
CaveGFP<- cbind(CaveolinPN4, GFP)
# Filtering and Normalization
# compute the effective library size by using TMM normalization
#Analysis for Cavin1
exp_design = data.frame(row.names = colnames(CaveGFP),
                            condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                            libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_design$condition, "Cont")
levels(exp_design)
head(exp_design)
ddsCav=DESeqDataSetFromMatrix(countData=CaveGFP, colData=exp_design, design=~condition)
ddsCav=DESeq(ddsCav)
resCav= results(ddsCav, alpha=0.05)

#comparing the increase between cavins. 
aPADJ= rownames(CaveGFP)[which(resCav$padj<0.1)]
write.csv(resCav[aPADJ, ], "CaveolinHEKtoPC3GFP.csv")

#Comparing HEK GFP to PC3 cavins. 
GFPPC3<- rawdata[c(1, 5, 9)]
Cav1PC3<- rawdata[c(2, 6, 10)]
Cav2PC3<- rawdata[c(3, 7, 11)]
Cav3PC3<- rawdata[c(4, 8, 12)]
HEKCav1<- rd[c(2, 7, 12)]
HEKCav2<- rd[c(1, 6, 11)]
HekCav3<- rd[c(5, 10, 15)]
HEKCaveolin<- rd[c(3, 8, 13)]
TestGFP<- cbind(HEKCaveolin, GFPPC3)
exp_design = data.frame(row.names = colnames(TestGFP),
                        condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                        libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_design$condition, "Cont")
levels(exp_design)
head(exp_design)
dds=DESeqDataSetFromMatrix(countData=TestGFP, colData=exp_design, design=~condition)
dds=DESeq(dds)
res= results(dds, alpha=0.05)
PADJ<- rownames(TestGFP)[which(res$padj<0.1)]
UpregNo<- which((res$log2FoldChange)>0)
write.csv(res[PADJ, ], "HEKCavevsPC3GFP.csv")

#HEK GFP to GFP PC3
HGFPpCav1<- cbind(Cav1PC3, HekGFP)
exp_design = data.frame(row.names = colnames(HGFPpCav1),
                        condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                        libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_design$condition, "Cont")
levels(exp_design)
head(exp_design)
dds=DESeqDataSetFromMatrix(countData=HGFPpCav1, colData=exp_design, design=~condition)
dds=DESeq(dds)
res= results(dds, alpha=0.05)
PADJ<- rownames(HGFPpCav1)[which(res$padj<0.1)]
