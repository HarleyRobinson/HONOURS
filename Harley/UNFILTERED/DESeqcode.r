#library(DESeq2)
rawdata <- read.csv("miRNA_raw_counts_HEK_Pellet.csv", header= TRUE, row.names=1)
rawdataPC3<- read.csv("miRNA_raw_counts_Correct_Name_ReOrdered.csv", header= TRUE, row.names=1)
rawdataPC3<- rawdataPC3[c(31:42)]
Cav2PN4= rawdataPC3[c(3, 7, 11 ,1, 5, 10)]
Cav1PN4= rawdataPC3[c(2, 6, 10, 1, 5, 10)]
#CaveolinPN4= rawdata[c(3, 8, 13, 1, 5, 14)]
Cav3PN4= rawdataPC3[c(4, 8, 12, 1, 5, 10)]
# Filtering and Normalization
# compute the effective library size by using TMM normalization

exp_design = data.frame(row.names = colnames(Cav3PN4),
                        condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                        libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_design$condition, "Test")
levels(exp_design)
head(exp_design)

dds=DESeqDataSetFromMatrix(countData=Cav3PN4, colData=exp_design, design=~condition)
dds=DESeq(dds)

res=results(dds, alpha=0.05)
OrderedRes=res[order(res$padj),]
PvalOrdered=res[order(res$pvalue, na.last= NA), ]
#na.last=true means that all NAs will be sortedto end
write.csv(as.data.frame(OrderedRes), file="PN4_Vs_Cav3_all_results.csv")
Sig0.1res=subset(OrderedRes, padj<0.1)
#^^ finds only values of 0.1 or higher based on the Benjamini-hochberg adjustment p.value. 
write.csv(as.data.frame(Sig0.1res), file="PN4_Cav3_0.1padj.csv")
Sig0.5pval=subset(OrderedRes, pvalue <0.05)
write.csv(as.data.frame(Sig0.5pval), file="PN4_Cav3_Pval0.05.csv")
rlogoddsdds=rlogTransformation(dds)

png("plotPCA_rlogoddsCav3.png")
plotPCA(rlogoddsdds)
dev.off()

#graph below gives best mirna given high BHpval annnnd pvalue. 
png("plotMAPN4_Cav3_0.1minpadj.png")
DESeq2::plotMA(res, ylim=c(-3, 3))
topGene<- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

png("plotMAPN4_Cav3_0.05Pvalue.png")
DESeq2::plotMA(res, ylim=c(-3, 3))
topGene<- rownames(head(PvalOrdered))
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()
