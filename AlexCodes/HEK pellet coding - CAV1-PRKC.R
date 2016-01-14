library(DESeq2)
d=read.delim("miRNA_raw_counts_HEK_pellet_C_Pn_P_S_Pr_FILTERED.txt", header=T, row.names=1)
prkccav1=cbind(d[13:15], d[1:3])
exp_design=data.frame(row.names=colnames(prkccav1), 
                      condition=c("Test", "Test", "Test", "Cont", "Cont", "Cont"), 
                      libType=c("X", "X", "X", "X", "X", "X"))
condition=relevel(exp_design$condition, "Cont")
dds=DESeqDataSetFromMatrix(countData=prkccav1, colData=exp_design, design=~condition)
dds=DESeq(dds)

svg('PRKC_vs_CAV1_plotMA.svg')
plotMA(dds, ylim = c(-3,3))
dev.off()

png('PRKC_vs_CAV1_plotMA.png')
plotMA(dds, ylim = c(-3,3))
dev.off()

svg('PRKC_vs_CAV1_plotPCA_dds.svg')
plotPCA(dds)
dev.off()

png('PRKC_vs_CAV1_plotPCA_dds.png')
plotPCA(dds)
dev.off()

res = results(dds)

res_ordered = res[order(res$padj),]

write.csv(as.data.frame(res_ordered), file="PRKC_vs_CAV1_results.csv")

res_sig = subset(res_ordered, padj < 0.1)

write.csv(as.data.frame(res_sig), file="PRKC_vs_CAV1_results_padj01.csv")

res_sig = subset(res_ordered, pvalue < 0.05)

write.csv(as.data.frame(res_sig), file="PRKC_vs_CAV1_results_pvalue005.csv")

#rld = rlog(dds)
rld = rlogTransformation(dds)

svg('PRKC_vs_CAV1_results_plotPCA_rld.svg')
plotPCA(rld)
dev.off()

png('PRKC_vs_CAV1_results_plotPCA_rld.png')
plotPCA(rld)
dev.off()
