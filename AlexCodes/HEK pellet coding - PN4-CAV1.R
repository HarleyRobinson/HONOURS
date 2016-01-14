library(DESeq2)
d=read.delim("miRNA_raw_counts_HEK_pellet_Pn_P_S_C_Pr_FILTERED.txt", header=T, row.names=1)
cav1pn4=cbind(d[10:12], d[1:3])
exp_design=data.frame(row.names=colnames(cav1pn4), 
                      condition=c("Test", "Test", "Test", "Cont", "Cont", "Cont"), 
                      libType=c("X", "X", "X", "X", "X", "X"))
condition=relevel(exp_design$condition, "Cont")
dds=DESeqDataSetFromMatrix(countData=cav1pn4, colData=exp_design, design=~condition)
dds=DESeq(dds)

svg('CAV1_vs_PN4_plotMA.svg')
plotMA(dds, ylim = c(-3,3))
dev.off()

png('CAV1_vs_PN4_plotMA.png')
plotMA(dds, ylim = c(-3,3))
dev.off()

svg('CAV1_vs_PN4_plotPCA_dds.svg')
plotPCA(dds)
dev.off()

png('CAV1_vs_PN4_plotPCA_dds.png')
plotPCA(dds)
dev.off()

res = results(dds)

res_ordered = res[order(res$padj),]

write.csv(as.data.frame(res_ordered), file="CAV1_vs_PN4_results.csv")

res_sig = subset(res_ordered, padj < 0.1)

write.csv(as.data.frame(res_sig), file="CAV1_vs_PN4_results_padj01.csv")

res_sig = subset(res_ordered, pvalue < 0.05)

write.csv(as.data.frame(res_sig), file="CAV1_vs_PN4_results_pvalue005.csv")

#rld = rlog(dds)
rld = rlogTransformation(dds)

svg('CAV1_vs_PN4_results_plotPCA_rld.svg')
plotPCA(rld)
dev.off()

png('CAV1_vs_PN4_results_plotPCA_rld.png')
plotPCA(rld)
dev.off()
