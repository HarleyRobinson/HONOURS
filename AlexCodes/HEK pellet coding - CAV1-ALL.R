library(DESeq2)
d=read.delim("miRNA_raw_counts_HEK_pellet_C_Pn_P_S_Pr_FILTERED.txt", header=T, row.names=1)
allcav=cbind(d[4:15], d[1:3])
exp_design=data.frame(row.names=colnames(allcav),
                      condition=c("PN", "PN", "PN", "PTRF", "PTRF", "PTRF", "SDPR", "SDPR", "SDPR", "PRKC", "PRKC", "PRKC", "Cont", "Cont", "Cont"),
                      libType=c("X", "X", "X", "X", "X", "X","X", "X", "X", "X", "X", "X", "X", "X", "X"))

condition=relevel(exp_design$condition, "Cont")
dds=DESeqDataSetFromMatrix(countData=allcav, colData=exp_design, design=~condition)
dds=DESeq(dds)

svg('ALL_vs_CAV_plotMA.svg')
plotMA(dds, ylim = c(-3,3))
dev.off()

png('ALL_vs_CAV_plotMA.png')
plotMA(dds, ylim = c(-3,3))
dev.off()

svg('ALL_vs_CAV_plotPCA_dds.svg')
plotPCA(dds)
dev.off()

png('ALL_vs_CAV_plotPCA_dds.png')
plotPCA(dds)
dev.off()

res = results(dds)

res_ordered = res[order(res$padj),]

write.csv(as.data.frame(res_ordered), file="ALL_vs_CAV_results.csv")

res_sig = subset(res_ordered, padj < 0.1)

write.csv(as.data.frame(res_sig), file="ALL_vs_CAV_results_padj01.csv")

res_sig = subset(res_ordered, pvalue < 0.05)

write.csv(as.data.frame(res_sig), file="ALL_vs_CAV_results_pvalue005.csv")

#rld = rlog(dds)
rld = rlogTransformation(dds)

svg('ALL_vs_CAV_results_plotPCA_rld.svg')
plotPCA(rld)
dev.off()

png('ALL_vs_CAV_results_plotPCA_rld.png')
plotPCA(rld)
dev.off()
