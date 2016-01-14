library(DESeq2)
#Always filter and normalise the data first (FilteringandSubsetting.r)
rawdata<- read.csv(file="miRNA_HEK_pellet_filtered.csv", header= TRUE, row.names=1)
Cav2PN4= rawdata[c(1,6, 11, 4, 9, 14)]
Cav1PN4= rawdata[c(2, 7, 12, 4, 9, 14)]
#CaveolinPN4= rawdata[c(3, 8, 13, 1, 5, 14)]
Cav3PN4= rawdata[c(4, 8, 15, 4, 9, 14)]
# Filtering and Normalization
# compute the effective library size by using TMM normalization
#Analysis for Cavin1
exp_designCav1 = data.frame(row.names = colnames(Cav1PN4),
                            condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                            libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_designCav1$condition, "Test")
levels(exp_designCav1)
head(exp_designCav1)
ddsCav1=DESeqDataSetFromMatrix(countData=Cav1PN4, colData=exp_designCav1, design=~condition)
ddsCav1=DESeq(ddsCav1)
resCav1= results(ddsCav1, alpha=0.05)

#for Cavin 2
exp_designCav2 = data.frame(row.names = colnames(Cav2PN4),
                            condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                            libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_designCav2$condition, "Test")
levels(exp_designCav2)
head(exp_designCav2)
ddsCav2=DESeqDataSetFromMatrix(countData=Cav2PN4, colData=exp_designCav2, design=~condition)
ddsCav2=DESeq(ddsCav2)
resCav2= results(ddsCav2, alpha=0.05)

#For Cavin3
exp_designCav3 = data.frame(row.names = colnames(Cav3PN4),
                        condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                        libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_designCav3$condition, "Test")
levels(exp_designCav3)
head(exp_designCav3)
ddsCav3=DESeqDataSetFromMatrix(countData=Cav3PN4, colData=exp_designCav3, design=~condition)
ddsCav3=DESeq(ddsCav3)
resCav3= results(ddsCav3, alpha=0.05)

#comparing the increase between cavins. 
a<- resCav1$log2FoldChange
a1 = which(a>0)
aSig= rownames(rawdata)[which(resCav1$pvalue<0.05)]
aPADJ= rownames(rawdata)[which(resCav1$padj<0.1)]
a<-rownames(rawdata)[a1]
b<- resCav2$log2FoldChange
b1 = which(b>0)
b<-rownames(rawdata)[b1]
bSig= rownames(rawdata)[which(resCav2$pvalue<0.05)]
bPADJ<- rownames(rawdata)[which(resCav2$padj<0.1)]
c<- resCav3$log2FoldChange
c1 = which(c>0)
c<-rownames(rawdata)[c1]
cSig= rownames(rawdata)[which(resCav3$pvalue<0.05)]
cPADJ<- rownames(rawdata)[which(resCav3$padj<0.1)]
Upregshared<- unique(c[c%in%a[a%in%b]])
Sig0.1<- unique(cPADJ[cPADJ%in%aPADJ[aPADJ%in%bPADJ]])
#Sig0.1Upreg<- unique(Upregshared[Upregshared%in%Sig0.1])
#lapply(Sig0.1Upreg, write, "AllCavinPC3PelletPval0.1Up.txt", append=TRUE, ncolumns=1000)

#downregulated
a<- resCav1$log2FoldChange
a1 = which(a<0)
a<-rownames(rawdata)[a1]
b<- resCav2$log2FoldChange
b1 = which(b<0)
b<-rownames(rawdata)[b1]
c<- resCav3$log2FoldChange
c1 = which(c<0)
c<-rownames(rawdata)[c1]
Downregshared<- unique(c[c%in%a[a%in%b]])
Sig0.05Downreg<- unique(Downregshared[Downregshared%in%Sig0.05])
lapply(Sig0.05Downreg, write, "AllCavinPC3PalletPval0.05DOWN.txt", append=TRUE, ncolumns=1000)
