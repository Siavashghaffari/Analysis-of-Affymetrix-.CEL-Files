#Dowload all microarray raw data (CEL files) of different conditions and put them in a folder (we don't have to unzip them), 
#Change our working directory to the directory which contains CEL files 
setwd("C:/Siavash/MYH9/Human Data/AffyMatrix/23")
#Load the Affy libary
library(affy)
#read all CEL files in the folder
data <- ReadAffy()
#Normalize the data
eset <- rma(data)
#Save the data to an output file to be used by other programs(Data will be log2 transformed and normalized)
GeneExp = exprs(eset)
#Load annotation Library
library(hgu133plus2.db)
#Only keep the ENTREZID, GENE NAME and GENE Symbol
tab = select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = c("ENTREZID","GENENAME","SYMBOL"))
#Getting mean of multiple affy Id mapping to the same GENE
exprSet = t(sapply(split(tab[,1], tab[,4]), function(ids){
  colMeans(GeneExp[ids,,drop=FALSE])
}))
# Check out all our conditions
C_name <- colnames(exprSet)
# Calculate the mean of 40 cases of healthy and 40 cases of Atherosclerosis 
Healthy_mean = apply(exprSet[, C_name[grep("MAMMARY", C_name)]], 1, mean)
Athero_mean = apply(exprSet[, C_name[grep("STANS", C_name)]], 1, mean)

#Reorder data and assmlble to create organized dataset
Healthy_dataset = exprSet[, C_name[grep("MAMMARY", C_name)]]
Athero_dataset = exprSet[, C_name[grep("STANS", C_name)]]
dataset = cbind(Healthy_dataset, Athero_dataset)

#Calculate p_values of statistical significance test using paried T-test to check if healthy and athero difference is significant
p_value_all_genes = apply(dataset, 1, function(x) {
  t.test(x[1:40], x[41:80],"two.sided", paired = TRUE)$p.value
})

#Generate FDR corrected p-values 
FDR =  p.adjust(p_value_all_genes, method = "BH")

#Find athero-prone genes as genes which has greater expression in athero patients with FDR lower than 0.005
A_prone = names(FDR[(Athero_mean>Healthy_mean)&(FDR<0.005)])
#The differential expression final list of athero-prone genes
DE_probeset = dataset[A_prone,]

#change the gene names to theri symbols (if needed)
tmp = rownames(DE_probeset)
tab_Match=cbind(unique(tab$ENTREZID), unique(tab$SYMBOL))
eid=which(tab_Match[, 1] %in% tmp)
tmp = tab_Match[eid,2]

#three ways to create heta maps
heatmap(DE_probeset)
hr <- hclust(as.dist(1-cor(t(DE_probeset), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(DE_probeset, method="spearman")), method="complete")
heatmap(DE_probeset, Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc), scale="col",ColSideColors=heat.colors(length(hc$labels)), scale="col")

library(pheatmap)
pheatmap(DE_probeset)

library(gplots)
heatmap.2(DE_probeset)

#Clustering
disty = dist(DE_probeset)
hc = hclust(disty)
plot(hc)

#Save the final list of athero_prones for furthur analysis in other programs
write.csv(DE_probeset,file="Heatmap.csv")
