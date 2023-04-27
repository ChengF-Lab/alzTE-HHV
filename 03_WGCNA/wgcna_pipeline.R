#!/usr/bin/Rscript

library(WGCNA)
library(reshape2)
library(stringr)
library(dplyr)

options(stringsAsFactors = FALSE)
enableWGCNAThreads() ####change this to allowWGCNAThreads() if you run in Rstudio

exprMat <- "locusTE_gene_FPKM.txt"  #### Each line is a locus-TE expression or a gene expression, columns are all the samples analyzed

type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, quote="", comment="", check.names=F)

dataExpr <- dataExpr %>% mutate_all(~as.numeric(as.character(.)))
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 0),]

dim(dataExprVar)

dataExpr <- as.data.frame(t(dataExprVar))
dataExpr<-log2(dataExpr + 1)
gsg = goodSamplesGenes(dataExpr, verbose = 3)
print(gsg$allOK)

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

sampleTree = hclust(dist(dataExpr), method = "average")


pdf("./WGCNA/hclust.pdf",width = 15,height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
	 ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	 main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
	 ylab="Mean Connectivity", type="n",
	 main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

dev.off()

power = sft$powerEstimate
print(power)


####### One-step network construction and module detection##
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
					   reassignThreshold = 0, mergeCutHeight = 0.25,
					   numericLabels = TRUE, pamRespectsDendro = FALSE,
					   saveTOMs=TRUE, corType = corType,
					   maxPOutliers=maxPOutliers, loadTOMs=TRUE,
					   saveTOMFileBase = paste0(exprMat, ".tom"),
					   verbose = 3)
print(table(net$colors))

# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

# Plot the dendrogram and the module colors underneath
pdf("./WGCNA/wgcna_01.pdf",width = 12,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


###### Module colors and gene correspondence
print(table(moduleLabels))
print(table(moduleColors))


# module eigengene,
MEs = net$MEs

MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)


pdf("./WGCNA/eigengene_heatmap.pdf",width = 8,height = 12)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
					  marHeatmap = c(3,4,2,2), 
					  plotDendrograms = T,
					  xLabelsAngle = 90)
dev.off()

##### TOM plot

TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)


## Loading objects: TOM
TOM <- as.matrix(TOM)
dissTOM = 1-TOM

##### Transform dissTOM with a power to make moderately strong 
##### connections more visible in the heatmap
plotTOM = dissTOM^7

##### Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

##### Call the plot function, it will take such as long time
pdf("./WGCNA/TOMplot.pdf",width = 8,height = 8)
TOMplot(plotTOM,net$dendrograms,moduleColors,main = "Network heatmap plot, all genes")
dev.off()


##### Connect to phenotype
trait <- "./ROSMAP_infected_Metadata.csv" #### Line is sample info and column is traits

##### Read phenotype files
if(trait != "") {
	  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
	                          check.names=FALSE, comment='',quote="")
	  sampleName = rownames(dataExpr)
	  traitData = traitData[match(sampleName, rownames(traitData)), ]
}


##### module correlation with phenotype
if (corType=="pearsoon") {
	  modTraitCor = cor(MEs_col, traitData, use = "p")
	  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
	  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
	  modTraitCor = modTraitCorP$bicor
	  modTraitP   = modTraitCorP$p
}


# signif
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("./WGCNA/module_trait_relationships.pdf",width = 10,height = 8)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
               yLabels = colnames(MEs_col),
			   cex.lab = 0.8,
			   ySymbols = colnames(MEs_col), colorLabels = FALSE,
			   colors = blueWhiteRed(50),
			   textMatrix = textMatrix, setStdMargins = FALSE, 
			   cex.text = 0.8, zlim = c(-1,1),
			   main = paste("Module-trait relationships"))
dev.off()


#### Analyze MEmagenta with specific trait
#### Calculate correlation matrix between module and gene
if (corType=="pearsoon") {
	  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
	  MMPvalue = as.data.frame(corPvalueStudent(
		         as.matrix(geneModuleMembership), nSamples))
} else {
	  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
	  geneModuleMembership = geneModuleMembershipA$bicor
	  MMPvalue   = geneModuleMembershipA$p
}


#### Calculate correlation matrix between trait and gene

if (corType=="pearsoon") {
	geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
	geneTraitP = as.data.frame(corPvalueStudent(
	as.matrix(geneTraitCor), nSamples))
} else {
	geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
	geneTraitCor = as.data.frame(geneTraitCorA$bicor)
	geneTraitP   = as.data.frame(geneTraitCorA$p)
}


###### Combine the two matrix
module = "turquoise"  #### change this to your interest module
pheno = "apoe_genotype"   #### change this to your interest trait
modNames = substring(colnames(MEs_col), 3)

module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))

moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))


pdf("./WGCNA/specific_module_key_gene.pdf",width = 8,height = 8)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
				   xlab = paste("Module Membership in", module, "module"),
				   ylab = paste("Gene significance for", pheno),
				   main = paste("Module membership vs. gene significance\n"),
				   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()


##### Select hub genes
hubs = chooseTopHubInEachModule(dataExpr, colorh=moduleColors, power=power, type=type)

#### extract gene and module information for each trait

traitNames<-names(traitData)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM.", modNames, sep="")
names(MMPvalue) = paste("p.MM.", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr, traitData, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(traitData), sep="")
names(GSPvalue) = paste("p.GS.", names(traitData), sep="")


for (trait in traitNames){
	col_trait<-match(trait,traitNames)
	info0<-data.frame(feature=names(dataExpr),
	modulecolor=net$colors,
	geneTraitSignificance[,col_trait],
	GSPvalue[,col_trait])
	names(info0)<-c("feature","moduleColor",paste("GS.",trait,sep=""),paste("p.GS.",trait,sep=""))
	
	modOrder<-order(abs(modTraitCor[,col_trait]),decreasing=T)
							      
	for (mod in 1:ncol(geneModuleMembership)){
		oldNames = names(info0)
		info0 = data.frame(info0, geneModuleMembership[, modOrder[mod]],MMPvalue[, modOrder[mod]])
		names(info0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),paste("p.MM.", modNames[modOrder[mod]], sep=""))
	}
	    
		n<-match(paste("GS.",trait,sep=""),names(info0))
		featureOrder = order(info0$moduleColor, -abs(info0[,n]))
		info = info0[featureOrder, ]
		
		write.table(info,paste(trait,"_summary_info.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
}


###### results

print(table(net$colors))
print(table(moduleLabels))
print(table(moduleColors))


print("You are all set, congretulations!!!!!")

######## End


