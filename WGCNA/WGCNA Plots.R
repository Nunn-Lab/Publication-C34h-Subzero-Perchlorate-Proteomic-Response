setwd("../")# this is the get back to the root folder in case you were one folder in from a previous code
setwd("./WGCNA/")

# Load WGCNA and flashClust libraries
library(WGCNA)
library(flashClust)
library(gridExtra)
library(cowplot)

#############load files###################################
datExpr = read.csv("datExprPchl.csv")
datTraits = read.csv("datTraitsPchl.csv")
load("TOM.PchlSP6.RData")

softpwr = 6 # tested in the previous code to prep for WGCNA
minModuleSize = 30 # minimum number of peptides clustered per module
MEDissThres = 0.0 #unmerged

############Custom for this dataset########################
datTraits$X = NULL
datTraits = datTraits[,-1] # this removes sample
#datTraits = datTraits[,-4] # this removes the phases
colnames(datTraits)[2] = "WCL Treatment"
colnames(datTraits)[4] = "Growth Phase"
#enableWGCNAThreads()

#############Set Up Files#################################
row.names(datExpr) = datExpr$X
datExpr$X = NULL
datExpr = as.data.frame(t(datExpr)) # transpose

if (dim(datTraits)[1] > dim(datExpr)[1]) {
  datTraits = as.matrix(datTraits)[rownames(datTraits) %in% rownames(datExpr),]
  datTraits = as.data.frame(unclass(datTraits))
}

if (dim(datTraits)[1] < dim(datExpr)[1]) {
  datExpr = datExpr[which(rownames(datExpr) %in% rownames(datTraits)),]
}

dissTOM = 1-TOM
rm(TOM)

##############clustered gene tree########################
geneTree = flashClust(as.dist(dissTOM), method="average")

################Set up Module Size#######################
dynamicMods = cutreeDynamic(dendro=geneTree, distM=dissTOM, pamRespectsDendro=FALSE, minClusterSize=minModuleSize) # deepSplit=2
moduleColors= labels2colors(dynamicMods)
MEs= moduleEigengenes(datExpr, colors= moduleColors,softPower=softpwr)$eigengenes
MEs = orderMEs(MEs)

MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
datTraits = as.data.frame(sapply(datTraits, as.character))

for (i in 1:ncol(datTraits))
  datTraits[,i] = as.numeric(as.factor(datTraits[,i]))

moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Correlate traits --------------------------------------------------------
datTraits = as.data.frame(sapply(datTraits, as.character))

for (i in 1:ncol(datTraits))
  datTraits[,i] = as.numeric(as.factor(datTraits[,i]))
  moduleTraitCor = cor(MEs, datTraits, use= "p")
  textMatrix= paste(signif(moduleTraitCor, 3), "\n(",signif(moduleTraitPvalue, 2), ")", sep = "")

##############create module files for proteins##################
  gene.names = colnames(datExpr)
  
  module_colors= setdiff(unique(moduleColors), "grey")
  # for (color in module_colors){
  #   module=gene.names[which(moduleColors==color)]
  #   write.table(module, paste("./Data/WGCNA/module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
  # }
prot_names = names(datExpr)[moduleColors == "blue"]  
write.table(prot_names, paste("module_blue_prot.txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

prot_names = names(datExpr)[moduleColors == "turquoise"]  
write.table(prot_names, paste("module_turquoise_prot.txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

prot_names = names(datExpr)[moduleColors == "brown"]  
write.table(prot_names, paste("module_brown_prot.txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

############ GS & MM for proteins ##############################
perchlorate = as.data.frame(datTraits$`WCL Treatment`)
modNames = substring(names(MEs),3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="s")
  
geneTraitSignificance = as.data.frame(cor(datExpr, perchlorate, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
names(geneTraitSignificance) = paste("GS.", names(perchlorate), sep="")
names(GSPvalue) = paste("p.GS.", names(perchlorate), sep="")


###############Dendogram#########################################
plot(METree, main= "Clustering of module eigenpeptides", xlab= "", sub= "")

##############Dendogram with Modules###########################
plotDendroAndColors(geneTree,moduleColors,c("Dynamic Tree Cut"), dendroLabels= F,
    hang=0.03, addGuide= TRUE, guideHang=0.05,main=paste("Cluster Dendrogram (Thld=", MEDissThres, ")"))

########### Correlation Heat Map #################
tiff("WGCNA Heatmap.tiff",res = 300, width = 4, height = 3, units = "in")
module_names = names(MEs)
module_names = gsub("^.{1,2}","",module_names)
module_names = paste(module_names,"Module", sep = " ")
module_names = str_to_title(module_names)
module_names = str_wrap(module_names,width = 5)

variable_names = colnames(datTraits)
variable_names = str_wrap(variable_names,width = 5)

labeledHeatmap(Matrix= moduleTraitCor, 
    xLabels= variable_names, xLabelsAngle = 0, xLabelsAdj = 0.5, cex.lab.x = 0.5,
    yLabels= names(MEs),cex.lab.y = 0.5,
    ySymbols= module_names, 
    plotLegend = T, keepLegendSpace = F, legendLabel = "Module Correlation Value (R)",cex.legendLabel = 0.4, cex.lab = 0.4,
    yColorLabels = T,yColorWidth = 0.03,
    #colorLabels= F, 
    colors= blueWhiteRed(50), 
    textMatrix= textMatrix, cex.text= 0.5, 
    setStdMargins = F, 
    zlim= c(-1,1), 
    main= "Module-Trait Correlation",cex.main = 0.6, adj = 0.3
    )
dev.off()
  