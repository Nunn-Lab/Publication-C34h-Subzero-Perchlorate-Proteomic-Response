setwd("../")# this is the get back to the root folder in case you were one folder in from a previous code
setwd("./WGCNA/")

## libraries ###################################################################
require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(plotly)
require(WGCNA)
require(flashClust)

## functions ###################################################################

# median location normalization 
normalize.median=function (df){
  df_med = apply(as.matrix(df), 2, function(x) median(x, na.rm = T))
  df_med_loc = df_med - median(df_med)
  df.median = t(t(df)-df_med_loc)
  return(df.median)
}

############sample specifics ####################
samplenum = 45
commonname = "X2022_FEB_07_NASA_PCHL"
pre_prot = "Q"

abacus.dat<-read.csv('ABACUS_output.csv', header=T, row.names=1)
meta.dt <- read.csv("metadata.csv") # used to be read_csv

#isolate adjnsaf columns
adjnsaf<-select(abacus.dat, contains('ADJNSAF')) # only use adjusted values for NMDS, 2908 proteins
nsaf.uniq<-cbind(adjnsaf, abacus.dat$ALL_NUMPEPSUNIQ)

#remove contaminants from abacus ( both using the pre-prot and the length of the rowname) & edit the column names
nsaf.uniq<-subset(nsaf.uniq, grepl(paste(pre_prot, collapse="|"), rownames(nsaf.uniq))) #2876 proteins
nsaf.uniq = nsaf.uniq[nchar(rownames(nsaf.uniq)) == 6,]

colnames(nsaf.uniq)<-sub(commonname, "P", colnames(nsaf.uniq))
colnames(nsaf.uniq)<-sub("_ADJNSAF", "", colnames(nsaf.uniq))

# isolate proteins with 2 unique peptides
col.prot<-subset(nsaf.uniq, select= 1:samplenum, nsaf.uniq[,samplenum + 1]>1) # 2755 proteins

# remove the last three samples because they are initial samples
col.prot = col.prot[,-c(43:45)]
meta.dt = meta.dt[-c(43:45),]

## import data #################################################################
df = col.prot
df[df == 0] = NA
df.temp = as.data.frame(df)

## missing/0 data ##############################################################

row.plot <- df.temp %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val))%>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha=0.8) +
  scale_fill_manual(name = "",
                    values = c('steelblue', 'tomato3'),
                    labels = c("Present", "Missing"))+
  labs(x = "Sample",
       y = "Row Number", title = "Missing values") +
  theme(axis.text.x=element_text(angle = -65, hjust = 0, vjust = 0, size = rel(0.8)))

row.plot

## log2 transform ##############################################################

df.log = log2(df)
df.median = normalize.median(df.log)

par(mfrow = c(3,1))
boxplot(df, main = "raw")
boxplot(df.log, main = "log2 scale")
boxplot(df.median, main = "log2 scale + median normalization")

## WGCNA Prep ##################################################################

# expression matrix
datExpr = df.median
datExpr = as.data.frame(t(datExpr))
dat.t = t(df.log)
rownames(dat.t) = meta.dt$Sample
#write.table(dat.t,"./Perchlorate Project/Raw Data/logData.txt")

# check variance=0 and outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #FALSE

if(!gsg$allOK)
{
    if(sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr) [!gsg$goodGenes], collapse=", ")));
    if(sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr) [!gsg$goodSamples], collapse=", ")));
    datExpr=datExpr[gsg$goodSamples, gsg$goodGenes]
}

# trait data
datTraits = meta.dt
rownames(datTraits) = datTraits$Sample

# check name alignment
table(rownames(datTraits)==rownames(datExpr)) 

# fix trait data order
datTraits = datTraits[order(rownames(datTraits)), ]
rownames(datTraits) = datTraits$Sample
table(rownames(datTraits)==rownames(datExpr))

# export data/trait
write.csv(t(datExpr), file = "datExprPchl.csv", row.names = T)
write.csv(datTraits, file = "datTraitsPchl.csv", row.names = T)

# check outlier
A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
cnames = names(datTraits)
rnames = rownames(datTraits)
datTraits = as.data.frame(sapply(datTraits, as.character))
for (i in 1:ncol(datTraits))
  datTraits[,i] = as.numeric(as.factor(datTraits[,i]))
datTraits[is.na(datTraits == T)]=0
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")


# soft power threshold 
powers = c(seq(from =1, to=30, by=1)) 
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, 
                        networkType="signed") #signed
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab= "Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit, signed R^2", type= "n", 
     main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", 
     ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# adjacency into topological overlap matrix
softPower = 6
adjacency = adjacency(datExpr, power = softPower, type = "signed") 
TOM = TOMsimilarity(adjacency, TOMType="signed") 
save(TOM, file="TOM.PchlSP6.RData", compress = F)


