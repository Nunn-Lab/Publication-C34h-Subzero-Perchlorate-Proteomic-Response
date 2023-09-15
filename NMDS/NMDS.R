setwd("../")# this is the get back to the root folder in case you were one folder in from a previous code
setwd("./NMDS/")

#load libraries
library(vegan) # this is for WGCNA other functions?
library(dplyr)
library(plotfunctions) # gradientLegend()
library(stringr) #str_split
library(reshape2) #melt
library(ggplot2) #ggplot
library(RColorBrewer) #brewer.pal

#load biostats package source code
source('../biostats.R')
source("../plainstat.R")

#############read in files######################
abacus<-read.csv('ABACUS_output.csv', header=T, row.names=1)
metdata = read.csv("metadata.csv", header = T)

pre_prot = "Q"
samplenum = 45 # this does not include the inital samples
commonname = "X2022_FEB_07_NASA_PCHL"

########### edit metadata file ####################
metdata$Time = gsub("Day","",as.character(metdata$Time))
metdata$Time = gsub("T","",as.character(metdata$Time))
metdata$Temperature = gsub("Minus_","",as.character(metdata$Temperature))
metdata$Temperature = (as.numeric(metdata$Temperature))*(-1)
metdata$Temperature = ifelse(metdata$Group == "Initial",NA,metdata$Temperature)

##############set up NMDS #########################
adjnsaf<-select(abacus, contains('ADJNSAF')) # only use adjusted values for NMDS, 2908 proteins
nsaf.uniq<-cbind(adjnsaf, abacus$ALL_NUMPEPSUNIQ)

#remove contaminants from abacus ( both using the pre-prot and the length of the rowname) & edit the column names
nsaf.uniq<-subset(nsaf.uniq, grepl(paste(pre_prot, collapse="|"), rownames(nsaf.uniq))) #2876 proteins
nsaf.uniq = nsaf.uniq[nchar(rownames(nsaf.uniq)) == 6,]

colnames(nsaf.uniq)<-sub(commonname, "P", colnames(nsaf.uniq))
colnames(nsaf.uniq)<-sub("_ADJNSAF", "", colnames(nsaf.uniq))

# isolate proteins with 2 unique peptides
twopeps<-subset(nsaf.uniq, select= 1:samplenum, nsaf.uniq[,samplenum + 1]>1) # 2755 proteins

#transpose the data
Cp.t<-t(twopeps[,1:samplenum])

#data transformation: Log(x+1)
Cp.tra<-data.trans((Cp.t+1), method='log', plot=F) # this normalizes the data

#do NMDS ordination
Cp.nmds<-metaMDS(Cp.tra, distance='bray', k=2, trymax=100, autotransform=F)

#------------set up colors ( 2 sets, 1 for time and 1 for the temperature tested) -------------
grps= metdata$Group

metdata5 = subset.data.frame(metdata, Temperature == "-5") # select only the metadata @-5C
metdata5 = colorgradient(metdata5,"Time","Blues",9,4,9)

metdata1 = subset.data.frame(metdata, Temperature == "-1") # select only the metadata @-1C
metdata1 = colorgradient(metdata1,"Time","Reds",9,4,9)

metdataInt = subset.data.frame(metdata, Group == "Initial")
metdataInt$rowcols = "grey"

metdata = rbind(metdata5,metdata1,metdataInt)

#plot figure
tiff("NMDS.tiff",res = 300, width = 6, height = 6, units = "in")
fig.Cp<-ordiplot(Cp.nmds, choices=c(1,2), type= "none", display='sites')
plot(fig.Cp$sites, pch =  16 ,col = metdata$rowcols, cex = 2) # pch = 16
ordihull(fig.Cp, groups=grps, draw='lines', label=F, col= 'black') # group by treatment

#legend('topleft', legend= paste(unique(metdata$Temperature), "?C"), pch = unique(metdata$pch))
gradientLegend(valRange = as.numeric(metdata5$Time), color = metdata5$rowcols ,pos = 0.14,side = 3,inside = T,dec = 0,cex = 0.75) # working values: pos = 0.14, side = 3, inside = T, dec = 0
gradientLegend(valRange = as.numeric(metdata1$Time), color = metdata1$rowcols ,pos = 0.45,side = 3,inside = T,dec = 0,cex = 0.75) # working values: pos = 0.45, side = 3, inside = T, dec = 0
text(x = -0.050, y = 0.018, cex = 0.75,labels = "Time (Days) @ -5°C")
text(x = -0.005, y = 0.018, cex = 0.75,labels = "Time (Days) @ -1°C")
text(x = -0.06, y = 0.01, cex = 1,labels = "Control")
text(x = -0.034, y = 0.003, cex = 1, labels = "Initial")
text(x = 0.05, y = 0.01,cex = 1, labels = "WCL")
dev.off()

