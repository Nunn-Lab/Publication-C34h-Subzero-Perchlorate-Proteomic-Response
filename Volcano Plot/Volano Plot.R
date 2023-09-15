setwd("../")# this is the get back to the root folder in case you were one folder in from a previous code
setwd("./Volcano Plot/")

#This code is for creating the volcano plot of differential abundance between Perchlorate and Control with Module Colors
library(ggrepel)
library(dplyr)

##################load in files#####################
Qspec_CW = read.csv("QSPEC_CW.txt_qspec_fdr",sep = "\t", header = T,row.names = 1)
abacus = read.csv('ABACUS_output.csv', header=T, row.names = 1)
module_proteins = list.files("./",pattern = "prot.txt",full.names = T)
KEGG_34H = read.csv('34h_annotations.csv', header=T, row.names=1) # original KEGGs for entire C. psychrerythraea 34H

logfold_limit = 0.5
zstat_limit = 2

module_colors = c("blue","turquoise","brown")

################# modify files ############################
# Module Protein Names
file_names = basename(module_proteins); file_names = gsub("_prot.txt","",file_names)
module_proteins = lapply(module_proteins,read.table,header = F); names(module_proteins) = file_names
module_proteins = mapply(cbind,module_proteins, "Module" = gsub("module_","",names(module_proteins)),SIMPLIFY = F)
rm(file_names)

#QSpec
Qspec_CW$Protein = rownames(Qspec_CW)
Qspec_CW["WGCNA group"] = "Other Protein"
Qspec_CW[grepl(paste(module_proteins[["module_turquoise"]]$V1,collapse = "|"),rownames(Qspec_CW)),"WGCNA group"] = "turquoise"
Qspec_CW[grepl(paste(module_proteins[["module_blue"]]$V1,collapse = "|"),rownames(Qspec_CW)),"WGCNA group"] = "blue"
Qspec_CW[grepl(paste(module_proteins[["module_brown"]]$V1,collapse = "|"),rownames(Qspec_CW)),"WGCNA group"] = "brown"

Qspec_CW$p.value = 2*pnorm(-abs(Qspec_CW$Zstatistic))
Qspec_CW$log.p.value = -log10(Qspec_CW$p.value)
Qspec_CW$ABACUS_DES = abacus$DEFLINE[match(Qspec_CW$Protein,rownames(abacus))]
Qspec_CW["PSM_Sum"] = as.numeric(apply(Qspec_CW[,2:43],1,sum))
Qspec_CW["PSM_Avg"] = as.numeric(apply(Qspec_CW[,2:43],1,mean))
Qspec_CW[,2:43] = NULL
Qspec_CW$ABACUS_DES = gsub("OS=.*","",Qspec_CW$ABACUS_DES)

Qspec_CW = merge.data.frame(Qspec_CW,KEGG_34H,by.x = "Protein",by.y = "Protein", all.x = T)

# adding this clause to modify the log-p-value of those whose z-statistic cause the p-value to be 0 which make the log inifinite. I am just making the log-p-value of those that
# have infinite value to become the maximum value found in the dataset. 37 is the threshold
for( i in 1:nrow(Qspec_CW)){
  if(abs(Qspec_CW$Zstatistic[i]) > 37){
    Qspec_CW$log.p.value[i] = max(Qspec_CW$log.p.value[is.finite(Qspec_CW$log.p.value)])
  }
}

QSPEC_CW_sig = Qspec_CW[abs(Qspec_CW$LogFoldChange) > logfold_limit & abs(Qspec_CW$Zstatistic) > zstat_limit, ]

#########################Statistics##################################
Inc_sig = filter(QSPEC_CW_sig,QSPEC_CW_sig$LogFoldChange > 0)
Dec_Sig = filter(QSPEC_CW_sig,QSPEC_CW_sig$LogFoldChange < 0)

#calculate the percentage of blue module proteins in the significantly increased proteins
Inc_sig %>% group_by(`WGCNA group`) %>% summarise(total_count = n(),.groups = 'drop') 


################QSPEC Analysis#######################################
df3 = Qspec_CW 

# add a column group and fill it with not-significant 
df3["group"] <- "NotSignificant"

#create groups to plot points based off the log fold change the log of the p-value. The log of p-value is better for volcano plots as it creates a better spread
df3[which(abs(df3['LogFoldChange']) >= 0.5 & between (df3$log.p.value, 1.301, 2) & df3["WGCNA group"] == "turquoise"),"group"] = "Turquoise_0.05" # between 0.05 and 0.01
df3[which(abs(df3['LogFoldChange']) >= 0.5 & df3$log.p.value > 2 & df3["WGCNA group"] == "turquoise"),"group"] = "Turquoise_0.01" # between 0.01 and 

df3[which(abs(df3['LogFoldChange']) >= 0.5 & between (df3$log.p.value, 1.301, 2) & df3["WGCNA group"] == "blue"),"group"] = "Blue_0.05" # between 0.05 and 0.01
df3[which(abs(df3['LogFoldChange']) >= 0.5 & df3$log.p.value > 2 & df3["WGCNA group"] == "blue"),"group"] = "Blue_0.01" # between 0.01 and 

df3[which(abs(df3['LogFoldChange']) >= 0.5 & between (df3$log.p.value, 1.301, 2) & df3["WGCNA group"] == "Other Protein"),"group"] = "Other_0.05" # between 0.05 and 0.01
df3[which(abs(df3['LogFoldChange']) >= 0.5 & df3$log.p.value > 2 & df3["WGCNA group"] == "Other Protein"),"group"] = "Other_0.01" # between 0.01 and 

#add the colors for the outlines/ fill circles 
cols<- c(NotSignificant="#FFFFFF", Other_0.05 = "#a9a9a9",Other_0.01 = "#656565",Turquoise_0.05="#9fefe7", Turquoise_0.01="#42E0D0", Blue_0.05 = "#7CB9E8", Blue_0.01 = "#3457D5")
outline<- c(NotSignificant="gray", Other_0.05 = "#555555", Other_0.01 = "#555555" ,Turquoise_0.05="#008081", Turquoise_0.01="#009092", Blue_0.05 = "#003d80", Blue_0.01 = "#003d80")

#specify the proteins to label ( top 5 log fold change, and bottom 5 log fold change)
high = df3[order(df3$LogFoldChange,decreasing = T),][1:5,]
low = df3[order(df3$LogFoldChange,decreasing = F),][1:5,]
sig_prot = df3[df3$log.p.value > 115,]
A5N.only = rbind.data.frame(high,low,sig_prot)
bad_name = c("protein", "family")
A5N.only$ABACUS_DES = gsub(paste(bad_name, collapse = "|"), "", A5N.only$ABACUS_DES)
A5N.only$ABACUS_DES = str_wrap(A5N.only$ABACUS_DES,width = 20)
A5N.only$ABACUS_DES = gsub(",.*","",A5N.only$ABACUS_DES)

tiff("QSPEC_label_module3.tiff",res = 300, width = 15, height = 7, units = "in") # width = 15, height = 7
# create the volcano plot
ggplot(df3,aes(x = LogFoldChange, y = log.p.value, cex = 2, fill = factor(group), color = factor(group))) + 
  geom_point(shape = 21) +
  ylim(0,250) +
  scale_colour_manual(values = outline) +
  scale_fill_manual(values = cols) +
  scale_shape_identity()+
  xlab("Log Fold Change") +
  ylab(" -Log10 (P-value)")+
  theme_classic(base_size = 8) +
  geom_hline(yintercept = 0, linetype= "dashed", colour="gray")+
  geom_vline (xintercept = 0.5, linetype='dashed', colour="gray")+
  geom_vline (xintercept = -0.5, linetype='dashed', colour="gray") +
  geom_text_repel(data=A5N.only, #separate dataset for specific labels
                  aes(label= ABACUS_DES, size = 2), #call protein name, font size?
                  force        = 2.5, #how far point is from bubble
                  nudge_x      = 0.2, #push along x axis
                  nudge_y = -0.5,
                  direction    = "both", #bubbles move x and y directions
                  hjust        = 1, # orientation to the right
                  show.legend = F,
                  box.padding = 0.9 #repel boxes from each other
                  
  ) + 
  theme(legend.position = "bottom", text = element_text(size=15),legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size = 5,shape = 21))) # before shape = 21
dev.off()

