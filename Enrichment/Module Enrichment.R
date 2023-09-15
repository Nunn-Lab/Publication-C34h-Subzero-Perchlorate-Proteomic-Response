setwd("../")# this is the get back to the root folder in case you were one folder in from a previous code
setwd("./Enrichment/")

library(ggplot2)
library(forcats)
library(gridExtra)
library(clusterProfiler)
library(ggpubr)
library(tidyr)
library(stringr)
library(ggh4x)

########Read in Files####################
module_proteins = list.files("./",pattern = "prot.txt",full.names = T) # WGCNA protein names
module_geneID = list.files("./",pattern = "geneID.tsv",full.names = T) # uniprot ames of the module proteins to be used for KEGG
module_compgo_ratio = read.table("Ratio Comp Summary.txt",header = T, sep = "\t") #CompGO protein ratio for each enriched pathway per module 
end_terms = c("GO:0005737", "GO:0006418","GO:1901607","GO:0009073","GO:0006637","GO:0004812","GO:0005524","GO:0009279","GO:0009063")
module_colors = c("blue","turquoise","brown")


########edit files###################
# Module Protein Names - from WGCNA
file_names = basename(module_proteins); file_names = gsub("_prot.txt","",file_names)
module_proteins = lapply(module_proteins,read.table,header = F); names(module_proteins) = file_names
module_proteins = mapply(cbind,module_proteins, "Module" = gsub("module_","",names(module_proteins)),SIMPLIFY = F)
rm(file_names)

#read in the list of dataframes that has the uniprot names/IDs for each protein ID(original) for each module
file_names = basename(module_geneID);file_names = gsub("_geneID.tsv","",file_names)
module_geneID = lapply(module_geneID,read.table,header = T, sep = "\t"); names(module_geneID) = file_names
module_geneID = module_geneID %>% lapply(. %>% mutate(Gene.Names = strsplit(as.character(Gene.Names)," ")) %>% unnest(Gene.Names))
module_geneID = lapply(module_geneID,function(x) x[(names(x) %in% c("From","Gene.Names"))])
rm(file_names)

# modify CompGO ratio file
module_compgo_ratio$Ratio = str_trim(module_compgo_ratio$Ratio)
module_compgo_ratio$Ratio = sapply(module_compgo_ratio$Ratio,function(x) eval(parse(text = as.character(x))),USE.NAMES = F)

#####################CompGO Enrichment##########################

#subset CompGO enrichment by just the leaves
module_compgo_ratio = module_compgo_ratio[module_compgo_ratio$GO %in% end_terms,]

# edit the column names
module_compgo_ratio = module_compgo_ratio[,c("GO","Description","Type","P.value","Module","Ratio")]
colnames(module_compgo_ratio) = c("ID","Description","Type","pvalue","Module","Protein Ratio")
rownames(module_compgo_ratio) = NULL


############KEGG Enrichment###################################
# enrich the proteins from turquoise module - KEGG
geneList = module_geneID$turquoise_uniprot$Gene.Names
turquoise_KEGG_enrich = as.data.frame(enrichKEGG(gene = geneList,organism = "cps",keyType = "kegg",pvalueCutoff = 0.05))
turquoise_KEGG_enrich$Module = "turquoise"

#enrich the proteins from the blue module - KEGG
geneList = module_geneID$blue_uniprot$Gene.Names
blue_KEGG_enrich = as.data.frame(enrichKEGG(gene = geneList,organism = "cps",keyType = "kegg",pvalueCutoff = 0.05))
blue_KEGG_enrich$Module = "blue"

enrich_modules_KEGG = rbind(turquoise_KEGG_enrich,blue_KEGG_enrich)
enrich_modules_KEGG$Description = sub(" - Colwellia psychrerythraea", "",enrich_modules_KEGG$Description)
rm(geneList,blue_KEGG_enrich,turquoise_KEGG_enrich)

enrich_modules_KEGG$Type = "KEGG Enrichment"
enrich_modules_KEGG = enrich_modules_KEGG[,c("ID","Description","Type","pvalue","Module","GeneRatio")]
enrich_modules_KEGG$GeneRatio = str_trim(enrich_modules_KEGG$GeneRatio)
enrich_modules_KEGG$GeneRatio = sapply(enrich_modules_KEGG$GeneRatio,function(x) eval(parse(text = as.character(x))),USE.NAMES = F)
rownames(enrich_modules_KEGG) = NULL
colnames(enrich_modules_KEGG)[6] = "Protein Ratio"

#############plot Enrichment#################################
df = rbind(module_compgo_ratio,enrich_modules_KEGG) # merge the two dataframes together

df$Type = gsub("biological","GO BP",df$Type)
df$Type = gsub("Cellular","GO CC",df$Type)
df$Type = gsub("Molecular","GO MF",df$Type)

df$Module = str_to_title(df$Module)

tiff("enrich_module.tiff",res = 300, width = 10, height = 12, units = "in")
ggplot(data = df)+
  geom_bar(mapping = aes(x = reorder(Description,`Protein Ratio`),y = `Protein Ratio`,fill = Module,group = Module),stat = "identity",position = position_dodge()) +
  facet_wrap(.~factor(Type,levels = c("GO BP","GO CC","GO MF","KEGG Enrichment")),nrow = 4,scale = "free",strip.position = "right") +
  scale_fill_identity(guide = "legend", labels = c("Blue Module","Turquoise Module")) +
  ylim(0,0.4) +
  coord_flip() + 
  labs(x = "Enrichment Funciton", y = "Protein Ratio") +
  theme_light() +
  force_panelsizes(rows = c(1,0.5,0.5,3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", axis.text.y = element_text(size = 15), strip.text.y = element_text(size = 15))
dev.off()

