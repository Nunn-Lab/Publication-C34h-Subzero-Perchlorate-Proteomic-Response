setwd("../")# this is the get back to the root folder in case you were one folder in from a previous code
setwd("./Cell Count/")

library(gridExtra)
library(tidyverse)
library(growthcurver)

source("../plainstat.R")


##### This script is to analyze and plot the raw count/OD data from the perchlorate project
data = read_excel_allsheets("Direct_Count_and_OD.xlsx")
cont_col = "black"
perch_col = "lightgreen"

#######edit the files##################

#create the OD data
data[[1]]$Temperature = "-1"
data[[1]]$Group = "Control"
data[[1]][grepl("WC",data[[1]]$`Sample Type`),"Group"] = "Perchlorate"
data[[2]]$Temperature = "-5"
data[[2]]$Group = "Control"
data[[2]][grepl("WC",data[[2]]$`Sample Type`),"Group"] = "Perchlorate"
OD_data = rbind(data[[1]],data[[2]])

#Create the Count data
data[[3]]$Group = "Control"
data[[4]]$Group = "Perchlorate"
colnames(data[[3]])[3] = "Count"
colnames(data[[4]])[3] = "Count"
Count_data = rbind(data[[3]],data[[4]])

#edit the OD data
OD_data$`Incubation day` = gsub("Day","",as.character(OD_data$`Incubation day`)); OD_data$`Incubation day` = as.numeric(OD_data$`Incubation day`)
OD_data$Temperature = as.numeric(OD_data$Temperature); OD_data = OD_data[!grepl("Uninoculated",OD_data$`Sample Type`),] # remove uninoculated rows - because OD same as T0
OD_data = OD_data[!grepl("Uninoculted",OD_data$`Sample Type`),]; colnames(OD_data)[2] = "Day"; OD_data$`Sample Type` = NULL

#edit the count data
Count_data$Day = gsub("Day","",as.character(Count_data$Day)); Count_data$Day = as.numeric(Count_data$Day)
Count_data$Temperature = gsub("Minus_","",as.character(Count_data$Temperature)); Count_data$Temperature = as.numeric(Count_data$Temperature)*(-1)

# create the average data by replicate
Count_data_avg = Count_data %>% group_by(Group,Temperature,Day) %>% summarise ("Count_Avg" = mean(Count), "Count_SD" = sd(Count))
OD_data_avg = OD_data %>% group_by(Group,Temperature,Day) %>% summarise (x = mean(OD_600)); colnames(OD_data_avg)[4] = "OD"
Avg_data = right_join(Count_data_avg,OD_data_avg,by=c('Day' = 'Day','Temperature' = 'Temperature','Group' = "Group"))

#using abundance measurements
Neg1_Cont_Ab = subset.data.frame(Count_data, Group == "Control" & Temperature == -1); #Neg1_Cont$`Incubation day` = Neg1_Cont$`Incubation day` * 24
Neg5_Cont_Ab = subset.data.frame(Count_data, Group == "Control" & Temperature == -5); #Neg5_Cont$`Incubation day` = Neg5_Cont$`Incubation day` * 24
Neg1_WC_Ab = subset.data.frame(Count_data, Group == "Perchlorate" & Temperature == -1); #Neg1_WC$`Incubation day` = Neg1_WC$`Incubation day` * 24
Neg5_WC_Ab = subset.data.frame(Count_data, Group == "Perchlorate" & Temperature == -5); #Neg5_WC$`Incubation day` = Neg5_WC$`Incubation day` * 24

Count_data$logCount = log10(Count_data$Count)

##########get growth curve to then plot##############
Neg1_Cont_fit_Ab = SummarizeGrowth(Neg1_Cont_Ab$Day,Neg1_Cont_Ab$Count); 
Neg5_Cont_fit_Ab = SummarizeGrowth(Neg5_Cont_Ab$Day,Neg5_Cont_Ab$Count); 
Neg1_WC_fit_Ab = SummarizeGrowth(Neg1_WC_Ab$Day,Neg1_WC_Ab$Count); 
Neg5_WC_fit_Ab = SummarizeGrowth(Neg5_WC_Ab$Day,Neg5_WC_Ab$Count);
model_list = list("Neg1_Cont_fit_Ab" = Neg1_Cont_fit_Ab,"Neg5_Cont_fit_Ab" =Neg5_Cont_fit_Ab,"Neg1_WC_fit_Ab" = Neg1_WC_fit_Ab,"Neg5_WC_fit_Ab" = Neg5_WC_fit_Ab)

#set up a table for each model 
log_params = data.frame()
for( i in 1:length(model_list)){
  #i = 1
  model = model_list[[i]]
  log_params[i,"Model"] = names(model_list[i])
  log_params[i,"Group"] = ifelse(grepl("Cont",names(model_list[i])),"Cont","Perchlorate")
  log_params[i,"Temperature"] = ifelse(grepl("1",names(model_list[i])),"-1", "-5")
  log_params[i,"K"] = model$vals$k
  log_params[i,"N0"] = model$vals$n0
  log_params[i,"R"] = model$vals$r
}
Neg1_data = subset(Count_data,Temperature == -1)
Neg1_params = log_params[grepl("-1",log_params$Temperature),]
Neg1_data$col = ifelse(Neg1_data$Group == "Control",cont_col,perch_col)

Neg5_data = subset(Count_data,Temperature == -5)
Neg5_params = log_params[grepl("-5",log_params$Temperature),]
Neg5_data$col = ifelse(Neg5_data$Group == "Control",cont_col,perch_col)

Neg1_data_OD = subset(OD_data,Temperature == -1)
Neg5_data_OD = subset(OD_data,Temperature == -5)

tiff("./Plots/Abundance Curves.tiff",res = 300, width = 9, height = 9, units = "in")
par(mfrow = c(2,2))
#plot -1 Abundance data with Estimated Growth Curve - Growthcurver package
plot(Neg1_data$Day,Neg1_data$Count,col = Neg1_data$col, pch = 16,xlab = "Time (Day)",ylab = "CFU/ml",main = "(A) Cell Abundance at -1C (Estimated Growth Curve)", ylim = c(0,7e+08))
curve((Neg1_params$K[1]) / (1 + (((Neg1_params$K[1]) - (Neg1_params$N0[1])) / (Neg1_params$N0[1])) * exp(-(Neg1_params$R[1]) * x)),xlim = c(0,15),add = T, col = cont_col)
curve((Neg1_params$K[2]) / (1 + (((Neg1_params$K[2]) - (Neg1_params$N0[2])) / (Neg1_params$N0[2])) * exp(-(Neg1_params$R[2]) * x)),xlim = c(0,15),add = T,col = perch_col)
symbols(x = c(7,9,11), y = c(494092083,314850416,605496666),circles = c(1.5,0.8,0.5),add = T,inches = F) # day 7,9,11
symbols(x = c(7,9,11), y = c(113367083,159682083,164457500),circles = c(0.8,0.8,0.5),add = T,inches = F,fg = c(perch_col,perch_col,perch_col)) # day 7,9,11

#plot -5 Abundance data with Estimated Growth Curve - Growthcurver package
plot(Neg5_data$Day,Neg5_data$Count,col = Neg5_data$col, pch = 16,xlab = "Time (Day)",ylab = "CFU/ml",main = "(B) Cell Abundance at -5C (Estimated Growth Curve)", ylim = c(0,7e+08))
curve((Neg5_params$K[1]) / (1 + (((Neg5_params$K[1]) - (Neg5_params$N0[1])) / (Neg5_params$N0[1])) * exp(-(Neg5_params$R[1]) * x)),xlim = c(0,15),add = T, col = cont_col)
curve((Neg5_params$K[2]) / (1 + (((Neg5_params$K[2]) - (Neg5_params$N0[2])) / (Neg5_params$N0[2])) * exp(-(Neg5_params$R[2]) * x)),xlim = c(0,15),add = T,col = perch_col)
symbols(x = c(5,9,11,13), y = c(231025500,387724583,493634167,415068750),circles = c(0.8,0.7,0.9,0.8),add = T,inches = F) # day 7,9,11
symbols(x = c(5,9,11,13), y = c(33558750,116899583,85695833,244265833),circles = c(0.6,0.7,0.6,0.5),add = T,inches = F,fg = c(perch_col,perch_col,perch_col)) # day 7,9,11
legend(-5,-1.5e+08,legend = c("Control","Perchlorate"),col = c(cont_col,perch_col),pch = 16,xpd = NA)

#plot -1 Abundance data with Smoothing 
plot(Neg1_data$Day,Neg1_data$Count,col = Neg1_data$col, pch = 16,xlab = "Time (Day)",ylab = "CFU/ml",main = "(C) Cell Abundance at -1C (Loess Smoothing)", ylim = c(0,7e+08))
smoothLine = loess.smooth(Neg1_data$Day[Neg1_data$Group == "Control"], Neg1_data$Count[Neg1_data$Group == "Control"])
lines(smoothLine$x,smoothLine$y,col = cont_col)
smoothLine = loess.smooth(Neg1_data$Day[Neg1_data$Group == "Perchlorate"], Neg1_data$Count[Neg1_data$Group == "Perchlorate"])
lines(smoothLine$x,smoothLine$y,col = perch_col)
symbols(x = c(7,9,11), y = c(494092083,314850416,605496666),circles = c(1.5,0.8,0.5),add = T,inches = F) # day 7,9,11
symbols(x = c(7,9,11), y = c(113367083,159682083,164457500),circles = c(0.8,0.8,0.5),add = T,inches = F,fg = c(perch_col,perch_col,perch_col)) # day 7,9,11

#plot -5 Abundance data with Smoothing 
plot(Neg5_data$Day,Neg5_data$Count,col = Neg5_data$col, pch = 16,xlab = "Time (Day)",ylab = "CFU/ml",main = "(D) Cell Abundance at -5C (Loess Smoothing)", ylim = c(0,7e+08))
smoothLine = loess.smooth(Neg5_data$Day[Neg5_data$Group == "Control"], Neg5_data$Count[Neg5_data$Group == "Control"])
lines(smoothLine$x,smoothLine$y,col = cont_col)
smoothLine = loess.smooth(Neg5_data$Day[Neg5_data$Group == "Perchlorate"], Neg5_data$Count[Neg5_data$Group == "Perchlorate"])
lines(smoothLine$x,smoothLine$y,col = perch_col)
symbols(x = c(5,9,11,13), y = c(231025500,387724583,493634167,415068750),circles = c(0.8,0.7,0.9,0.8),add = T,inches = F) # day 7,9,11
symbols(x = c(5,9,11,13), y = c(33558750,116899583,85695833,244265833),circles = c(0.6,0.7,0.6,0.5),add = T,inches = F,fg = c(perch_col,perch_col,perch_col)) # day 7,9,11
dev.off()
