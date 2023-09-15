# Standard functions to use

#------import an excel file with multiple sheets into a list of dataframes ----------------
library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#--------create a gradient color pallette based off numeric values ------------------------

# file has to be a numeric vector either as a row,column, or vector
# colscheme = color scheme to input via brewer.pal, character
# colnum = number of colors to pull from the colscheme gradient
colorgradient = function(file,column,colscheme,colnum,mincol,maxcol){
  # adding optional arguments for function
  if (missing(mincol)){
    mincol = 0
  }else {
    mincol = mincol
  }
  
  if (missing(maxcol)){
    maxcol = colnum
  } else {
    maxcol = maxcol
  }
  
  if(missing(column)){
    file1 = file
  } else {
    file1 = as.numeric(file[[column]]) # this will add quotation marks to end end of column name, cannot do file$column
  }

  #set up colors
  colpal1 = brewer.pal(colnum,colscheme)[mincol:maxcol]
  cols = rev(colpal1) 
  pal = colorRampPalette(cols) # Use the following line with RColorBrewer
  
  # #Rank variable for colour assignment
  rowcols = as.numeric(file1) * (-1)# this is the numeric value as a vector
  rc.order = findInterval(rowcols, sort(rowcols))
  rowcols = pal(length(sort(rowcols)))[rc.order]
  ordcol = cbind(file,rowcols)
  return(as.data.frame(ordcol))
}

#-------------convert a KEGG pathway output from a text to a readable dataframe - GHOST KOALA Output ----------------

KEGGPathGHOSTConv = function(filepath){
  #read in libraries
  library(stringr)
  
  KEGG_filt_Des = read.table(filepath, header=F, sep = "\t")
  
  KEGG_filt_Des$Pathway = NA
  
  for ( i in 1:nrow(KEGG_filt_Des)){ # finds the headings( types of pathways) and creates fills the new column with them
    
    heading = ifelse(startsWith(KEGG_filt_Des[i,1],"k") == T, KEGG_filt_Des[i,1], heading) # headings start with k, non-headings start with a space
    KEGG_filt_Des[[i,2]] = ifelse(startsWith(KEGG_filt_Des[i,1],"k") == T, NA, heading)
  }
  KEGG_filt_Des = KEGG_filt_Des[complete.cases(KEGG_filt_Des),] # remove rows that have NAs in them
  
  # split the KEGG description into 4 columns
  for (i in 1:nrow(KEGG_filt_Des)){
    
    # modify the KEGG column
    test_str = KEGG_filt_Des[[i,1]]
    test_str = substring(trimws(test_str),4)#remove the white spaces and the "ko:" at the beginning of the string
    
    segA = strsplit(test_str,split = ";")[[1]][1] # split initial by ; save it
    segB = strsplit(test_str,split = "; ")[[1]][2] # split initial by ; save it
    KEGG_filt_Des$KEGG_ID[[i]] = strsplit(segA, " ")[[1]][1]# split the first segment by space
    KEGG_filt_Des$Prot_Name[[i]] = strsplit(segA, " ")[[1]][2]# split the first segment by space
    KEGG_filt_Des$Prot_Des[[i]] = strsplit(segB,"\\[")[[1]][1]
    KEGG_filt_Des$Enzyme_ID[[i]] = strsplit(segB,"\\[")[[1]][2]
  }
  KEGG_filt_Des$V1 = NULL
  KEGG_filt_Des$KEGG_ID = as.character(KEGG_filt_Des$KEGG_ID)
  KEGGPath = KEGG_filt_Des
  return(KEGGPath)
}


#-------convert BLAST KOALA ouput ( copied into a text file) into a usable dataframe --------------------
# The file path is where the text file is located and Prot_Prefix is the prefix associated to all of your proteins. Can have multiple different prefixes. 
# Prot_prefix format is c("") if only one and c("","", ...) for more than one

# You can save the output of the file through "=" to another variable, or just run the command and it will automatically save into the 
# global environment as "Pathway_DF"
KEGGPathKOALAConv = function(filepath,Prot_Prefix){
  
  #filepath = './Membrane Preperations/July 2023_Coculture/KEGG/KO_Blast_Pathway_raw.txt'
  #Prot_Prefix = c("RX","Q4","R1","A0","tr|","gi|","sp|","PRTC","NM")
  ##read in libraries
  #require(library(stringr))

  #read in clauses
  PROT_Name = Prot_Prefix
  KEGG_Pathway = read.table(filepath, header=F, sep = "\t")
  General_Pathway = c("Metabolism","Genetic Information Processing", "Environmental Information Processing", "Cellular Processes", "Organismal Systems", "Human Diseases")
  
  # isolate the SubClass types
  SubClass = KEGG_Pathway$V1[!KEGG_Pathway$V1 %in% General_Pathway]
  SubClass = SubClass[!startsWith(SubClass,"K") == T] # removes all KO numbers
  SubClass = SubClass[!startsWith(SubClass,"0") == T] # removes all types
  
  #in case there are multiple names to subset by
  for ( i in 1:length(PROT_Name)){
    SubClass = SubClass[!startsWith(SubClass,PROT_Name[[i]]) == T]
  }
  
  #remove any other empty words
  SubClass = SubClass[SubClass != ""]
  
  #isolate the Types
  Types = KEGG_Pathway$V1[startsWith(KEGG_Pathway$V1,"0") == T]
  
  # create new dataframe to fill
  columns = c("Protein", "KO", "Type", "SubClass", "General")
  Tot_df = data.frame(matrix(nrow = 1,ncol = length(columns)))
  colnames(Tot_df) = columns
  
  for ( i in 1:nrow(KEGG_Pathway)){ # finds the headings( types of pathways) and creates fills the new column with them
    # save the value of the general KEGG pathway
    if ((KEGG_Pathway[[i,1]] %in% General_Pathway) == T){ GFdata = KEGG_Pathway[[i,1]]}
    
    # save the value of the SubClass
    if ((KEGG_Pathway[[i,1]] %in% SubClass) == T){SCdata = KEGG_Pathway[[i,1]]}
    
    # fill the column with the type
    if ((KEGG_Pathway[[i,1]] %in% Types) == T){ 
      Tdata = KEGG_Pathway[[i,1]] # all Types have the number of products at the end
      
      amount = as.numeric(sub("\\).*", "", sub(".*\\(", "", Tdata))) #this will calculate the length of the dataframe based on the last number in the line
      df = KEGG_Pathway$V1[(i +1) : (i +(amount*2))] # create a dataframe with all KO and protein under that type
      
      # extract KO and Protein within that dataframe
      for ( j in 1 : length(df)){
        #j = 1
        df_nm = df[[j]]
        
        if (startsWith(df_nm,"K") == T){
          Protein_nm = df[[j + 1]]
          
          if (grepl(",",Protein_nm) == T){
            # create a list of proteins seperated by a comma
            ProtList = scan(text = Protein_nm,sep = ",", what = "")
            
            # creates a new little dataframe by the number of proteins associated with the same KO number and saves it
            Kdata_list = rep(df_nm,length(ProtList))
            Tdata_list = rep(Tdata, length(ProtList))
            SCdata_list = rep(SCdata, length(ProtList))
            Gdata_list = rep(GFdata, length(ProtList))
            
            new_df = cbind.data.frame(ProtList,Kdata_list,Tdata_list,SCdata_list,Gdata_list)
            colnames(new_df) = c("Protein", "KO", "Type","SubClass", "General")
            Tot_df = rbind.data.frame(Tot_df,new_df)
          } else {
            new_df = cbind.data.frame(Protein_nm,df_nm,Tdata,SCdata,GFdata)
            colnames(new_df) = c("Protein", "KO", "Type","SubClass", "General")
            Tot_df = rbind.data.frame(Tot_df,new_df)
          }
        }
      }
    }
  } 
  Tot_df = Tot_df[-1,] # delete the first row
  row.names(Tot_df) = 1:nrow(Tot_df) # reset the row numbers
  
  # clean up the dataframe
  Tot_df$Type = str_replace(Tot_df$Type, " \\s*\\([^\\)]+\\)", "") # remove the number in parentheses at the end of each type
  Tot_df$Type <- gsub("^.{0,5}", "", Tot_df$Type) # remove the numbers before the Type
  Pathway_df = Tot_df
  #return(Pathway_df)
  #assign("Pathway_DF",Pathway_df,envir = .GlobalEnv)
}
  
################arranging 2 ggplots with the same legend###################
library("grid")
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }






