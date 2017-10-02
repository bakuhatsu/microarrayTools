##############################
## Sven Nelson              ##
## 3/12/2015                ##
## Function: countTFs       ##
##############################

countTFs <- function(geneListA, geneListB, A, B, labA, labB, title, latexTable, noPlot) {
	# TFplot(TFcounts(DvsWTdryset$UP), TFcounts(DvsWTdryset$DN))
	TFplot(TFcounts(geneListA), TFcounts(geneListB),A, B, labA, labB, title, latexTable, noPlot)
}

TFcounts <- function(geneList) {
	# geneList is a list of AT numbers (not case sensitive)
  geneList <- toupper(geneList) # removes case-sensitivity
  ## Example code: UP and DOWN
  # ARvsD12h_UP_TFs <- TFcounts(ARvsD12hset$UP)
  # ARvsD12h_DN_TFs <- TFcounts(ARvsD12hset$DN)

  # This version takes a geneList of differentially expressed genes and returns 
  # a dataframe with number of hits per TF type for easy plotting
  # Future versions may include the ability to return the list of genes in a given type
  
  # Please ensure that TFsWithType.tsv is in your working directory or TFsWithType exists
  if(is.null(TFsWithType)) {
    TFsWithType <- read.table("TFsWithType.tsv", header=TRUE, sep="\t")
  }
  if(is.null(TFsWithType)) {
    writeLines("\nUnable to locate 'TFsWithType' file.")
  }
  
  # Subset TFsWithType to a dataframe that only included genes in geneList with TF family
  TFhits <- subset(TFsWithType, Protein_ID %in% geneList)
  
  # Lists the different TF families present in geneList
  familiesRepresented <- unique(TFhits$Family)
  

	# Start by creating a data frame with 1 column for each TF family category
  TFfamilyCounts <- data.frame(matrix(NA, nrow = 1, ncol = length(familiesRepresented)))
  colnames(TFfamilyCounts) <- familiesRepresented
  
  # Fill in the dataframe with the number of hits in each TF family
  for(i in 1:length(TFfamilyCounts[1,])) { 
    #TFfamilyCounts[,i] <- countHits(TFhits, familiesRepresented[i])
    
  	TFfamilyCounts[,i] <- nrow(subset(TFhits, Family %in% familiesRepresented[i])) 
  
  	# Here is would be easy to modify this code to return a table of hits
  	# subset(TFhits, familiesRepresented[i] %in% Family) # make use of this code
  }
    
  # returns a data.frame of counts for each TF family
  return(TFfamilyCounts)
}

### Takes two TFcounts objects and plots them using ggplot2
TFplot <- function(TFcountsA, TFcountsB, A = "UPreg", B = "DOWNreg", labA = A, labB = B, title = "TF families", latexTable=F, noPlot=F) {
  ## Example code: UP and DOWN
  # TFplot(ARvsD12h_UP_TFs, ARvsD12h_DN_TFs)
  ## Example code: diffA and diffB
  # TFplot(ARvsD12h_TFs, WTvsD12h_TFs,A="ARvsD12h",B="WTvsD12h")
  
  #### Preparing data ####
  # Use TFcounts to define TFcounts:
  # TFcountsA => UP, TFcountsB => DOWN
  
  # Make combined list of colnames (unique), ordered alphabetically
  
  combinedNames <- sort(unique(c(colnames(TFcountsA),colnames(TFcountsB)))) # combined list
  # Create a combined dataframe
  TFcluster.df <- data.frame(matrix(NA, nrow = length(combinedNames), ncol = 3))
  rownames(TFcluster.df) <- combinedNames
  colnames(TFcluster.df) <- c("Family","A","B")
  TFcluster.df$Family <- combinedNames
  
  # Fill in the dataframe with the number of hits in each TF family
  for(i in 1:length(combinedNames)) { 
    #print(combinedNames)
    if (combinedNames[i] %in% colnames(TFcountsA)) {
      #print("enteredA")
      TFcluster.df$A[i] <- TFcountsA[1,combinedNames[i]]
    }  
    if (combinedNames[i] %in% colnames(TFcountsB)) {
      #print("enteredB")
      TFcluster.df$B[i] <- TFcountsB[1,combinedNames[i]] 
    }
  }
  #print(TFcluster.df)
  # Add code here to output a latex table (replace NAs with 0s)
  if(latexTable==TRUE) { # p-values get cut off, need to be rounded...
    TFclust <- TFcluster.df
    TFclust$A[is.na(TFclust$A)] <- 0
    TFclust$B[is.na(TFclust$B)] <- 0
    TFtable <- TFclust[,2:3]
    
    require(xtable)
    writeLines("\\documentclass[border={(0.5pt) (0.8pt) (1pt) (1pt)}]{standalone}

\\begin{document}
\\SweaveOpts{concordance=TRUE}\n")
    print.xtable(xtable(TFtable,digits = c(0,0,0),floating=FALSE))
    writeLines("\\end{document}")
    #print(tab,type="html")
  }
  
  #if (returnTable) { 
  #  return(TFcluster.df)
  #} # nothing after this point will be run if table was returned
  if(!noPlot) {
    require(reshape) # for melt
    
    TFcluster.long <- melt(TFcluster.df,
                           ## ID variables: 
                           # variables to keep but not split apart on
                           id.vars="Family",
                           # Measure variables: the source columns
                           measure.vars=c("B","A"),
                           # Name of the destination column that 
                           # will identify the original
                           # column that the measurement came from
                           variable_name="Comparison"
    )
    
    # Reorder the data by Ontology
    TFcluster.long$Family <- factor(TFcluster.long$Family, levels = combinedNames)
    
    #require(scales) # For percent_format() in plot
    require(ggplot2)
    
    if(A == "UPreg" & B == "DOWNreg") { # UP => pinkish and DOWN => blueish
      colorA <- "#E82A76" # UP => pinkish
      colorB <- "#3B2DD6" # DOWN => blueish
    } else {
      colorA <- "#E69F00" # A => yellow
      colorB <- "#999999" # B => gray
    }
    
    #### Now for the TAGGIT plot ####
    ggplot(data=TFcluster.long, aes(x=Family, y=value, fill=factor(Comparison))) +
    	geom_bar(position='dodge',stat='identity', width=0.8) +
      #scale_y_continuous(labels = percent_format()) +
      scale_y_continuous() +
      coord_flip() +
      xlab("") + # Set x-axis label
      ylab("Number of TFs") +
      theme_bw() +
      ggtitle(title) +
      scale_x_discrete(breaks=combinedNames, labels=combinedNames) +
      #theme(legend.position = c(.700, .250), legend.background = element_rect(fill = "transparent"), legend.text.align=0) + 
      theme(legend.position = c(.800, .900), legend.background = element_rect(fill = "transparent"), legend.text.align=0) +
      labs(fill = NULL) +
      # colorA and colorB are reversed because this is a horizontal plot
      scale_fill_manual(breaks=c("A", "B"),values=c(colorB,colorA),labels = c(labA, labB))
    
    # For publications exported EPS at 400x511 resolution (??) 
  }
}