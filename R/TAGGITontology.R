##############################
## Sven Nelson              ##
## 9/15/2014                ##
## Function: TAGGITontology ##
##############################

# Original TAGGIT paper: (Carrera et al., 2007)
# This function is hardcoded for use with Arabidopsis ATH1 microarray data
# To use with data from another organism, you will need make a few edits.

### Requires myAnnot to exist
# require(ath1121501.db)
# myAnnot <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "), stringsAsFactors = FALSE)

# myAnnot[x, ]

# For testing, can use: GOI1_ATs
# [1] "AT4G25420" "AT5G51810" "AT5G07200" "AT1G30040" "AT1G15550" "AT5G15230" "AT4G36930"
# [8] "AT4G18350" "AT4G19170" "AT3G24220" "AT1G78390" "AT4G19230" "AT2G29090" "AT5G45340"

### I set up another dataset for testing using example genesets (see end of file)

#' TAGGITontology
#'
#' Original TAGGIT paper: (Carrera et al., 2007)
#' This function is hardcoded for use with Arabidopsis ATH1 microarray data
#' To use with data from another organism, you will need make a few edits.

#' @param geneList A vector of genes as AGI identiers or Affymetrix probe_ids.
#' @param useSearchTerms Set to \code{FALSE} to prevent search term searches in gene descriptions (faster, but results will be based only on genes specified in \code{taggitAGIs} file)
#' @param outputFileName Path to output excel file containing lists of hits with descriptions.  File is output into the current working directory as \code{"TAGGITontologyHits.xlsx"} by default.  Set to \code{"none"} for no output.
#' @param annotationFile Use default unless you wish to provide your own updated \code{myAnnot} file. See description for instructions for creating this file.
#' @param taggitAGIs Generally use default.  Can also pass the path of your own taggitAGIs.tsv (tab-separated) file defining specific AGI identifiers for each TAGGIT term.
#' @param taggitSearchTerms Generally use default.  Can also pass the path of your own taggitSearchTerms.tsv (tab-separated) file defining TAGGIT term-specific search terms to search for within gene descriptions.
#'
#' @return A dataframe containing lists of hits for each TAGGIT term.  Also can save an excel spreadsheet containing these hits with descriptions separated by TAGGIT term.
#' @export
#'
#' @examples
#' # TBD
TAGGITontology <- function(geneList, useSearchTerms = TRUE, outputFileName = "TAGGITontologyHits.xlsx", annotationFile = "default", taggitAGIs = "default", taggitSearchTerms = "default") { # geneList is a list of AT numbers
  # geneList can also be a list of probe_ids (or even *experimental* a list of short names)
  # outputFileName = "none" to prevent saving excel spreadsheet of hits (faster)
  # initial version will take a geneList of differentially expressed genes and return a dataframe with number of hits per TAGGIT ontology term for easy plotting

  ## Check to see if myAnnot exists
  ## myAnnot is a data frame that allows quick mapping between:
  ## a) probe_id, b) gene name (AT number), c) symbol (abbreviation), and d) description
  ## probe_ids are the rownames and columns are: ACCNUM, SYMBOL, and DESC
  ## ACCNUM = gene name, SYMBOL = abbreviated name, DESC = description

  ## Load myAnnot ##
  if (annotationFile == "default") {
    #load(file = "myAnnot.RData")
    data("myAnnot")
  } else {
    myAnnot = annotationFile
  }

  # Create and populate a dataframe with the AT num, symbol, and description
  geneFrame <- myAnnot[getProbeID(geneList),] # see below for example of organization
  print(geneList[1:4])
  print(getProbeID(geneList)[1:4])
  print(geneFrame[1:3,1])
  #              ACCNUM        SYMBOL  DESC
  #248961_at     AT5G45650     NA      subtilase family protein

  # getProbeIDs gives warnings about multiple probe_ids per AT #.
  # Should be fine, but to get rid of the warnings could just pass probe_ids directly.

  if(taggitAGIs == "default") {
    data("taggitAGIs")
  } else {
    TAGGITguideAGIs <- read.table(taggitAGIs, header=TRUE, sep="\t", stringsAsFactors = F)
  }
  if(!exists("TAGGITguideAGIs")) {
    writeLines("\nUnable to locate 'TAGGITguideAGIs' file.")
  }

  if(taggitSearchTerms == "default") {
    data("taggitSearchTerms")
  } else {
    TAGGITguideSearchTerms <- read.table(taggitSearchTerms, header=TRUE, sep="\t", stringsAsFactors = F) # Strings, not factors
  }
  if(!exists("TAGGITguideSearchTerms")) {
    writeLines("\nUnable to locate 'TAGGITguideSearchTerms' file.")
  }

  # Create progress bar:
  total <- length(TAGGITguideAGIs[1,])
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  # How many: overlap between TAGGITguideAGI[column i] and geneFrame$ACCNUM
  # Order matters only for intersect (Larger dataset must be first)
  # I tested this recently and it seems like this has been fixed so order doesn't matter.
  ontologyCounts <- data.frame(matrix(NA, nrow = 1, ncol = 27))
  colnames(ontologyCounts) <- colnames(TAGGITguideAGIs)
  for(i in 1:length(TAGGITguideAGIs[1,])) { # 1 to 27
    #if (length(BUP)>length(AUP)) {
    searchTermHits <- 0
    AGImatches <- c()
    searchResultProbes <- c()

    # for progress bar
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
    if (length(TAGGITguideAGIs[,i][startsWith(TAGGITguideAGIs[,i], pattern = "A")])>length(geneFrame$ACCNUM)) { #requires gdata for startsWith()
      AGImatches <- intersect(TAGGITguideAGIs[,i], geneFrame$ACCNUM)
    } else {
      AGImatches <- intersect(geneFrame$ACCNUM, TAGGITguideAGIs[,i])
    }
    AGIremaining <- setdiff(geneFrame$ACCNUM, AGImatches)
    newFrame <- geneFrame[getProbeID(AGIremaining),]
    # print(TAGGITguideAGIs[1:4,i]) # ok
    # print(AGImatches[1:3]) # NAs
    # print(geneFrame$ACCNUM[1:3]) # NAs
    # make this into a number of hits
    if (useSearchTerms) {

      searchTerms <- NULL
      if (!(is.na(TAGGITguideSearchTerms[1,i]) | TAGGITguideSearchTerms[1,i]=="")) {
        searchTerms <- TAGGITguideSearchTerms[!(is.na(TAGGITguideSearchTerms[,i]) | TAGGITguideSearchTerms[,i]==""),i]
      }

      if (!is.null(searchTerms)) {
        # newFrame is geneFrame of those not already categorized based on AGI match
        # rowname is probe_id, then ACCNUM, SYMBOL, DESC
        # searchTermHits <- searchDesc(newFrame, searchTerms) # returns searchTermHits
        if (i==27) {
          ## SPECIAL CASE FOR "Unnanotated" must exactly equal search term
          if (outputFileName != "none") {
            searchResultProbes <- c(row.names(newFrame[newFrame$DESC=="expressed protein",]), row.names(newFrame[newFrame$DESC=="hypothetical protein",]), row.names(newFrame[newFrame$DESC=="unknown protein",]))
          } else {
            searchTermHits <- sum(newFrame$DESC == "expressed protein") + sum(newFrame$DESC == "hypothetical protein") + sum(newFrame$DESC == "unknown protein")
          }
        } else if (i!=27) {
          if (class(newFrame[,3])=="character" & class(searchTerms)=="character") {
            if (outputFileName != "none") {
              searchResults <- searchDescReturnHits(newFrame[,3], searchTerms) # C++ function with Rcpp
              #searchResultAGIs <- newFrame[searchResults,1]
              searchResultProbes <- row.names(newFrame[searchResults + 1,]) # + 1 because C++ counts from index 0, whereas R counts from index 1
            } else {
              # Source searchTermHits here?
              #print(newFrame[,3])
              print(searchTerms)
              searchTermHits <- searchDesc(newFrame[,3], searchTerms) # C++ function with Rcpp
              #searchTermHits <- searchDesc(tolower(newFrame[,3]), tolower(searchTerms)) # C++ function with Rcpp
            }
          } else {
            writeLines("\nERROR: Please ensure that searchTerms and newFrame[,3] are character vectors.  Skipping search term query in descriptions.")
          }


        }
      }
    }

    ## Return Hits ##
    if (outputFileName != "none") {
      listOfHits <- c(getProbeID(AGImatches), searchResultProbes)
      listOfHits <- listOfHits[!is.na(listOfHits)]
      if (i == 1) {
        Dormancy.related <- myAnnot[listOfHits,]
      } else if (i == 2) {
        Germination.related <- myAnnot[listOfHits,]
        #print(Germination.related)
      } else if (i == 3) {
        ABA <- myAnnot[listOfHits,]
      } else if (i == 4) {
        Auxin <- myAnnot[listOfHits,]
      } else if (i == 5) {
        Brassinosteroid <- myAnnot[listOfHits,]
      } else if (i == 6) {
        Cytokinin <- myAnnot[listOfHits,]
      } else if (i == 7) {
        Ethylene <- myAnnot[listOfHits,]
      } else if (i == 8) {
        Gibberellin <- myAnnot[listOfHits,]
      } else if (i == 9) {
        Jasmonic.acid <- myAnnot[listOfHits,]
      } else if (i == 10) {
        Seed.storage.proteins.LEAs <- myAnnot[listOfHits,]
      } else if (i == 11) {
        Inhibition.protein.degrad <- myAnnot[listOfHits,]
      } else if (i == 12) {
        Protein.degradation <- myAnnot[listOfHits,]
      } else if (i == 13) {
        Heat.Shock <- myAnnot[listOfHits,]
      } else if (i == 14) {
        Cell.wall.modification <- myAnnot[listOfHits,]
      } else if (i == 15) {
        Cell.cycle.related <- myAnnot[listOfHits,]
      } else if (i == 16) {
        Cytoskeleton <- myAnnot[listOfHits,]
      } else if (i == 17) {
        Translation.associated <- myAnnot[listOfHits,]
      } else if (i == 18) {
        DNA.repair <- myAnnot[listOfHits,]
      } else if (i == 19) {
        Respiration <- myAnnot[listOfHits,]
      } else if (i == 20) {
        Electron.Transport <- myAnnot[listOfHits,]
      } else if (i == 21) {
        Pentose.phosphate.pathway <- myAnnot[listOfHits,]
      } else if (i == 22) {
        Glycolysis.and.gluconeogenesis <- myAnnot[listOfHits,]
      } else if (i == 23) {
        Krebs.cycle <- myAnnot[listOfHits,]
      } else if (i == 24) {
        Beta.oxidation <- myAnnot[listOfHits,]
      } else if (i == 25) {
        Stress <- myAnnot[listOfHits,]
      } else if (i == 26) {
        Photosynthesis.chloroplast <- myAnnot[listOfHits,]
      } else if (i == 27) {
        Unannotated <- myAnnot[listOfHits,]
      }
      ontologyCounts[,i] <- length(AGImatches) + length(searchResultProbes)
    } else {
      ontologyCounts[,i] <- length(AGImatches) + searchTermHits
    }
    if (exists("AGImatches")) {
      rm(AGImatches)
    }
    if (exists("searchResultProbes")) {
      rm(searchResultProbes)
    }
  }
  close(pb) # Progress bar complete
  if (outputFileName != "none") {
    #require(WriteXLS)
    # Write to a multisheet excel file.  Each data frame will make one sheet.
    WriteXLS::WriteXLS(c("Dormancy.related", "Germination.related", "ABA", "Auxin", "Brassinosteroid", "Cytokinin", "Ethylene", "Gibberellin", "Jasmonic.acid", "Seed.storage.proteins.LEAs", "Inhibition.protein.degrad", "Protein.degradation", "Heat.Shock", "Cell.wall.modification", "Cell.cycle.related", "Cytoskeleton", "Translation.associated", "DNA.repair", "Respiration", "Electron.Transport", "Pentose.phosphate.pathway", "Glycolysis.and.gluconeogenesis", "Krebs.cycle", "Beta.oxidation", "Stress", "Photosynthesis.chloroplast", "Unannotated"), ExcelFileName=outputFileName, row.names=FALSE)
    #writeLines("\n")
    writeLines(sprintf("Analysis complete, results have been output to an excel spreadsheet: '%s'", outputFileName)) # outputFileName can include the path
  }
  # returns a data.frame of counts for each TAGGIT category
  return(ontologyCounts)
}

### I set up another dataset for testing using example genesets (see end of file)
# rm(GeneSet)
# rm(GeneSet_UP)
# rm(GeneSet_DOWN)
#
# GeneSet_UP <- ARvsD12hset$UP[1:200]
# GeneSet_DOWN <- ARvsD12hset$DN[1:200]
#
# #save(x, y, file = "xy.RData")
# save(GeneSet_UP, GeneSet_DOWN, file = "GeneSetdata.RData")

# ## To use:
# data("GeneSetdata") # loads GeneSet_UP and GeneSet_DOWN objects which contain example data.
#
# # Make a list containing two lists: 1) up-regulated genes, and 2) down-regulated genes
# GeneSet <- list()
# GeneSet$UP <- GeneSet_UP # A list of upregulated genes
# GeneSet$DN <- GeneSet_DOWN # A list of down-regulated genes
#
# ## Create TAGGITontology objects and save it to a variable
# #  This step may take a long time especially with list of more than 1000 genes
# GeneSet_TAGGIT_UP <- TAGGITontology(GeneSet$UP) # ignore any warning messages
# GeneSet_TAGGIT_DN <- TAGGITontology(GeneSet$DN) # ignore any warning messages
#
# ## Plot the results of the TAGGIT analysis using ggplot2 via the TAGGITplot function
# # Comparing UP and DOWN regulated genesets
# TAGGITplot(GeneSet$UP, GeneSet$DN, GeneSet_TAGGIT_UP, GeneSet_TAGGIT_DN, title = "")
#
# ## Export the image: for best results export in EPS (vector) format
# ## For output like Figure 4, export at 400x511 resolution
#
# TAGGITguideAGIs <- read.table("/Users/sven/Dropbox/00. Code files/TAGGITguideAGIs.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
# TAGGITguideSearchTerms <- read.table("/Users/sven/Dropbox/00. Code files/TAGGITguideSearchTerms.tsv", header=TRUE, sep="\t", stringsAsFactors = F) # Strings, not factors
# #save(x, y, file = "xy.RData")
# save(TAGGITguideAGIs, file = "taggitAGIs.RData")
# save(TAGGITguideSearchTerms, file = "taggitSearchTerms.RData")
