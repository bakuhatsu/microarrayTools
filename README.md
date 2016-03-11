# microarrayTools
Tools for data analysis of microarray and other transcriptome datasets

### Installing and loading the `microarrayTools` package:
```r
## Install and load devtools for loading packages from GitHub
install.packages("devtools") # to allow us to install packages from GitHub
library(devtools)

## Install microarrayTools package including TAGGITontology, TAGGITplot, and getProbeID
install_github("bakuhatsu/microarrayTools") # user name/library
library(microarrayTools)
```
The microarrayTools package gives you the TAGGITontology, TAGGITplot, getProbeID, and venndia functions

### Setting up lists of up-/down-regulated genes for `TAGGITontology()`:
```r
#### Comparing up-regulated and down-regulated genes for a single comparison ####
## Create a vector of probe_ids (ex 248961_at) or AGI identifiers (ex AT5G45650)
# UP and DOWN regulated genesets
# NOTE: "etc..." is used to indicate that the list continues, do not use in the real vector.
GeneSet_UP <- c("AT4G25420", "AT5G51810", "AT5G07200", "AT1G30040", etc...)
GeneSet_DOWN <- c("AT4G18350", "AT4G19170", "AT3G24220", "AT1G78390", etc...)

# If you have loaded microarray data into R, you can easily create a list from topTable
# Please see the Bioconductor documentation for instructions for use of topTable
# NOTE: Replace "YourData.rma" and "ComparisonCoef" with appropriate input from your data
GeneSet_UP <- row.names(subset(topTable(YourData.rma, coef=ComprisonCoef, adjust="fdr", 
                      sort.by="B", number=Inf, p.value = 0.05)), logFC>0) # only up-regulated genes
GeneSet_DOWN <- row.names(subset(topTable(YourData.rma, coef=ComprisonCoef, adjust="fdr", 
                      sort.by="B", number=Inf, p.value = 0.05)), logFC<0) # only down-regulated genes
```
### Test run using the provided datasets in `GeneSetdata`:
```r
# Example data is provided, use the following code to load it
data("GeneSetdata") # loads GeneSet_UP and GeneSet_DOWN objects which contain example data.

## Make a list containing two vectors: 1) up-regulated genes, and 2) down-regulated gene
GeneSet <- list()
GeneSet$UP <- GeneSet_UP # A list of upregulated genes
GeneSet$DN <- GeneSet_DOWN # A list of down-regulated genes

## Create TAGGITontology objects 
(may take a few minutes, but the built-in progress bar will keep you informed on the progress)
# Returns dataframe for plotting and outputs an excel sheet of hits to the working directory.
GeneSet_TAGGIT_UP <- TAGGITontology(GeneSet$UP, outputFileName = "TAGGITontologyHits_UP.xlsx")   
GeneSet_TAGGIT_DN <- TAGGITontology(GeneSet$DN), utputFileName = "TAGGITontologyHits_DN.xlsx")

## Plot the results of the TAGGIT analysis using ggplot2 via the TAGGITplot function
# Comparing UP and DOWN regulated genesets
TAGGITplot(GeneSet$UP, GeneSet$DN, GeneSet_TAGGIT_UP, GeneSet_TAGGIT_DN, title = "")
## Export the image: for best results export in EPS (vector) format
## For output like the Figure 4 in the publication, export at 400x511 resolution
```
### Comparing sets of overall differentially regulated genes:
(not separatting up- and down-regulated)
```r
#### Total differentially regulated genes (up and down combined) for two comparisons ####
## Create two lists as before: GeneSet_AvsB and GeneSet_CvsD

## Create TAGGITontology objects (the genesets below are not provided, so the code below is just an example)
GeneSet_TAGGIT_AvsB <- TAGGITontology(GeneSet_AvsB)
GeneSet_TAGGIT_CvsD <- TAGGITontology(GeneSet_CvsD)

## Plot the results of the TAGGIT analysis using ggplot2 via the TAGGITplot function
TAGGITplot(GeneSet_AvsB, GeneSet_CvsD, GeneSet_TAGGIT_AvsB, GeneSet_TAGGIT_CvsD, A = "AvsB", B = "CvsD", title = "")
## Export the image: for best results export in EPS (vector) format
## For output like Supplemental Fig. S3, export at 400x511 resolution
```
