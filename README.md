# microarrayTools
Tools for data analysis of microarray and other transcriptome datasets

While this package will work with the default Rgui application that comes with R, it is highly recommended to use RStudio for data analysis in R.  Using RStudio will make things like saving output to vector images (like EPS), tiff, or PDF, effortless, and provide code completion, as well as a wealth of other functionalities.  RStudio is freely available for Windows, Mac, or Linux from [https://www.rstudio.com](https://www.rstudio.com). 

### Installing and loading the `microarrayTools` package:
Installing through GitHub requires the `devtools` package, so instructions for installing and loading that package are provided below.
```r
## Install and load devtools for loading packages from GitHub
install.packages("devtools") # to allow us to install packages from GitHub
library(devtools)
```
Since GitHub packaged are compiled on your machine to run, you may be prompted to "install build tools" or something similar.  Follow the instructions, which will install the tools needed to compile this package.

`microarrayTools` also relies on the `ath1121501.db` package from Bioconductor. Bioconductor packages cannot be automatically installed like other R dependencies, so you need to install `ath1121501.db` manually.  Install via the code below (if you do not already have bioconductor the first install takes a while, so be prepared).
```r
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## If you haven't previously installed Bioconductor, run the next line.
biocLite() # Installs Bioconductor base packages, will take a long time for a fresh install.  
## To install the ath1121501.db package
biocLite("ath1121501.db")
```
Now you should be ready to install and then load the `microarrayTools` package
```r
## Install microarrayTools package including TAGGITontology, TAGGITplot, getProbeID, and venndia
install_github("bakuhatsu/microarrayTools") # syntax for installing from GitHub: username/library
library(microarrayTools) # To load the package
```
The microarrayTools package gives you the TAGGITontology, TAGGITplot, getProbeID, and venndia functions.  TAGGITontology also makes use of a few custom C++ functions (requiring Rcpp) to dramatically speed up the analysis.

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

## Make a list containing two vectors: 1) up-regulated genes, and 2) down-regulated genes
GeneSet <- list()
GeneSet$UP <- GeneSet_UP # A vector of upregulated genes
GeneSet$DN <- GeneSet_DOWN # A vector of down-regulated genes

## Create TAGGITontology objects 
# (may take a few minutes, but the built-in progress bar will keep you informed on the progress)
# Returns dataframe for plotting and outputs an excel sheet of hits to the working directory.
GeneSet_TAGGIT_UP <- TAGGITontology(GeneSet$UP, outputFileName = "TAGGITontologyHits_UP.xlsx")   
GeneSet_TAGGIT_DN <- TAGGITontology(GeneSet$DN, outputFileName = "TAGGITontologyHits_DN.xlsx")

## Plot the results of the TAGGIT analysis using ggplot2 via the TAGGITplot function
# Comparing UP and DOWN regulated genesets
TAGGITplot(GeneSet$UP, GeneSet$DN, GeneSet_TAGGIT_UP, GeneSet_TAGGIT_DN, title = "")
## Export the image: for best results export in EPS (vector) format
## For output like Figure 4 in the publication, export at 400x511 resolution
```
### Comparing sets of overall differentially regulated genes:
(not separating up- and down-regulated)
```r
#### Total differentially regulated genes (up and down combined) for two comparisons ####
## Create two lists as before: GeneSet_AvsB and GeneSet_CvsD

## Create TAGGITontology objects (the genesets below are not provided, so the code below is just an example)
GeneSet_TAGGIT_AvsB <- TAGGITontology(GeneSet_AvsB)
GeneSet_TAGGIT_CvsD <- TAGGITontology(GeneSet_CvsD)

## Plot the results of the TAGGIT analysis using ggplot2 via the TAGGITplot function
TAGGITplot(GeneSet_AvsB, GeneSet_CvsD, GeneSet_TAGGIT_AvsB, GeneSet_TAGGIT_CvsD, A = "AvsB", B = "CvsD", 
                      title = "")
## Export the image: for best results export in EPS (vector) format
## For output like Supplemental Fig. S3, export at 400x511 resolution
```
