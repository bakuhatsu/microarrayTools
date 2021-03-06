% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venndia.R
\name{venndia}
\alias{venndia}
\title{venndia}
\usage{
venndia(A, B, titleA = "", titleB = "", getdata = FALSE,
  getcounts = FALSE, border = NA, NoAandB = FALSE, highLow = "redGreen",
  noPlot = FALSE, font = 1, backgrd = c(1, 2), ...)
}
\arguments{
\item{A}{A geneset object, which is a list containing two vectors $UP and $DN which are lists of up-regulated and downregulated genes in a comparison.}

\item{B}{Another geneset object for identifying overlap with geneset A.}

\item{titleA}{A title for the left circle (dataset A)}

\item{titleB}{A title for the right cirlce (dataset B)}

\item{getdata}{Set to \code{TRUE} to return a geneset object for A-only, AB-overlap, and B-only.  (Use $A, $AB, or $B after the function call to specify which and you can feed this into another venndia call)}

\item{getcounts}{Set to \code{TRUE} to return a list the number of genes up- or down-regulated for A-only, AB-overlap, and B-only.  Use $A, $AB, or $B after the function call to specify which if desired. \code{getdata} overrides \code{getcounts} (you can't do both)}

\item{border}{Set \code{"black"} to add a black border around circles.}

\item{NoAandB}{Set \code{TRUE} to remove the "A" and "B" from the plot.}

\item{highLow}{Choose the font color scheme for numbers of up- and down-regulate genes: \code{redGreen} or \code{pinkBlue}}

\item{noPlot}{Set \code{FALSE} to prevent plot from displaying (useful when getdata is set to \code{TRUE} and geneset is passed to another venndia call).}

\item{font}{Select font size (default is 1).}

\item{backgrd}{Select the background colors for circles by passing a vector of length 2 containing two integers from 1 to 12.  Recommended: default = c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12).}

\item{...}{Possible other parameters.}
}
\value{
A Venn diagram visualization and/or a geneset object containg lists of up- and down-regulated genes.
}
\description{
venndia takes two objects (A and B) each with up- and down-regulated genelists. To make a list containing two lists: 1) up-regulated genes, and 2) down-regulated genes. venndia is currently case-sensitive, so please provide list in all capitals.
}
\examples{
## For the two genesets to compare:
#' # Example data is provided, use the following code to load it
data("GeneSetdata") # loads GeneSet_UP and GeneSet_DOWN objects which contain example data.

## Make a list containing two vectors: 1) up-regulated genes, and 2) down-regulated gene
## Create geneset A
GeneSetA <- list()
GeneSetA$UP <- GeneSet_UP[sample(1:200, 100)] # Random set of 100 genes up
GeneSetA$DN <- GeneSet_DOWN[sample(1:200, 100)] # Random set of 100 genes down

GeneSetB <- list()
GeneSetB$UP <- GeneSet_UP[sample(1:200, 100)] # Random set of 100 genes up
GeneSetB$DN <- GeneSet_DOWN[sample(1:200, 100)] # Random set of 100 genes down

## Test the overlap
venndia(GeneSetA, GeneSetB) # Outputs a graphic

## Get the results of the overlap between A and B
venndia(GeneSetA,GeneSetB, getdata = TRUE)$AB

## Compare this overlap with a third geneset (GeneSetC).
GeneSetC <- list()
GeneSetC$UP <- GeneSet_UP[sample(1:200, 100)] # Random set of 100 genes up
GeneSetC$DN <- GeneSet_DOWN[sample(1:200, 100)] # Random set of 100 genes down

overlapAB <- venndia(GeneSetA,GeneSetB, getdata = TRUE)$AB

venndia(overlapAB, GeneSetC)

## Prettier version
venndia(overlapAB, GeneSetC, NoAandB = TRUE, border = "black")
}
