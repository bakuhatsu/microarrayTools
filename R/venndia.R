##############################
## Sven Nelson              ##
## 4/23/2012                ##
## Function: venndia        ##
##############################

# Also includes a simple function for drawing circles, called: circle
#' circle
#'
#' A function to draw circles.
#' @param x x coordinate
#' @param y y coordinate
#' @param r radius
#' @param ...
#'
#' @return A graphic of a circle
#'
#' @examples
#' # No example
circle <- function(x, y, r, ...) {
  ang <- seq(0, 2*pi, length = 100)
  xx <- x + r * cos(ang)
  yy <- y + r * sin(ang)
  polygon(xx, yy, ...)
}

## venndia takes two objects (A and B) each with up- and down-regulated genelists
## To make a list containing two lists: 1) up-regulated genes, and 2) down-regulated genes
#GeneSet <- list()
#GeneSet$UP <- GeneSet_UP # A list of upregulated genes (AGI identifiers all caps)
#GeneSet$DN <- GeneSet_DOWN # A list of down-regulated genes (AGI identifiers all caps)

## venndia is currently case-sensitive, so please provide list in all capitals
#' venndia
#'
#' venndia takes two objects (A and B) each with up- and down-regulated genelists. To make a list containing two lists: 1) up-regulated genes, and 2) down-regulated genes. venndia is currently case-sensitive, so please provide list in all capitals.
#' @param A A geneset object, which is a list containing two vectors $UP and $DN which are lists of up-regulated and downregulated genes in a comparison.
#' @param B Another geneset object for identifying overlap with geneset A.
#' @param titleA A title for the left circle (dataset A)
#' @param titleB A title for the right cirlce (dataset B)
#' @param getdata Set to \code{TRUE} to return a geneset object for A-only, AB-overlap, and B-only.  (Use $A, $AB, or $B after the function call to specify which and you can feed this into another venndia call)
#' @param getcounts Set to \code{TRUE} to return a list the number of genes up- or down-regulated for A-only, AB-overlap, and B-only.  Use $A, $AB, or $B after the function call to specify which if desired. \code{getdata} overrides \code{getcounts} (you can't do both)
#' @param border Set \code{TRUE} to add a black border around circles.
#' @param NoAandB Set \code{TRUE} to remove the "A" and "B" from the plot.
#' @param highLow Choose the font color scheme for numbers of up- and down-regulate genes: \code{redGreen} or \code{pinkBlue}
#' @param noPlot Set \code{FALSE} to prevent plot from displaying (useful when getdata is set to \code{TRUE} and geneset is passed to another venndia call).
#' @param font Select font size (default is 1).
#' @param backgrd Select the background colors for circles by passing a vector of length 2 containing two integers from 1 to 12.  Recommended: default = c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12).
#' @param ... Possible other parameters.
#'
#' @return A Venn diagram visualization and/or a geneset object containg lists of up- and down-regulated genes.
#' @export
#'
#' @examples
#' # TBD
venndia <- function(A, B, titleA="", titleB="", getdata=FALSE, getcounts=FALSE, border=NA, NoAandB = FALSE, highLow = "redGreen", noPlot=FALSE, font=1, backgrd=c(1,2), ...){
  # Alternative UP/DOWN-regulation colors (better for color blind): highLow = "pinkBlue"
  ## Background colors of circles: backgrd = c(left circle color, right circle color)
  # 1 = bluish, 2 = orange, 3 = blue, 4 = red, 5 = green, 6 = yellow

  attach(A)
  AUP <- UP
  ADN <- DN
  detach(A)

  attach(B)
  BUP <- UP
  BDN <- DN
  detach(B)

  # UP
  unionAB_UP <- union(AUP, BUP)
  uniqueA_UP <- setdiff(AUP, BUP) # A but not B
  uniqueB_UP <- setdiff(BUP, AUP) # B but not A
  # Order matters only for intersect (Larger dataset must be first)
  # This is no longer the case after a recent update, but left this code as it does no harm
  if (length(BUP)>length(AUP)) {
    intersAB_UP <- intersect(BUP, AUP)
  } else {
    intersAB_UP <- intersect(AUP, BUP)
  }

  nA_UP <- length(uniqueA_UP)
  nB_UP <- length(uniqueB_UP)
  nAB_UP <- length(intersAB_UP)

  # DOWN
  unionAB_DN <- union(ADN, BDN)
  uniqueA_DN <- setdiff(ADN, BDN)
  uniqueB_DN <- setdiff(BDN, ADN)
  # Order matters only for intersect (Larger dataset must be first)
  # This is no longer the case after a recent update, but left this code as it does no harm
  if (length(BDN)>length(ADN)) {
    intersAB_DN <- intersect(BDN, ADN)
  } else {
    intersAB_DN <- intersect(ADN, BDN)
  }

  nA_DN <- length(uniqueA_DN)
  nB_DN <- length(uniqueB_DN)
  nAB_DN <- length(intersAB_DN)

  if (!noPlot) {
    par(mar=c(2, 2, 2, 2))
    plot(-10, -10, ylim=c(0, 9), xlim=c(0, 9), axes=FALSE, ...)
    circleBackgrounds <- c(
      c(rgb(0.2,1,1,.5), rgb(1,.67,0,.5)), # left = blueish, right = orange
      c(rgb(0,0,1,.5), rgb(1,0,0,.5)), # blue, red
      c(rgb(0,.5,.1,.5), rgb(1,1,0,.5)), # green, yellow
      c(rgb(0.6,1,1,.5), rgb(0.6,1,0.8,.5)), # light blueish,
      c(rgb(0.2,1,0.6,.5), rgb(1,0.86,0.46,.5)), # , light orangish
      c(rgb(1,0,1,.5), rgb(0.961,1,0.98,.5)) #
    )
    circle(x=3, y=4.5, r=3, col=circleBackgrounds[backgrd[1]], border=border) #left color
    circle(x=6, y=4.5, r=3, col=circleBackgrounds[backgrd[2]], border=border) #right color orange
    text( x=c(2.825,6.075), y=c(8,8), c(titleA,titleB), cex=2, col="black" )
    if (NoAandB == FALSE) {
      text( x=c(1.2, 7.7), y=c(6.3, 6.3), c("A", "B"), cex=3, col="gray100" )
    }

    # Set colos scheme for up and downregulated gene numbers
    if (highLow == "redGreen"){
      colors <- c("red","green4")
    } else if (highLow == "pinkBlue") {
      colors <- c("#E82A76","#3B2DD6")
    }

    # UP text
    text(
      x=c(2, 7, 4.5),
      y=c(5.5, 5.5 , 5.5),
      c(nA_UP, nB_UP, nAB_UP),
      cex=2,
      col=colors[1],
      font=font
    )
    # DOWN text
    text(
      x=c(2, 7, 4.5),
      y=c(3.5, 3.5 , 3.5),
      c(nA_DN, nB_DN, nAB_DN),
      cex=2,
      col=colors[2],
      font=font
    )
    par(mar=c(5.1,4.1,4.1,2.1))
  }

  if(getdata){
    list(A=list(UP=uniqueA_UP, DN=uniqueA_DN), B=list(UP=uniqueB_UP, DN=uniqueB_DN),
         AB=list(UP=intersAB_UP, DN=intersAB_DN)
    )
  } else if(getcounts){ # getdata overrides getcounts (you can't do both)
    list(A=list(UP=nA_UP, DN=nA_DN), B=list(UP=nB_UP, DN=nB_DN),
         AB=list(UP=nAB_UP, DN=nAB_DN)
    )
  }
}
