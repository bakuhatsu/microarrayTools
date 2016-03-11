##############################
## Sven Nelson              ##
## 11/14/2014               ##
## Function: TAGGITplot     ##
##############################

# Original TAGGIT paper: (Carrera et al., 2007)

### Takes a TAGGITontology object and plots it using ggplot2

#' TAGGITplot
#'
#' Takes a TAGGITontology object and plots it using ggplot2
#' @param geneListA A vector of genes as AGI identiers or Affymetrix probe_ids.
#' @param geneListB Another vector of genes as AGI identiers or Affymetrix probe_ids.
#' @param ontologyCountsA The dataframe resulting from a \code{TAGGITontology()} call on \code{geneListA}.
#' @param ontologyCountsB The dataframe resulting from a \code{TAGGITontology()} call on \code{geneListB}.
#' @param A UPreg by default (colors will be red/blue), change if not comparing up- to down-regulated geneLists and colors will also change.
#' @param B DOWNreg by default (colors will be red/blue), change if not comparing up- to down-regulated geneLists and colors will also change.
#' @param labA Label for A (defaults to whatever was passed to \code{A}).
#' @param labB Label for B (defaults to whatever was passed to \code{B}).
#' @param title Default is \code{"TAGGIT ontology"}.  Set to \code{""} for no title.
#'
#' @return A plot
#' @export
#'
#' @examples
#' # TBD
TAGGITplot <- function(geneListA, geneListB, ontologyCountsA, ontologyCountsB, A = "UPreg", B = "DOWNreg", labA = A, labB = B, title = "TAGGIT ontology") {
  ## Example code: UP and DOWN
  # TAGGITplot(ARvsD12hset$UP, ARvsD12hset$DN, ARvsD12h_UP_TAGGIT, ARvsD12h_DN_TAGGIT)
  ## Example code: diffA and diffB
  # TAGGITplot(c(ARvsD12hset$UP,ARvsD12hset$DN), c(WTvsD12hset$UP,WTvsD12hset$DN), ARvsD12h_TAGGIT, WTvsD12h_TAGGIT,A="ARvsD12h",B="WTvsD12h")

  #### Preparing data ####
  # Use TAGGITontology to define ontologyCounts:
  # ontologyCountsA => UP, ontologyCountsB => DOWN

  # Create dataframe with data to plot: convert counts into percentage
  TAGGITcluster.df <- c()

  TAGGITcluster.df$Ontology <- names(ontologyCountsA)
  TAGGITcluster.df$A <- unlist(ontologyCountsA)/length(geneListA)
  TAGGITcluster.df$B <- unlist(ontologyCountsB)/length(geneListB)


  TAGGITcluster.df <- as.data.frame(TAGGITcluster.df)

  #require(reshape) # for melt (actually melt.data.frame, it seems)

  TAGGITcluster.long <- reshape::melt.data.frame(TAGGITcluster.df,
                             ## ID variables:
                             # variables to keep but not split apart on
                             id.vars="Ontology",
                             # Measure variables: the source columns
                             measure.vars=c("B","A"),
                             # Name of the destination column that
                             # will identify the original
                             # column that the measurement came from
                             variable_name="Comparison"
  )

  ## Code for testing: Prints out the order of bars currently ##
  #levels(TAGGITcluster12hARvsDandWTvsDp.long$Ontology)

  # Reorder the data by Ontology
  TAGGITcluster.long$Ontology <- factor(TAGGITcluster.long$Ontology, levels = c("Unannotated",  "Photosynthesis.chloroplast.related", "Stress", "Beta.oxidation", "Krebs.cycle",  "Glycolysis.and.gluconeogenesis", "Pentose.phosphate.pathway", "Electron.Transport",  "Respiration",  "DNA.repair", "Translation.associated", "Cytoskeleton", "Cell.cycle.related", "Cell.wall.modification", "Heat.Shock", "Protein.degradation", "Inhibition.of.protein.degradation", "Seed.storage.proteins.Late.Embryogenesis.Abundant", "Jasmonic.acid", "Gibberellin", "Ethylene", "Cytokinin", "Brassinosteroid", "Auxin", "ABA",  "Germination.related",  "Dormancy.related"))

  #levels(TAGGITcluster0hand12hARvsDp.long$Ontology) # Prints order again

  #require(scales) # For percent_format() in plot
  #require(ggplot2)

  if(A == "UPreg" & B == "DOWNreg") { # UP => pinkish and DOWN => blueish
    colorA <- "#E82A76" # UP => pinkish
    colorB <- "#3B2DD6" # DOWN => blueish
  } else {
    colorA <- "#E69F00" # A => yellow
    colorB <- "#999999" # B => gray
  }

  # To remove the "Unnanotated" category, includes genes annotated as "unknown protein"
  TAGGITcluster.long <- TAGGITcluster.long[TAGGITcluster.long$Ontology!="Unannotated",]

  #### Now for the TAGGIT plot ####
  ggplot2::ggplot(data=TAGGITcluster.long, ggplot2::aes(x=Ontology, y=value, fill=factor(Comparison))) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 7.5, xmax = 8.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 13.5, xmax = 14.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 15.5, xmax = 16.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 17.5, xmax = 18.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 19.5, xmax = 20.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 21.5, xmax = 22.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 23.5, xmax = 24.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_rect(ggplot2::aes(xmin = 25.5, xmax = 26.5, ymin = -Inf, ymax = Inf), fill = "gray90") +
    ggplot2::geom_bar(position='dodge',stat='identity', width=0.8) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::coord_flip() +
    ggplot2::xlab("") + # Set x-axis label
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_discrete(breaks=c("Dormancy.related", "Germination.related", "ABA",  "Auxin",  "Brassinosteroid",  "Cytokinin",  "Ethylene",  "Gibberellin",  "Jasmonic.acid",  "Seed.storage.proteins.Late.Embryogenesis.Abundant",  "Inhibition.of.protein.degradation",	"Protein.degradation",	"Heat.Shock",	"Cell.wall.modification",	"Cell.cycle.related",	"Cytoskeleton",	"Translation.associated",	"DNA.repair",	"Respiration",	"Electron.Transport",	"Pentose.phosphate.pathway",	"Glycolysis.and.gluconeogenesis",	"Krebs.cycle",	"Beta.oxidation",	"Stress",	"Photosynthesis.chloroplast.related"), labels=c("Dormancy related", "Germination related", "ABA", "Auxin", "Brassinosteroid", "Cytokinin", "Ethylene", "Gibberellin", "Jasmonic acid", "LEAs/seed storage prot.", "Inhib. protein degradation", "Protein degradation", "Heat Shock", "Cell-wall modification", "Cell cycle related", "Cytoskeleton", "Translation associated", "DNA repair", "Respiration", "Electron Transport", "Pentose phos. pathway", "Glycol. & gluconeogen.", "Krebs cycle", "Beta oxidation", "Stress", "Photosynthesis")) + # removed "Unannotated"
    ggplot2::theme(legend.position = c(.700, .250), legend.background = ggplot2::element_rect(fill = "transparent"), legend.text.align=0) +
    ggplot2::labs(fill = NULL) +
    # colorA and colorB are reversed because this is a horizontal plot
    ggplot2::scale_fill_manual(breaks=c("A", "B"),values=c(colorB,colorA),labels = c(labA, labB))

  # For publications exported EPS at 400x511 resolution
}
