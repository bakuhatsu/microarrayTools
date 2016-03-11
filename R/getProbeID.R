##############################
## Sven Nelson              ##
## 1/18/2013                ##
## Function: getProbeID     ##
##############################

# Input is a gene name: short name (GID1b), AT number (At5g60100), or probe_id
# Input is case-insensitive unless you input the probe_ids (must end in "_at")
# Output is the probe_id for that gene on ATH1 Arabidopsis microarray

#' getProbeID
#'
#' Input is a gene name: short name (GID1b), AT number (At5g60100), or probe_id. Input is case-insensitive unless you input the probe_ids (must end in "_at"). Output is the probe_id for that gene on ATH1 Arabidopsis microarray.
#' @param gene A string or vector of strings specifying a gene or list of genes to convert to Affymetrix probe_id format.  Can take AGI identifiers (eg At1g12345), Affymetrix probe_ids (end in _at), or short names such as ACT1, UBQ11 (not recommended since multiple genes can share the same short name).
#'
#' @return A string or vector of strings representing the Affymetrix probe_id unique identifiers of the gene(s) input (where multiple probe_ids exist for a single gene, the first will be returned)
#' @export
#'
#' @examples
#' # Get the probe_ids of single genes
#' getProbeID("At5g60100")
#' getProbeID("SLY1")
#' getProbeID("245790_at")
#'
#' # Get probe_ids from a list of genes
#' x <- c("AT4G25420", "AT5G51810", "AT5G07200", "AT1G30040","AT4G18350", "AT4G19170", "AT3G24220", "AT1G78390")
#' getProbeID(x)
getProbeID <- function(gene) { # Can take a single gene or list of genes
  require(ath1121501.db)
  for(i in 1:length(gene)) {
    tryCatch(
      if (grepl("[Aa][Tt][012345MmCc][GgTt][0123456789Ee][0123456789]", gene[i])) {
        # make all uppercase
        gene[i] <- toupper(gene[i])
        # Change AT name into probe_id (uses first if mapped to multiple probe_ids)
        gene[i] <- AnnotationDbi::get(gene[i], revmap(ath1121501.db::ath1121501ACCNUM))[1]
      } else if (!grepl("_at", gene[i])) {
        # make all uppercase
        gene[i] <- toupper(gene[i])
        # Change short name to probe_id (uses first if mapped to multiple probe_ids)
        gene[i] <- AnnotationDbi::get(gene[i], revmap(ath1121501.db::ath1121501SYMBOL))[1]
      },
      error=function(e) NULL # Allows direct passthrough of NA or unknown values
    )
  }
  return(gene)
}
