#library(roxygen2); # Read in the roxygen2 R package uncomment these two lines to modify the help files for the package
#roxygenise();      # Builds the help files
#' MicroRNA profiles of 84 colon adenocarcinomas 
#'
#' Data from a non-coding RNA profiling by array. RNA extraction from pairs of colon adenocarcinoma tissue and adjacent nontumorous tissue.
#' RNA from pairs of tissues were hybridized on the same day with up to 12 pairs of tissues being hybridized on a given day. 
#' To identify differentially expressed microRNAs, paired class comparison in BRB array tools was used. 
#' Samples were paired in the individual from where the tissues originated. 
#' Therefore, every tumor tissue has it's own nontumorous reference for comparison.
#' Each pair of columns contains the nontumorous and tumor tissue, respectively for a biomarker. Each of the 84 rows is from a different individual.
#' 
#'
#' @docType data
#'
#' @usage data(colondat)
#'
#' @format An object of class \code{"data.frame"} with 84 rows and 1510 variables
#' \describe{
#'   \item{biomarkernamen}{nontumorous tissue for the respective biomarker}
#'   \item{biomarkernameT}{tumor tissue for the respective biomarker}
#'   }
#' @keywords datasets
#'
#' @references Schetter et al. (2008) JAMA 299:425-1036
#' (\href{https://pubmed.ncbi.nlm.nih.gov/18230780/}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7828}{Gene Expression Omnibus}
#'
#' @examples
#' data(colondat)
#' araB1b <- data.frame(nontumorous=colondat$`araB 1bn`,tumor=colondat$`araB 1bT`) # nontumorous and tumor tissue partially matched pair samples for microRNA arab1b
#' p.values.pooling <- (araB1b$nontumorous,araB1b$tumor,method = "liptak",alternative = "greater")
"colondat"
