#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a BrainGene, that is, a TrenaGene with standard gene expression data
#'
#' @include TrenaGene.R
#' @import methods
#'
#' @name BrainGene
#'

.BrainGene <- setClass("BrainGene", contains="TrenaGene")
#------------------------------------------------------------------------------------------------------------------------
BrainGene <- function(geneSymbol)
{
   expressionData <- list()

} # constructor
#------------------------------------------------------------------------------------------------------------------------
