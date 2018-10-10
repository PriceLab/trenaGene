#' @import trena
#' @importFrom DBI   dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#'
#' @title TrenaGene
#------------------------------------------------------------------------------------------------------------------------
#' @name TrenaGene-class
#' @rdname TrenaGene-class
#' @aliases TrenaGene
#'
#' @import methods

.TrenaGene <- setClass ("TrenaGene",
                        representation = representation(
                            geneSymbol="character",
                            genome="character",
                            transcripts="data.frame",
                            expressionDataDirectory="character",
                            quiet="logical"
                            )
                         )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getTranscriptsTable',      signature='obj', function(obj) standardGeneric ('getTranscriptsTable'))
setGeneric('getExpressionMatrixNames', signature='obj', function(obj) standardGeneric ('getExpressionMatrixNames'))
setGeneric('loadExpressionData',       signature='obj', function(obj, datasetName) standardGeneric ('loadExpressionData'))
setGeneric('getEnhancers',             signature='obj', function(obj) standardGeneric ('getEnhancers'))
setGeneric('getEncodeDHS',             signature='obj', function(obj) standardGeneric ('getEncodeDHS'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class Trena
#'
#' @description
#' TrenaGene and its (projected) subclasses provide convenient containers in which to collect
#'  trena-related aggregation of a gene's (a hierarchy of classes) including expression data,
#' transcript and variant info, genomic and epigenomic context, trena models and/or the means to create them
#'
#' @rdname TrenaGene-class
#'
#' @param geneSymbol  A character string in standard HUGO nomenclature
#' @param genomeName A string indicating the genome used by the Trena object.
#'                  Currently, only human and mouse ("hg38","mm10") are supported
#' @parame expressionDataDirectory A string pointing to a collection of RData expression matrices
#' @param quiet A logical indicating whether or not the Trena object should print output
#'
#' @return An object of the TrenaGene class
#'
#' @export
#'
#' @examples
#' # Create a Trena object using the human hg38 genome
#' cola1a <- TrenaGene("COL1A1", "hg38")
#'
TrenaGene <- function(geneSymbol, genomeName, expressionDataDirectory=NA_character_, quiet=TRUE)
{
   stopifnot(genomeName %in% c("hg38"))
   tbl.transcripts <- .getCodingTranscripts(geneSymbol, genomeName)
   expressionData <- list()

   .TrenaGene(geneSymbol=geneSymbol,
              genome=genomeName,
              transcripts=tbl.transcripts,
              expressionDataDirectory=expressionDataDirectory,
              quiet=quiet)

} # ctor
#------------------------------------------------------------------------------------------------------------------------
.getCodingTranscripts <- function(geneSymbol, genomeName)
{
   driver <- RPostgreSQL::PostgreSQL()
   genome.db <- DBI::dbConnect(driver, user= "trena", password="trena", dbname="hg38", host="khaleesi")
   expected.tables <- c("gtf")
   stopifnot(all(expected.tables %in% DBI::dbListTables(genome.db)))
   query.0 <- sprintf("select * from gtf where gene_name='%s'", geneSymbol)
   query.1 <- "and gene_biotype='protein_coding' and moleculetype in ('gene', 'transcript')"
   query <- paste(query.0, query.1)
   tbl.transcripts <- DBI::dbGetQuery(genome.db, query)
   stopifnot(nrow(tbl.transcripts) > 0)

   return(tbl.transcripts)


} # .getCodingTranscripts
#------------------------------------------------------------------------------------------------------------------------
#' Get the transcripts for the gene
#'
#' @rdname getTranscriptsTable
#' @aliases getTranscriptsTable
#'
#' @param obj An object of class TrenaGene
#'
#' @export

setMethod('getTranscriptsTable',  'TrenaGene',

   function(obj) {
        return(obj@transcripts)
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get the names of the expression matrices - their names with directory and .RData suffix stripped out
#'
#' @rdname getExpressionMatrixName
#' @aliases getExpressionMatrixName
#'
#' @param obj An object of class TrenaGene
#'
#' @export


setMethod('getExpressionMatrixNames',  'TrenaGene',

   function(obj) {
      if(is.na(obj@expressionDataDirectory))
         return(list())
      all.files <- list.files(obj@expressionDataDirectory)
      rdata.filenames <- grep(".RData$", all.files, value=TRUE)
      sub(".RData", "", rdata.filenames, fixed=TRUE)
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the a specifc expression matrix
#'
#' @rdname getExpressionMatrixName
#' @aliases getExpressionMatrixName
#'
#' @param obj An object of class TrenaGene
#' @param datasetName A numeric matrix
#'
#' @export

setMethod('loadExpressionData',  'TrenaGene',

    function(obj, datasetName){
       if(datasetName %in% getExpressionMatrixNames(obj)){
          filename <- sprintf("%s.RData", datasetName)
          full.path <- file.path(obj@expressionDataDirectory, filename)
          stopifnot(file.exists(full.path))
          eval(parse(text=paste("mtx <- ", load(full.path))))
          invisible(mtx)
          }
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname getEnhancers
#' @aliases getEnhancers
#'
#' @param obj An object of class TrenaGene
#'
#' @export

setMethod('getEnhancers',  'TrenaGene',

     function(obj){
        full.path <- file.path(obj@expressionDataDirectory, "genomicFeatures", "geneHancer.v4.7.allGenes.RData")
        stopifnot(file.exists(full.path))
        load(full.path)
        subset(tbl.enhancers, geneSymbol == obj@geneSymbol)
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the dnase hypersensitivity regions in the expansive region covered by the enhancer
#'
#' @rdname getEncodeDHS
#' @aliases getEncodeDHS
#'
#' @param obj An object of class TrenaGene
#'
#' @export


setMethod('getEncodeDHS',   'TrenaGene',

    function(obj){
       hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                             geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
       tbl.enhancers <- getEnhancers(obj)
       chrom <- tbl.enhancers$chrom[1]
       loc.min <- min(tbl.enhancers$start)
       loc.max <- max(tbl.enhancers$end)
       tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, loc.min, loc.max)
       return(tbl.dhs)
       })

#------------------------------------------------------------------------------------------------------------------------


