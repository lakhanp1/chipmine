


#' Transcripts GRanges object from TxDB
#'
#' @param txdb \code{TxDB} object which will be used for annotation
#' @param excludeType Types of transcripts to exclude from annotation. Should be a
#' character vector. Default: \code{c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA")}.
#' This argument work only when \code{TxDB} object has \code{TXTYPE} column with
#' appropriate transcripy type values.
#' @param txIds A vector of transcript IDs to be used specifically in the annotation
#' process instead of full transcript set. These should be internal tx_ids from \code{TxDB}
#' object. This is useful feature to exclude tRNA, rRNA transcripts while annotating
#' the regions. Default: NULL
#'
#' @return GRanges object from TxDB
#'
#' @examples NA
get_txdb_transcripts_gr <- function(txdb, excludeType = NULL, txIds = NULL){

  if(exists("transcriptsGr", envir=txdbEnv, inherits=FALSE)) {

    transcriptsGr <- get(x = "transcriptsGr", envir = txdbEnv)

  } else{

    ## transcript to gene map
    txToGene <- suppressMessages(
      AnnotationDbi::select(
        x = txdb, keys = AnnotationDbi::keys(x = txdb, keytype = "TXID"),
        columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
      dplyr::mutate(TXID = as.character(TXID)) %>%
      dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

    assign(x = "txToGene", value = txToGene, envir = txdbEnv)

    ## decide which transcript types to select for the annotation
    allTxTypes <- unique(txToGene$txType)
    selectType <- allTxTypes[which(!allTxTypes %in% excludeType)]

    ## txIds filter list
    txFilter <- NULL
    txFilter$tx_type <- selectType
    if(!is.null(txIds)){
      txFilter$tx_id <- txIds
    }

    ## extract transcript GRanges
    transcriptsGr <- GenomicFeatures::transcripts(
      x = txdb, filter = txFilter,
      columns = c("tx_id", "tx_name", "tx_type", "gene_id")
    )

    mcols(transcriptsGr)$tx_id <- as.character(mcols(transcriptsGr)$tx_id)

    assign(x = "transcriptsGr", value = transcriptsGr, envir = txdbEnv)

  }

  return(transcriptsGr)
}


##################################################################################


#' 5UTR GRanges object from TxDB
#'
#' @inheritParams get_txdb_transcripts_gr
#'
#' @return GRanges object from TxDB
#'
#' @examples NA
get_txdb_fiveUtr_gr <- function(txdb, txIds = NULL){

  if(exists("fiveUtrGr", envir=txdbEnv, inherits=FALSE)) {
    fiveUtrGr <- get(x = "fiveUtrGr", envir = txdbEnv)

  } else{
    fiveUtrGrl <- GenomicFeatures::fiveUTRsByTranscript(txdb)

    if(!is.null(txIds)){
      fiveUtrGrl <- fiveUtrGrl[names(fiveUtrGrl) %in% txIds]
    }

    ## remember to combine the multi-exon UTRs from UTR GRangesList
    fiveUtrGr <- unlist(range(fiveUtrGrl))
    mcols(fiveUtrGr)$tx_id <- names(fiveUtrGr)
    assign(x = "fiveUtrGr", value = fiveUtrGr, envir = txdbEnv)
  }

  return(fiveUtrGr)
}

##################################################################################

#' 3UTR GRanges object from TxDB
#'
#' @inheritParams get_txdb_transcripts_gr
#'
#' @return GRanges object from TxDB
#'
#' @examples NA
get_txdb_threeUtr_gr <- function(txdb, txIds = NULL){
  if(exists("threeUtrGr", envir=txdbEnv, inherits=FALSE)) {
    threeUtrGr <- get(x = "threeUtrGr", envir = txdbEnv)

  } else{
    threeUtrGrl <- GenomicFeatures::threeUTRsByTranscript(txdb)

    if(!is.null(txIds)){
      threeUtrGrl <- threeUtrGrl[names(threeUtrGrl) %in% txIds]
    }

    ## remember to combine the multi-exon UTRs from UTR GRangesList
    threeUtrGr <- unlist(range(threeUtrGrl))
    mcols(threeUtrGr)$tx_id <- names(threeUtrGr)
    assign(x = "threeUtrGr", value = threeUtrGr, envir = txdbEnv)
  }

  return(threeUtrGr)
}


##################################################################################


#' Exons GRanges object from TxDB
#'
#' @inheritParams get_txdb_transcripts_gr
#'
#' @return GRanges object from TxDB
#'
#' @examples NA
get_txdb_exons_gr <- function(txdb, txIds = NULL){
  if(exists("exonsGr", envir=txdbEnv, inherits=FALSE)) {
    exonsGr <- get(x = "exonsGr", envir = txdbEnv)
  } else{
    exonsGrl <- GenomicFeatures::exonsBy(x = txdb, by = "tx")

    if(!is.null(txIds)){
      exonsGrl <- exonsGrl[names(exonsGrl) %in% txIds]
    }

    exonsGr <- unlist(exonsGrl)
    mcols(exonsGr)$tx_id <- names(exonsGr)
    assign(x = "exonsGr", value = exonsGr, envir = txdbEnv)
  }

  return(exonsGr)
}


##################################################################################


#' Introns GRanges object from TxDB
#'
#' @inheritParams get_txdb_transcripts_gr
#'
#' @return GRanges object from TxDB
#'
#' @examples NA
get_txdb_introns_gr <- function(txdb, txIds = NULL){
  if(exists("intronsGr", envir=txdbEnv, inherits=FALSE)) {
    intronsGr <- get(x = "intronsGr", envir = txdbEnv)
  } else{
    intronsGrl <- GenomicFeatures::intronsByTranscript(x = txdb)

    if(!is.null(txIds)){
      intronsGrl <- intronsGrl[names(intronsGrl) %in% txIds]
    }

    intronsGr <- unlist(intronsGrl)
    mcols(intronsGr)$tx_id <- names(intronsGr)
    assign(x = "intronsGr", value = intronsGr, envir = txdbEnv)
  }

  return(intronsGr)
}


##################################################################################


