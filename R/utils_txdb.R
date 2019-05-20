


get_txdb_transcripts_gr <- function(txdb, excludeType = NULL, tx = NULL){

  if(exists("transcriptsGr", envir=txdbEnv, inherits=FALSE)) {

    transcriptsGr <- get(x = "transcriptsGr", envir = txdbEnv)

  } else{

    ## transcript to gene map
    txToGene <- suppressMessages(
      AnnotationDbi::select(x = txdb, keys = AnnotationDbi::keys(x = txdb, keytype = "TXID"),
                            columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
      dplyr::mutate(TXID = as.character(TXID)) %>%
      dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

    ## decide which transcript types to select for the annotation
    allTxTypes <- unique(txToGene$txType)
    selectType <- allTxTypes[which(!allTxTypes %in% excludeType)]

    ## tx filter list
    txFilter <- NULL
    txFilter$tx_type <- selectType
    if(!is.null(tx)){
      txFilter$tx_id <- tx
    }

    ## extract transcript GRanges
    transcriptsGr <- GenomicFeatures::transcripts(x = txdb,
                                                  columns = c("tx_id", "tx_name", "tx_type"),
                                                  filter = txFilter)

    assign(x = "transcriptsGr", value = transcriptsGr, envir = txdbEnv)

  }

  return(transcriptsGr)
}


##################################################################################


get_txdb_fiveUtr_grl <- function(txdb){

  if(exists("fiveUtrGrl", envir=txdbEnv, inherits=FALSE)) {
    fiveUtrGrl <- get(x = "fiveUtrGrl", envir = txdbEnv)
  } else{
    fiveUtrGrl <- GenomicFeatures::fiveUTRsByTranscript(txdb)
    assign(x = "fiveUtrGrl", value = fiveUtrGrl, envir = txdbEnv)
  }

  return(fiveUtrGrl)
}

##################################################################################

get_txdb_threeUtr_grl <- function(txdb){
  if(exists("threeUtrGrl", envir=txdbEnv, inherits=FALSE)) {
    threeUtrGrl <- get(x = "threeUtrGrl", envir = txdbEnv)
  } else{
    threeUtrGrl <- GenomicFeatures::threeUTRsByTranscript(txdb)
    assign(x = "threeUtrGrl", value = threeUtrGrl, envir = txdbEnv)
  }

  return(threeUtrGrl)
}


##################################################################################


get_txdb_exons_gr <- function(txdb){
  if(exists("exonsGr", envir=txdbEnv, inherits=FALSE)) {
    exonsGr <- get(x = "exonsGr", envir = txdbEnv)
  } else{
    exonsGr <- unlist(GenomicFeatures::exonsBy(x = txdb, by = "tx"))
    mcols(exonsGr)$tx_id <- names(exonsGr)
    assign(x = "exonsGr", value = exonsGr, envir = txdbEnv)
  }

  return(exonsGr)
}


##################################################################################


get_txdb_introns_gr <- function(txdb){
  if(exists("intronsGr", envir=txdbEnv, inherits=FALSE)) {
    intronsGr <- get(x = "intronsGr", envir = txdbEnv)
  } else{
    intronsGr <- unlist(GenomicFeatures::intronsByTranscript(x = txdb))
    mcols(intronsGr)$tx_id <- names(intronsGr)
    assign(x = "intronsGr", value = intronsGr, envir = txdbEnv)
  }

  return(intronsGr)
}


##################################################################################








