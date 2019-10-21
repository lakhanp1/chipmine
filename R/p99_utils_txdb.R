


get_txdb_transcripts_gr <- function(txdb, excludeType = NULL, tx = NULL){

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

    ## tx filter list
    txFilter <- NULL
    txFilter$tx_type <- selectType
    if(!is.null(tx)){
      txFilter$tx_id <- tx
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


get_txdb_fiveUtr_gr <- function(txdb, tx = NULL){

  if(exists("fiveUtrGr", envir=txdbEnv, inherits=FALSE)) {
    fiveUtrGr <- get(x = "fiveUtrGr", envir = txdbEnv)

  } else{
    fiveUtrGrl <- GenomicFeatures::fiveUTRsByTranscript(txdb)

    if(!is.null(tx)){
      fiveUtrGrl <- fiveUtrGrl[names(fiveUtrGrl) %in% tx]
    }

    ## remember to combine the multi-exon UTRs from UTR GRangesList
    fiveUtrGr <- unlist(range(fiveUtrGrl))
    mcols(fiveUtrGr)$tx_id <- names(fiveUtrGr)
    assign(x = "fiveUtrGr", value = fiveUtrGr, envir = txdbEnv)
  }

  return(fiveUtrGr)
}

##################################################################################

get_txdb_threeUtr_gr <- function(txdb, tx = NULL){
  if(exists("threeUtrGr", envir=txdbEnv, inherits=FALSE)) {
    threeUtrGr <- get(x = "threeUtrGr", envir = txdbEnv)

  } else{
    threeUtrGrl <- GenomicFeatures::threeUTRsByTranscript(txdb)

    if(!is.null(tx)){
      threeUtrGrl <- threeUtrGrl[names(threeUtrGrl) %in% tx]
    }

    ## remember to combine the multi-exon UTRs from UTR GRangesList
    threeUtrGr <- unlist(range(threeUtrGrl))
    mcols(threeUtrGr)$tx_id <- names(threeUtrGr)
    assign(x = "threeUtrGr", value = threeUtrGr, envir = txdbEnv)
  }

  return(threeUtrGr)
}


##################################################################################


get_txdb_exons_gr <- function(txdb, tx = NULL){
  if(exists("exonsGr", envir=txdbEnv, inherits=FALSE)) {
    exonsGr <- get(x = "exonsGr", envir = txdbEnv)
  } else{
    exonsGrl <- GenomicFeatures::exonsBy(x = txdb, by = "tx")

    if(!is.null(tx)){
      exonsGrl <- exonsGrl[names(exonsGrl) %in% tx]
    }

    exonsGr <- unlist(exonsGrl)
    mcols(exonsGr)$tx_id <- names(exonsGr)
    assign(x = "exonsGr", value = exonsGr, envir = txdbEnv)
  }

  return(exonsGr)
}


##################################################################################


get_txdb_introns_gr <- function(txdb, tx = NULL){
  if(exists("intronsGr", envir=txdbEnv, inherits=FALSE)) {
    intronsGr <- get(x = "intronsGr", envir = txdbEnv)
  } else{
    intronsGrl <- GenomicFeatures::intronsByTranscript(x = txdb)

    if(!is.null(tx)){
      intronsGrl <- intronsGrl[names(intronsGrl) %in% tx]
    }

    intronsGr <- unlist(intronsGrl)
    mcols(intronsGr)$tx_id <- names(intronsGr)
    assign(x = "intronsGr", value = intronsGr, envir = txdbEnv)
  }

  return(intronsGr)
}


##################################################################################


