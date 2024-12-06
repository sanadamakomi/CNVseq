#' @title Create grange into chunk
#' @param gr A grange object.
#' @param bin Bin size, default: 1E4.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @importFrom GenomicRanges trim
creat_chunk <- function(gr, bin=1E4){
    ck <- breakInChunks(width(gr), chunksize=bin)
    chunks <- GRanges(Rle(rep(seqnames(gr), length(ck))), IRanges(ck), seqinfo=seqinfo(gr))
    chunks <- GenomicRanges::shift(chunks, start(gr)-1)
    return(trim(chunks))
}
