#' @title Calculating read count in region
#' @param region A grange object of region to extract reads in BAM file.
#' @param bamPath A character string of the BAM path.
#' @param mapq.filter A non-negative integer specifying the minimum mapping
#'   quality to include. BAM reads with mapping qualities less than mapqFilter
#'   are discarded.
#' @param minoverlap Minimum overlap size for region and reads, default: 75L.
#' @export
#' @import bitops
#' @import Biostrings
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomicAlignments
#' @importFrom GenomicAlignments start
#' @importFrom GenomicAlignments end
#' @importFrom GenomicAlignments width
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools ScanBamParam
count_reads_of_region <- function(region, bamPath, mapq.filter = 0, minoverlap=75L) {
    what <- c("mapq", "flag")
    flag <- scanBamFlag(isPaired=TRUE, isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE, isSupplementaryAlignment=FALSE, isDuplicate=FALSE)
    param <- ScanBamParam(which = region, what = what, flag = flag, mapqFilter = mapq.filter)
    inBam <- readGAlignments(file = bamPath, index = bamPath, param = param)
    hit <- findOverlaps(region, inBam, minoverlap=minoverlap)
    rc <- sapply(1:length(region), function(x){
        if (x %in% queryHits(hit)) {
            return(length(which(queryHits(hit)==x)))
        }
        return(0)
    })
    mcols(region)$rc <- rc
    print(region)
    return(region)
}
