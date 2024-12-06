#' @title Check BAM file.
#' @description It will stop if BAM file is illegal.
#' @param x A character string or vector of BAM File path.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom Rsamtools BamFile
check_bam <- function(x) {
    a <- sapply(x, function(bam) {
        if (! file.exists(bam)) stop(paste0(bam, ' file is missing.'))
        bai <- BamFile(bam)$index
        if (is.na(bai)) stop(paste0(bam, '.bai index file is missing.'))
        bam_info <- file.info(bam)
        bai_info <- file.info(bai)
        dt <- difftime(bai_info$mtime, bam_info$mtime, units = 'secs')
        dt <- as.numeric(dt)
        if (dt < 0) stop(paste0(bam, ' index file is older than BAM file.'))
    })
}
