#' @title Split wgs region into bins
#' @param path Path of output, a BED file.
#' @param bin An integer of bin size, default: 1E4.
#' @param access_bed BED file include sequence-accessible region. If NULL it
#'   will use whole genome region.
#' @export
batch_split <- function(path, bin=1E4, access_bed="") {
    wgs_bin <- get_wgs_bin(as.integer(bin), access_bed)
    df <- as.data.frame(wgs_bin)
    if (is.null(path)) path <- "target.bed"
    cat(paste(as.character(as.vector(df[, "seqnames"])),
              as.numeric(as.vector(df[, "start"])), as.numeric(as.vector(df[, "end"])), sep = '\t'), sep = '\n', file = path)
}
