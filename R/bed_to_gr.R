#' @title Read BED file and output grange
#' @param file Path of BED file.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomicAlignments
#' @importFrom utils read.table
bed_to_gr <- function(file) {
	hg19_seqinfo <- get_hg19_seqinfo()
    indat <- read.table(file, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
                        fill = TRUE, stringsAsFactors = FALSE)
    if (nrow(indat) < 3) stop(paste0('Error format: ', file))
    chr <- as.character(as.vector(indat[,1]))
    chr <- gsub('chr', '', chr)
    idx <- which(chr %in% c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    if (nrow(indat) == 3) {
        gr <- GRanges(Rle(chr[idx]), IRanges(start = as.numeric(as.vector(indat[idx, 2])),
                                             end = as.numeric(as.vector(indat[idx, 3]))),
                      seqinfo=hg19_seqinfo)
    } else {
        gr <- GRanges(Rle(chr[idx]), IRanges(start = as.numeric(as.vector(indat[idx, 2])),
                                             end = as.numeric(as.vector(indat[idx, 3]))),
                      seqinfo=hg19_seqinfo, id = indat[idx, 4])
    }
}
