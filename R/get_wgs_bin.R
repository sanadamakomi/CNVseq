#' @title Split wgs region into bins
#' @param bin An integer of bin size, default: 1E4.
#' @param access_bed BED file include sequence-accessible region. If NULL it
#'   will use whole genome region.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @importFrom GenomicRanges trim
#' @importFrom utils read.table
get_wgs_bin <- function(bin=1E4, access_bed="") {
	hg19_seqinfo <- get_hg19_seqinfo()
    if (file.exists(access_bed)) {
        indat <- read.table(access_bed, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
        indat <- indat[which(indat[,1] %in% seqnames(hg19_seqinfo)),,drop=FALSE]
        wgs_gr <- GRanges(Rle(as.character(as.vector(indat[, 1]))),
                          IRanges(start = as.numeric(as.vector(indat[, 2])),
                                  end = as.numeric(as.vector(indat[, 3]))),
                          seqinfo = hg19_seqinfo)
        wgs_gr <- trim(wgs_gr)
    } else {
        wgs_gr <- GRanges(Rle(seqnames(hg19_seqinfo)),
                          IRanges(start = 1,
                                  end = seqlengths(hg19_seqinfo)),
                          seqinfo = hg19_seqinfo)
    }

    # split genomic range into bins
    wgs_gr_lst <- GRangesList()
    for (i in 1:length(wgs_gr)) {
        wgs_gr_lst[[i]] <- creat_chunk(wgs_gr[i], bin=bin)
    }
    return(unlist(wgs_gr_lst))
}
