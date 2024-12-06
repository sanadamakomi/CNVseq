#' @title Read annovation database file and return granges
#' @param prefix A prefix of database file, e.g. centromere_telomere for
#'   hg19_centromere_telomere.txt.
#' @param db_path Path of database directory.
#' @param p A integer for the column index of chromosome, default: 1.
#' @param header A bool for the header of database file, default: FALSE.
#' @param file_encode The file endode of database file, default: UTF-8.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomicAlignments
#' @importFrom utils read.table
#' @importFrom GenomicRanges trim
read_db_file <- function(prefix, db_path, p=1, header=FALSE, file_encode='UTF-8') {
    hg19_seqinfo <- get_hg19_seqinfo()
	file <- paste0(db_path, "/hg19_", prefix, ".txt")
    db_df <- read.table(file, header = header, sep = "\t", quote="", comment.char="#", na.strings= "NA",fill=TRUE, stringsAsFactors = FALSE, fileEncoding=file_encode)
    db_df <- as.data.frame(t(apply(db_df, 1, function(x){gsub("\\n", "", x)})))
    idx <- which(as.character(as.vector(db_df[,p])) %in% seqnames(hg19_seqinfo))
    db_df <- db_df[idx,,drop=FALSE]
    gr <- GRanges(Rle(as.character(as.vector(db_df[,p]))),
                  IRanges(start=as.numeric(as.vector(db_df[,p+1]))+1,
                          end=as.numeric(as.vector(db_df[,p+2]))),
                  seqinfo=hg19_seqinfo)
    gr <- trim(gr)
    mcols(gr) <- db_df[, (p+3):ncol(db_df), drop=FALSE]
    return(gr)
}
