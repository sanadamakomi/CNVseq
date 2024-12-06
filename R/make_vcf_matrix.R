#' @title Create VCF matrix
#' @param x A data frame of CNV result.
#' @export
make_vcf_matrix <- function(x) {
    CHROM <- as.character(as.vector(x[, "chrom"]))
    POS <- as.numeric(as.vector(x[, "loc.start"]))
    ID <- rep(".", nrow(x))
    REF <- rep("N", nrow(x))
    ALT <- paste0("<", as.character(as.vector(x[, "svtype"])), ">")

    x$gq <- as.character(as.vector(round(x[,"quality"],0)))
    x[which(is.na(x[,"quality"])), "gq"] <- "0"
    x[which(is.infinite(x[,"quality"])), "gq"] <- "0"

    if ("gq" %in% names(x)) {
        QUAL <- as.character(as.vector(x[,"gq"]))
    } else {
        QUAL <- rep(".", nrow(x))
    }
    if ("filter" %in% names(x)) {
        FILTER <- as.character(as.vector(x[, "filter"]))
    } else {
        FILTER <- rep(".", nrow(x))
    }
    x$svlen <- as.numeric(as.vector(x[, "loc.end"])) - as.numeric(as.vector(x[, "loc.start"])) + 1
    x$fc <- 2^as.numeric(as.vector(x[, "seg.mean"]))
    info.vec <- c("loc.end", "svtype", "svlen", "fc", "seg.mean", "num.mark", "predict")
    new.name <- c("END", "SVTYPE", "SVLEN", "FOLD_CHANGE", "FOLD_CHANGE_LOG", "PROBES", "PREDICT")
    info.df <- x[, info.vec, drop = FALSE]
    colnames(info.df) <- new.name
    INFO <- make_vcf_info(info.df)

    x$gt <- "0/1"
    x[which(x$svtype == "DEL" & x$gt < (x$ploidy/2)), "gt"] <- "1/1"
    x[which(x$svtype == "DUP" & x$gt > (x$ploidy+1)), "gt"] <- "1/1"

    FORMAT <- make_vcf_format(x[, c("gt", "gq", "cn"), drop=FALSE])
    cbind(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, rep("GT:GQ:CN", nrow(x)), FORMAT)
}
