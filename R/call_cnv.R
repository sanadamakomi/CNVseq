#' @title Calling CNV
#' @param data A data frame of log2ratio to do CBS. It must has three column:
#'   chromosome, end, log2.
#' @param s_id A charactor string of sample id.
#' @param is_male A bool value, TRUE when the gender is male.
#'   change-points, default: 0.05.
#' @param fraction Fraction of CNV >=100Kb and >=5Mb, default:c(0.5, 0.5).
#' @export
call_cnv <- function(data, s_id, is_male=TRUE, fraction=c(0.5, 0.5)) {
	hg19_seqinfo <- get_hg19_seqinfo()
    # fraction:>=100kb, >=5M
    # x y ploidy
    if (is_male) {
        x_ploidy <- 1
        y_ploidy <- 1
    } else {
        x_ploidy <- 2
        y_ploidy <- 0
    }

    # <5M: fractio >=50%，cn<=1.5 or cn >=2.5
    result <- do_cbs(data, s_id, alpha=0.05, min.width=2)
    result$ploidy <- 2
    result[which(as.character(as.vector(result$chrom)) == "X"), "ploidy"] <- x_ploidy
    result[which(as.character(as.vector(result$chrom)) == "Y"), "ploidy"] <- y_ploidy
    result$cn <- apply(result[,c("seg.mean", "ploidy")], 1, function(x) {
        round(x[2]*2^x[1],2)
    })
    cn_bottom <- sapply(result$ploidy, function(x){ploidy_to_cn(fraction[1], x, x-1)})
    cn_up <- sapply(result$ploidy, function(x){ploidy_to_cn(fraction[1], x, x+1)})
    result$svtype <- "CNV"
    result[which(result$cn < cn_bottom), "svtype"] <- "DEL"
    result[which(result$cn > cn_up), "svtype"] <- "DUP"

    result$positive <- FALSE
    result[which(result$cn < cn_bottom & (result$loc.end - result$loc.start) >= 1E5), "positive"] <- TRUE
    result[which(result$cn > cn_up & (result$loc.end - result$loc.start) >= 2E5), "positive"] <- TRUE

    result$predict <- 0
    result[which((result$cn <= (result$ploidy - 0.9) | result$cn >= (result$ploidy + 0.9)) & result$positive), "predict"] <- 1

    # >5M: fraction >= 30%，cn>=1.9 or cn <= 2.1
    data_2 <- merge_bin_df(data, group_size=50)
    result_2 <- do_cbs(data_2, s_id, alpha=1E-2, min.width=2)
    result_2$ploidy <- 2
    result_2[which(as.character(as.vector(result_2$chrom)) == "X"), "ploidy"] <- x_ploidy
    result_2[which(as.character(as.vector(result_2$chrom)) == "Y"), "ploidy"] <- y_ploidy
    result_2$cn <- apply(result_2[,c("seg.mean", "ploidy")], 1, function(x) {
        round(x[2]*2^x[1],2)
    })

    cn_bottom_2 <- sapply(result_2$ploidy, function(x){ploidy_to_cn(fraction[2], x, x-1)})
    cn_up_2 <- sapply(result_2$ploidy, function(x){ploidy_to_cn(fraction[2], x, x+1)})

    result_2$svtype <- "CNV"
    result_2[which(result_2$cn < cn_bottom_2), "svtype"] <- "DEL"
    result_2[which(result_2$cn > cn_up_2), "svtype"] <- "DUP"

    result_2$positive <- FALSE
    result_2[which(result_2$cn < cn_bottom_2 & (result_2$loc.end - result_2$loc.start) >= 5E6), "positive"] <- TRUE
    result_2[which(result_2$cn > cn_up_2 & (result_2$loc.end - result_2$loc.start) >= 5E6), "positive"] <- TRUE

    result_2$predict <- 0
    result_2[which((result_2$cn <= (result_2$ploidy - 0.9) | result_2$cn >= (result_2$ploidy + 0.9)) & result_2$positive), "predict"] <- 1

    cnv_5m <- result_2[which(result_2$positive),,drop=FALSE]

    # merge
    result <- merge_segment(result)
    if (nrow(cnv_5m) > 0) {
        # merge and exchange cnv_5m
        cnv_5m_gr <- GRanges(
            Rle(as.character(as.vector(cnv_5m[,"chrom"]))),
            IRanges(start=as.numeric(as.vector(cnv_5m[,"loc.start"])),
                    end=as.numeric(as.vector(cnv_5m[,"loc.end"]))),
            seqinfo=get_hg19_seqinfo)
        result_gr <- GRanges(
            Rle(as.character(as.vector(result[,"chrom"]))),
            IRanges(start=as.numeric(as.vector(result[,"loc.start"])),
                    end=as.numeric(as.vector(result[,"loc.end"]))),
            seqinfo=hg19_seqinfo)
        hit <- findOverlaps(result_gr, cnv_5m_gr)
        if (length(hit) > 0) {
            result <- result[setdiff(1:nrow(result), queryHits(hit)),]
        }
        result <- rbind(result, cnv_5m)
    }

    # sort
    result$chrom <- factor(result$chrom, levels = c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    result <- result[with(result, order(result$chrom, result$loc.start, result$loc.end)), ]

    # if ovelape chromosome > 99%, then output the full length
    for (i in nrow(result)) {
        if (!result[i, "positive"]) next
        seq_len <- seqlengths(hg19_seqinfo[as.character(result[i, "chrom"])])
        if ((result[i, "loc.end"] - result[i, "loc.start"])/seq_len > 0.99) {
            result[i, "loc.start"] <- 10000
            result[i, "loc.end"] <- seq_len
        }
    }
    return(result)
}
