#' @title Filtering transcripts
#' @param x A data.frame of CNV calling result.
#' @export
merge_segment <- function(x) {
    out_dat <- list()
    for (chr in unique(x[,"chrom"])) {
        idx <- which(x[,"chrom"] == chr)
        dat <- x[idx,,drop=FALSE]
        if (nrow(dat) > 2) {
            group_indices <- cumsum(c(TRUE, diff(as.numeric(factor(dat$svtype))) != 0))
            segments <- split(dat, group_indices)
            segment_counts <- sapply(segments, function(x) nrow(x))
            for (n in names(segments)) {
                dat_seg <- segments[[n]]
                if (nrow(dat_seg) > 1 & dat_seg[1,"svtype"] != "CNV") {
                    sv_size <- dat_seg[nrow(dat_seg), "loc.end"] - dat_seg[1, "loc.start"] + 1
                    mean_cn <- mean(dat_seg[, "cn"], na.rm=TRUE)
                    positive_flag <- FALSE
                    predict_flag <- 0
                    if (dat_seg[1,"svtype"] == "DEL" & sv_size >= 1E5) {
                        positive_flag <- TRUE
                        if (mean_cn <= (dat_seg[1, "ploidy"] - 0.9)) {
                            predict_flag <- 1
                        }

                    }
                    if (dat_seg[1,"svtype"] == "DUP" & sv_size >= 2E5) {
                        positive_flag <- TRUE
                        if (mean_cn >= (dat_seg[1, "ploidy"] + 0.9)) {
                            predict_flag <- 1
                        }
                    }

                    dat_seg_new <- data.frame(
                        "ID"=dat_seg[1, "ID"],
                        "chrom"=dat_seg[1, "chrom"],
                        "loc.start"=dat_seg[1, "loc.start"],
                        "loc.end"=dat_seg[nrow(dat_seg), "loc.end"],
                        "num.mark"=sum(dat_seg[, "num.mark"]),
                        "seg.mean"=mean(dat_seg[, "seg.mean"], na.rm=TRUE),
                        "seg.sd"=mean(dat_seg[, "seg.sd"], na.rm=TRUE),
                        "seg.median"=mean(dat_seg[, "seg.median"], na.rm=TRUE),
                        "seg.mad"=mean(dat_seg[, "seg.mad"], na.rm=TRUE),
                        "quality"=mean(dat_seg[, "quality"], na.rm=TRUE),
                        "cn"=mean_cn,
                        "ploidy"=dat_seg[1, "ploidy"],
                        "svtype"=dat_seg[1, "svtype"],
                        "positive"=positive_flag,
                        "predict"=predict_flag
                    )
                    segments[[n]] <- dat_seg_new
                }
            }
            out_dat[[chr]] <- do.call("rbind", segments)
        } else {
            out_dat[[chr]] <- dat
        }
    }
    out_df <- do.call("rbind", out_dat)
    rownames(out_df) <- NULL
    return(out_df)
}
