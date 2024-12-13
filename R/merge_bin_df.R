#' @title Merge log2ratio with a specific interval of bins
#' @param data A data frame of cnr file.
#' @param group_size A interval to merge, default: 50.
#' @export
merge_bin_df <- function(data, group_size=50) {
    data$group <- with(data, paste(data$chromosome, ceiling(seq_along(data$end) / group_size), sep = "_"))
    safe_median <- function(x) median(x[is.finite(x)], na.rm = TRUE)
    result <- aggregate(
        log2 ~ group + chromosome,
        data = data,
        FUN = safe_median
    )
    result$end <- sapply(result$group, function(g) {
        group_chr <- unlist(strsplit(g, "_"))[1]
        group_number <- as.numeric(unlist(strsplit(g, "_"))[2])
        end_position <- min(group_number * group_size, max(data$end[data$chromosome == group_chr]))
        if (end_position > nrow(data)) end_position <- nrow(data)
        data$end[end_position]
    })
    return(result[,c("chromosome", "end", "log2"),drop=FALSE])
}
