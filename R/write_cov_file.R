#' @title Write read count result file
#' @param gr A GRange with a column named rc.
#' @param path Path to write to.
#' @export
write_cov_file <- function(gr, path) {
    df <- as.data.frame(gr)
    cat(paste(as.character(as.vector(df[, "seqnames"])),
              as.numeric(as.vector(df[, "start"])), as.numeric(as.vector(df[, "end"])), as.numeric(as.vector(df[, "rc"])), sep = '\t'), sep = '\n', file = path)
}
