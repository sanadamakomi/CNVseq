#' @title Create VCF INFO column
#' @param x A data frame of CNV result.
#' @export
make_vcf_info <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { gsub(" ", "", as.character(as.vector(y))) })
    x <- matrix(data=x, ncol=n_col)
    sapply(1:nrow(x), function(i) {
        paste(paste0(toupper(col_n), "=", as.character(as.vector(x[i,,drop=FALSE]))), collapse=";")
    })
}
