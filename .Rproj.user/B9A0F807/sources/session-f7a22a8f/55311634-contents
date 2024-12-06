#' @title Count copy numer with different mosaicism level
#' @param f Fraction of variants.
#' @param ploidy the ploidy of chromosome.
#' @param cn Copy number.
#' @export
ploidy_to_cn <- function(f, ploidy, cn) {
    # Fraction to copy number
    if (cn < 0) {
        return(0)
    }
    if (ploidy == 0) {
        return(100)
    }
    all <- f*cn+(1-f)*ploidy
    lr <- log2(all/ploidy)
    return(round(ploidy*2^lr, 1))
}
