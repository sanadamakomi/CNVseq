#' @title Hg19 seqinfo
#' @export
#' @importFrom GenomeInfoDb Seqinfo
get_hg19_seqinfo <- function() {
    Seqinfo(
    seqnames=c(as.character(as.vector(seq(1, 22, 1))), "X", "Y"),
    seqlengths=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566),
    isCircular=rep(FALSE, 24),
    genome="hg19")
}
