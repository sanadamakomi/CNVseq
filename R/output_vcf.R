#' @title Output CNV calling result to VCF file
#' @param dat A data.frame of CNV calling result.
#' @param path Path of VCF file.
#' @export
output_vcf <- function(dat, path) {
    id <- get_id(path)
    header <- make_vcf_header()
    vcf_matrix <- make_vcf_matrix(dat)
    cat(header, file=path, sep="\n")
    cat(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", id), collapse = "\t"),
        file = path, sep = "\n", append = TRUE)
    cat(apply(vcf_matrix, 1, paste, collapse = "\t"), file = path, sep = "\n", append = TRUE)

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}
