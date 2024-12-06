#' @title Create VCF header
#' @export
make_vcf_header <- function() {
    fileformat <- "##fileformat=VCFv4.2"
    fileDate <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
    alt <- c('##ALT=<ID=DEL,Description="Deletion">',
             '##ALT=<ID=DUP,Description="Duplication">')
    info <- c('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
              '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
              '##INFO=<ID=FOLD_CHANGE,Number=1,Type=Float,Description="Fold change">',
              '##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description="Log fold change">',
              '##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of probes in CNV">',
              '##INFO=<ID=PREDICT,Number=1,Type=Integer,Description="Predict of True CNV">'
    )
    format <- c('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype.">',
                '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality.">',
                '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events.">')
    c(fileformat, fileDate, alt, info, format)
}
