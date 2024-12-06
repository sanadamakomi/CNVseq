#' @title Calculate read count in bams
#' @param bam_dirs Path of BAM directory, separated by comma(,).
#' @param bam_paths Path of BAM file, separated by comma(,).
#' @param bed_path Path of BED file.
#' @param thread An integer providing the number of thread, default: 4.
#' @param mapq.filter A non-negative integer specifying the minimum mapping
#'   quality to include. BAM reads with mapping qualities less than mapqFilter
#'   are discarded, default: 0.
#' @param out_dir A character string of directory to output coverage files.
#' @export
batch_cov <- function(bam_dirs, bam_paths, bed_path, out_dir, thread=4, mapq.filter=0) {
    if (is.null(bam_paths) & is.null(bam_dirs)) stop('Input --bamFiles or --bamDirs ')
    if (is.null(bed_path)) stop('Input --bedFile')

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir <- normalizePath(out_dir)

    if (!is.null(bam_paths)) {
        bam_path_lst <- unlist(strsplit(bam_paths, ','))
    } else {
        bam_path_lst <- c()
        bam_dir_lst <- unlist(strsplit(bam_dirs, ","))
        for (bam_dir in bam_dir_lst) {
            bam_path_lst <- c(bam_path_lst, list.files(path = bam_dir, pattern = ".bam$", all.files = TRUE, full.names = TRUE, recursive = TRUE))
        }
    }
    check_bam(bam_path_lst)

    bin_gr <- bed_to_gr(bed_path)

    for (bam_path in bam_path_lst) {
        count_reads_of_region_multicore(bin_gr, bam_path, thread=thread, batch=1000, tmpDir=out_dir, mapq.filter=mapq.filter, minoverlap=75L)
    }
}
