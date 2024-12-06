#' @title Merge CNN, VCF, CNR, CNS file.
#' @description Input cnn_paths or cnn_dirs, it will merge multiple samples'
#'   read count result files(.cnn) and create a merged file. Input auto_file and
#'   sex_file, it will merge autosome and sex chromosome result, the VCF, CNR,
#'   CNS file can be input.
#' @param cnn_paths Path of CNN files, separated by comma(,).
#' @param cnn_dirs Path of CNN directory, separated by comma(,).
#' @param auto_file path of VCF, CNR, CNS file.
#' @param sex_file path of VCF, CNR, CNS file.
#' @param out_path Output file path.
#' @export
batch_merge <- function(cnn_paths, cnn_dirs, auto_file, sex_file, out_path) {
    if (is.null(cnn_paths) & is.null(cnn_dirs) & (is.null(auto_file) | is.null(sex_file))) stop('Input --cnnFiles or --cnnDirs, or --autoFile and --sexFile')
    if (!is.null(cnn_paths)) {
        cnn_path_lst <- unlist(strsplit(cnn_paths, ','))

        if (is.null(out_path)) out_path <- "merge.cnm"
        if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
        merger_cov_files(cnn_path_lst, out_path)

    } else if (!is.null(cnn_dirs)) {
        cnn_path_lst <- c()
        cnn_dir_lst <- unlist(strsplit(cnn_dirs, ","))
        for (cnn_dir in cnn_dir_lst) {
            cnn_path_lst <- c(cnn_path_lst, list.files(path = cnn_dir, pattern = ".cnn$", all.files = TRUE, full.names = TRUE, recursive = TRUE))
        }

        if (is.null(out_path)) out_path <- "merge.cnm"
        if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
        merger_cov_files(cnn_path_lst, out_path)

    } else {
        if (is.null(out_path)) out_path <- basename(auto_file)
        if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
        merge_result_files(auto_file, sex_file, out_path)
    }
}
