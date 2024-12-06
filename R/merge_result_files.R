#' @title Merge autosome and sex chromosome result
#' @param auto_file path of VCF, CNR, CNS file.
#' @param sex_file path of VCF, CNR, CNS file.
#' @param out_path path of output file.
#' @export
merge_result_files <- function(auto_file, sex_file, out_path) {
    if (is.null(auto_file) | is.null(sex_file)) stop("Input --autoFile and --sexFile")
    if (!file.exists(auto_file)) stop(paste0("Error: File not exists:\n", auto_file))
    # write
    in_con <- file(auto_file, "r")
    out_con <- file(out_path, "w")

    # read and write line by line
    while (length(line <- readLines(in_con, n = 1, warn = FALSE)) > 0) {
        flag <- grepl("^X|^Y|^chrX|^chrY", line)
        if (!flag) {
            writeLines(line, out_con)
        }
    }
    close(in_con)

    if (file.exists(sex_file)) {
        in_con <- file(sex_file, "r")
        while (length(line <- readLines(in_con, n = 1, warn = FALSE)) > 0) {
            flag <- grepl("^X|^Y|^chrX|^chrY", line)
            if (flag) {
                writeLines(line, out_con)
            }
        }
        close(in_con)
    }

    close(out_con)
    write(paste0("Write to path: \n", normalizePath(out_path)), stdout())
}
