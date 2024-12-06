#' @title Calculate gender and total read from a coverage file.
#' @param cnn_paths Path of CNN files, separated by comma(,).
#' @param cnn_dirs Path of CNN directory, separated by comma(,).
#' @param out_dir Output directory path.
#' @export
#' @importFrom utils read.table
#' @importFrom utils write.table
batch_gender <- function(cnn_paths, cnn_dirs, out_dir) {
    if (is.null(cnn_paths) & is.null(cnn_dirs)) stop('Input --cnnFiles or --cnnDirs ')

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir <- normalizePath(out_dir)

    if (!is.null(cnn_paths)) {
        cnn_path_lst <- unlist(strsplit(cnn_paths, ','))
    } else {
        cnn_path_lst <- c()
        cnn_dir_lst <- unlist(strsplit(cnn_dirs, ","))
        for (cnn_dir in cnn_dir_lst) {
            cnn_path_lst <- c(cnn_path_lst, list.files(path = cnn_dir, pattern = ".cnn$", all.files = TRUE, full.names = TRUE, recursive = TRUE))
        }
    }

    for (cnn_path in cnn_path_lst) {
        df <- read.table(cnn_path, header = FALSE, sep = "\t", quote = "",
                         comment.char = "#", na.strings = "NA",
                         fill = TRUE, stringsAsFactors = FALSE)
        s_id <- get_id(cnn_path)
        out_path <- paste0(out_dir, "/", s_id, ".gender.txt")

        idx <- which(as.character(as.vector(df[,1])) == "Y" &
                     as.numeric(as.vector(df[,2])) <= 2655723 &
                     as.numeric(as.vector(df[,3])) >= 2654895)
        gender <- ifelse(sum(as.numeric(as.vector(df[idx, 4])), na.rm=TRUE) > 0, "Male", "Female")
        total_read <- sum(as.numeric(as.vector(df[,4])), na.rm=TRUE)
        out_df <- data.frame(
            "sample"=s_id,
            "sex"=gender,
            "total.read"=total_read
        )
        write.table(out_df, file = out_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        write(paste0("Write to path: \n", out_path), stdout())
    }
}
