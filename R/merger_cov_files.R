#' @title Merge coverage files.
#' @description Align and merge coverage files (<filename>.cnn) with chromosome,
#'   start and end position. Four required fields in a coverage file are
#'   chromosome name, start and end position, depth of coverage.
#' @param files A character vector contains several coverage files path.
#' @param path Path to write to.
#' @return A data frame, of which columns are chromosome, start position, end
#'   position, and depths in input coverage files.
#' @export
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom data.table data.table
merger_cov_files <- function(files, path = NULL) {
    if (is.null(path)) path <- 'mergedCov.txt'

    write('Start to merge coverage files...', stdout())

    coverageMatrix <- NULL
    for (i in files) {
        id <- get_id(i)
        indat <- read.table(i, header = FALSE, sep = "\t", quote = "", comment.char = "#",
                            na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
        colnames(indat) <- c("chr", "start", "end", id)
        indat2 <- data.table(indat, key = c("chr", "start", "end"))
        if (is.null(coverageMatrix)) {
            coverageMatrix <- indat
        } else {
            coverageMatrix <- merge(coverageMatrix, indat[,c("chr", "start", "end", id)], by = c("chr", "start", "end"), all=TRUE)
        }
    }
    coverageMatrix$chr <- factor(coverageMatrix$chr,
                                 levels = c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    coverageMatrix <- coverageMatrix[, c('chr', 'start', 'end',
                                         setdiff(colnames(coverageMatrix), c('chr', 'start', 'end')))]
    write.table(coverageMatrix[with(coverageMatrix,
                                    order(coverageMatrix$chr, coverageMatrix$start, coverageMatrix$end)), ],
                file = path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}
