#' @title Calculating read count in multiple thread
#' @param region A grange object of region to extract reads in BAM file.
#' @param bamPath A character string of the BAM path.
#' @param thread An integer providing the number of thread.
#' @param batch An integer giving how many GRanges are performed in a batch.
#' @param tmpDir A character string of directory to output coverage files
#'   (<sampleid>.cnn). Default is the current folder.
#' @param mapq.filter A non-negative integer specifying the minimum mapping
#'   quality to include. BAM reads with mapping qualities less than mapqFilter
#'   are discarded.
#' @param minoverlap Minimum overlap size for region and reads, default: 75L.
#' @export
#' @import parallel
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomicAlignments
#' @importFrom GenomicAlignments start
#' @importFrom GenomicAlignments end
#' @importFrom GenomicAlignments width
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools ScanBamParam
count_reads_of_region_multicore <- function(region, bamPath, thread, batch, tmpDir = NULL, mapq.filter = 0, minoverlap=75L) {
    timestamp <- paste0('depth_', as.character(as.integer(Sys.time())))
    if (is.null(tmpDir)) tmpDir <- '.'
    id <- get_id(bamPath)
    tmp <- normalizePath(tmpDir)
    tmp <- paste0(tmp, '/tmp_', id, '_', timestamp)

    if (! file.exists(tmp) ) dir.create(path = tmp, recursive = TRUE)

	batch.index <- IRanges(breakInChunks(totalsize = length(region),  chunksize = batch))
    chunks <- GRangesList()
    for (i in 1:length(batch.index)) {
        chunks[[i]] <- region[start(batch.index)[i]:end(batch.index)[i]]
    }
    if (thread <= 1){
        a <- lapply(chunks, function(gr){
            gr <- count_reads_of_region(gr, bamPath, mapq.filter = mapq.filter, minoverlap = minoverlap)
            save(gr, file = paste0(tmp, "/", id, '_', timestamp, "_",
                                   seqnames(gr)[1], "_", min(start(gr)), "_", max(end(gr)),".Rdata"))
        })
    } else{
        clsp <- makeCluster(thread, type = "FORK", outfile = paste0(tmp, "/log.txt"))
        a <- parLapply(clsp, chunks, function(gr, ...) {
            gr <-  count_reads_of_region(gr, bamPath, mapq.filter = mapq.filter, minoverlap = minoverlap)
            save(gr, file = paste0(tmp, "/", id, '_', timestamp, "_",
                                   seqnames(gr)[1], "_", min(start(gr)), "_", max(end(gr)),".Rdata"))
        })
        stopCluster(clsp)
    }
    allrdata <- list.files(path = tmp, pattern = ".Rdata", all.files = TRUE,
                           full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
    result <- GRangesList()
    for (i in 1:length(allrdata)) {
        gr <- get(load(allrdata[i]))
        result[[i]] <- gr
    }
    unlink(tmp, recursive = TRUE)
    result <- unlist(result)
    write_cov_file(result, path = paste0(tmpDir, "/",  id, ".cnn"))
    write(paste0("Write to path: \n", tmpDir, "/",  id, ".cnn"), stdout())
}
