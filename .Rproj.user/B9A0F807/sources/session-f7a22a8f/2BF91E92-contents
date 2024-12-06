#' @title Extract sample id from file path
#' @param file Path of file.
#' @param ptn A charactor for spliting string. By default NULL, it will split
#'   string by '.'.
#' @export
get_id <- function(file, ptn=NULL) {
    sampleid  <- rev(unlist(strsplit(file, "/")))[1]
    if(is.null(ptn)) {
        sampleid  <- unlist(strsplit(sampleid, '\\.'))[1]
    }else{
        sampleid  <- unlist(strsplit(sampleid, ptn))[1]
    }
    return(sampleid)
}
