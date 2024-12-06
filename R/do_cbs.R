#' @title Circular Binary Segmentation
#' @param data A data frame of log2ratio to do CBS. It must has three column:
#'   chromosome, end, log2.
#' @param s_id A charactor string of sample id.
#' @param alpha A numeric value of significance levels for the test to accept
#'   change-points, default: 0.05.
#' @param min.width An integer value of the minimum number of markers for a
#'   changed segment, default: 2.
#' @export
#' @import DNAcopy
do_cbs <- function(data, s_id, alpha=0.05, min.width=2) {
    cna_project <- CNA(cbind(data$log2),
                       data$chromosome,
                       data$end,
                       data.type = "logratio",
                       sampleid = s_id)
    # Smooth
    cna_project_smooth <- smooth.CNA(cna_project)

    # Segment
    cna_project_segment <- segment(cna_project_smooth, alpha=alpha, min.width=min.width, verbose=1)
    result <- segments.summary(cna_project_segment)
    # Output
    p_result <- segments.p(cna_project_segment)
    result$quality <- round(-log10(p_result$pval),2)
    result$cn <- round(2*2^result$seg.mean, 2)
    return(result)
}
