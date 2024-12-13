#' @title Call CNVs from read count files
#' @param cnm_path Path of merged CNN file.
#' @param cnn_paths Path of CNN files, separated by comma(,).
#' @param cnn_dirs Path of CNN directory, separated by comma(,).
#' @param out_dir Output directory path.
#' @param test_id Sample id to call CNV.
#' @param ref_ids Reference ids to create a baseline, separated by comma(,).
#' @param by_gender A bool value to call by gender.
#' @param fraction Fraction of CNV >=100Kb and >=5Mb, default:c(0.5, 0.5).
#' @export
#' @importFrom utils read.table
#' @importFrom utils write.table
batch_call_cnv <- function(cnm_path, cnn_paths, cnn_dirs, out_dir, test_id, ref_ids, by_gender=FALSE, fraction=c(0.5, 0.5)) {
    if (is.null(test_id)) stop('Input --testId')
    if (is.null(cnm_path) & is.null(cnn_paths) & is.null(cnn_dirs)) stop('Input --cnmFile or --cnnFiles or --cnnDirs')

    # create output dir
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir <- normalizePath(out_dir)

    # create merge file
    if (is.null(cnm_path)) {
        if (!is.null(cnn_paths)) {
            cnn_path_lst <- unlist(strsplit(cnn_paths, ','))
        } else {
            cnn_path_lst <- c()
            cnn_dir_lst <- unlist(strsplit(cnn_dirs, ","))
            for (cnn_dir in cnn_dir_lst) {
                cnn_path_lst <- c(cnn_path_lst, list.files(path = cnn_dir, pattern = ".cnn$", all.files = TRUE, full.names = TRUE, recursive = TRUE))
            }
        }
        cnm_path <- paste0(out_dir, "/merge.cnm")
        merger_cov_files(cnn_path_lst, cnm_path)
    }
    # call cnv
    cov_df <- read.table(cnm_path, header = TRUE, sep = "\t", quote = "",
                         comment.char = "#", na.strings = "NA",
                         fill = TRUE, stringsAsFactors = FALSE)
    all_id <- setdiff(colnames(cov_df), c("chr", "start", "end"))
    if (!test_id %in% all_id) stop(paste0("Sample [ ", test_id," ] not in file:", cnm_path))

    # call gender
    gender <- apply(cov_df[which(cov_df[,"chr"] == "Y" & cov_df[,"start"] < 2655723 & cov_df[,"end"] > 2654895 ), all_id, drop=FALSE], 2, function(x) {ifelse(sum(x) > 10, "Male", "Female")})
    names(gender) <- all_id

    # reference
    if (!is.null(ref_ids)) {
        ref_id_lst <- unlist(strsplit(ref_ids, ","))
        ref_id_lst <- setdiff(intersect(ref_id_lst, all_id), test_id)
    } else {
        ref_id_lst <- setdiff(all_id, test_id)
    }

    if (by_gender) {
        idx <- which(gender[ref_id_lst] == gender[test_id])
        if (length(idx) == 0) {
            write('No sample has the same gender as the test.', stdout())
            q()
        }
        ref_id_lst <- ref_id_lst[idx]
    }

    # rm other id
    all_id <- c(test_id, ref_id_lst)
    cov_df <- cov_df[,c("chr", "start", "end", all_id), drop=FALSE]

    # normalize
    total_read <- apply(cov_df[, all_id, drop=FALSE], 2, "sum")
    ratio <- min(total_read) / total_read
    names(ratio) <- all_id
    for(i in all_id) {
        cov_df[,i] <- as.integer(cov_df[,i] * ratio[i])
    }

    ref_bin_cov <- apply(cov_df[,ref_id_lst,drop=FALSE], 1, median, na.rm=TRUE)
    test_bin_cov <- cov_df[,test_id,drop=FALSE]

    # mask low coverage
    bin_ave_depth <- apply(cov_df[,ref_id_lst,drop=FALSE], 1, mean, na.rm=TRUE)
    total_ave_depth <- median(bin_ave_depth)
    # low_depth_idx <- which(bin_ave_depth < total_ave_depth*0.1)
    low_depth_idx <- which(bin_ave_depth < 10)

    # calculate log2 ratio
    test_bin_log_ratio <- apply(apply(test_bin_cov, 2, "/", ref_bin_cov), 2, "log2")
    data <- data.frame(
        "chromosome"=as.character(as.vector(cov_df[,"chr"])),
        "start"=as.numeric(as.vector(cov_df[,"start"])),
        "end"=as.numeric(as.vector(cov_df[,"end"])),
        "gene"="gene",
        "log2"=as.numeric(as.vector(test_bin_log_ratio[,test_id])),
        "depth"=test_bin_cov[,test_id],
        "weight"=1
    )

    # mask low depth
    if (length(low_depth_idx) > 0) {
        data[low_depth_idx, "weight"] <- 0
    }
    data <- data[which(data[, "weight"] == 1),,drop=FALSE]

    # save
    cnr_path <- paste0(out_dir, "/", test_id, ".cnr")
    write.table(data, file = cnr_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")

    # call cnv
    is_male <- ifelse(gender[test_id] == "Male", TRUE, FALSE)
    result_all <- call_cnv(data, test_id, is_male, fraction)

    # output all segment
    cns_path <- paste0(out_dir, "/", test_id, ".cns")
    write.table(result_all, file = cns_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")

    # output VCF
    result_positive <- result_all[which(result_all$positive),]
    vcf_path <- paste0(out_dir, "/", test_id, ".cnv.vcf")
    output_vcf(result_positive, vcf_path)
}
