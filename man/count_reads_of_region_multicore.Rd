% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_reads_of_region_multicore.R
\name{count_reads_of_region_multicore}
\alias{count_reads_of_region_multicore}
\title{Calculating read count in multiple thread}
\usage{
count_reads_of_region_multicore(
  region,
  bamPath,
  thread,
  batch,
  tmpDir = NULL,
  mapq.filter = 0,
  minoverlap = 75L
)
}
\arguments{
\item{region}{A grange object of region to extract reads in BAM file.}

\item{bamPath}{A character string of the BAM path.}

\item{thread}{An integer providing the number of thread.}

\item{batch}{An integer giving how many GRanges are performed in a batch.}

\item{tmpDir}{A character string of directory to output coverage files
(<sampleid>.cnn). Default is the current folder.}

\item{mapq.filter}{A non-negative integer specifying the minimum mapping
quality to include. BAM reads with mapping qualities less than mapqFilter
are discarded.}

\item{minoverlap}{Minimum overlap size for region and reads, default: 75L.}
}
\description{
Calculating read count in multiple thread
}
