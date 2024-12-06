#!/usr/bin/Rscript
#' version 1.0, 2024.12.05 compiled by R 4.4.0*
#' CNVseq


#' Blind warning
options(warn = -1)
#' Blind scientific notation
options(scipen=999)

#' Default
VERSION="1.0"
mode="help"
modeHelp="help"
inPath=NULL
outPath=NULL
binSize=1E4
accessBed=""
bedFile=NULL
bamFiles=NULL
bamDirs=NULL
thread=4
mapq=0
outDir="."
cnnFiles=NULL
cnnDirs=NULL
testId=NULL
refIds=NULL
cnmFile=NULL
autoFile=NULL
sexFile=NULL
vcfFile=NULL
smallCNV=0.5
bigCNV=0.5

#' Options
library(GetoptLong)
GetoptLong(
    "modeHelp=s", "Please input --modeHelp name to check options in each mode",
    "mode=s", "'split', 'cov', 'gender', 'merge', 'call'",
    "inPath=s", "Path of input",
    "outPath=s", "Path of output",
    "binSize=i", "Size of bin to calculate coverage",
    "accessBed=s", "BED file include sequence-accessible region",
    "bedFile=s", "Path of BED file",
    "bamFiles=s", "Path of BAM files, separated by comma(,)",
    "bamDirs=s", "Path of BAM directory, separated by comma(,)",
    "thread=i", "Number of thread, default:4",
    "mapq=i", "Mapq to filter reads when calculate coverage, default:0",
    "outDir=s", "Output directory path",
    "cnnFiles=s", "Path of CNN files, separated by comma(,)",
    "cnnDirs=s", "Path of CNN directory, separated by comma(,)",
    "cnmFile=s", "Path of merged CNN file",
    "autoFile=s", "Path of VCF, CNR, CNS file",
    "sexFile=s", "Path of VCF, CNR, CNS file",
    "vcfFile=s", "Path of VCF file",
    "testId=s", "Sample ID",
    "refIds=s", "Reference sample IDs, separated by comma(,)",
    "byGender", "Call CNV by gender",
    "smallCNV=f", "Fraction of CNV >=100Kb, default:0.5",
    "bigCNV=f", "Fraction of CNV >=5Mb, default:0.5",
    head = 'An example to show how to use the script',
    foot = 'Please contact chenzn for comments'
)

#' Mode help
if (modeHelp == "split") {
    write("Usage: split wgs region into bins and output a BED file.", stdout())
    write("Optional:", stdout())
    write("    --outPath\tOut path of BED file, default:target.bed", stdout())
    write("    --binSize\tBin size, default:1E4", stdout())
    write("    --accessBed\tBED file include sequence-accessible region. If NULL it will use whole genome region. Recommend input /tddata/reference/CNV_baseline/cnvkit_region/access.human_g1k_v37_decoy.bed", stdout())
    q()
} else if (modeHelp == "cov") {
    write("Usage: Calculate read count in BED region", stdout())
    write("    --bamFiles\tPath of BAM files, separated by comma(,)", stdout())
    write("    --bamDirs\tPath of BAM directory, separated by comma(,)", stdout())
    write("    --bedFile\tPath of BED file", stdout())
    write("Optional:", stdout())
    write("    --outDir\tOutput directory path, default:.", stdout())
    write("    --thread\tNumber of thread, default:4", stdout())
    write("    --mapq\tMapq to filter reads when calculate coverage, default:0", stdout())
    q()
} else if (modeHelp == "gender") {
    write("Usage: Calculate read count in BED region", stdout())
    write("    --cnnFiles\tPath of CNN files, separated by comma(,)", stdout())
    write("    --cnnDirs\tPath of CNN directory, separated by comma(,)", stdout())
    write("Optional:", stdout())
    write("    --outDir\tOutput directory path, default:.", stdout())
    q()
} else if (modeHelp == "merge") {
    write("Usage: Merge CNN files", stdout())
    write("    --cnnFiles\tPath of CNN files, separated by comma(,)", stdout())
    write("    --cnnDirs\tPath of CNN directory, separated by comma(,)", stdout())
    write("    --autoFile\tPath of VCF, CNR, CNS file", stdout())
    write("    --sexFile\tPath of VCF, CNR, CNS file", stdout())
    write("Optional:", stdout())
    write("    --outPath\tOut path of CNN file, default:merge.cnm", stdout())
    q()
} else if (modeHelp == "call") {
    write("Usage: Call CNV", stdout())
    write("    --testId\tSample ID", stdout())
    write("    --cnmFile\tMerged CNN file:merge.cnm", stdout())
    write("    --cnnFiles\tPath of CNN files, separated by comma(,)", stdout())
    write("    --cnnDirs\tPath of CNN directory, separated by comma(,)", stdout())
    write("Optional:", stdout())
    write("    --outDir\tOutput directory path, default:.", stdout())
    write("    --refIds\tReference sample IDs, separated by comma(,)", stdout())
    write("    --byGender\tCall by gender, default: FALSE", stdout())
    write("    --smallCNV\tFraction of CNV >=100Kb, default:0.5", stdout())
    write("    --bigCNV\tFraction of CNV >=5Mb, default:0.5", stdout())
    q()
} else {
    if (! modeHelp == "help") stop(paste0('Error modeHelp "', modeHelp, '"'))
    if (mode == "help" & modeHelp == "help") {
        write("Please input '--modeHelp modename' to see each mode's details.Mode name:\nsplit, cov, gender, merge, call", stdout())
        q()
    }
}

#' Import packages
library(CNVseq)

#' Performing
if (mode == "split") {
    batch_split(outPath, binSize, accessBed)
} else if (mode == "cov") {
    batch_cov(bamDirs, bamFiles, bedFile, outDir, thread, mapq)
} else if (mode == "gender") {
    batch_gender(cnnFiles, cnnDirs, outDir)
} else if (mode == "merge") {
    batch_merge(cnnFiles, cnnDirs, autoFile, sexFile, outPath)
} else if (mode == "call") {
    batch_call_cnv(cnmFile, cnnFiles, cnnDirs, outDir, testId, refIds, byGender, c(smallCNV, bigCNV))
} else {
    if (! mode == "help") stop(paste0('Error Mode "', mode, '"'))
}
