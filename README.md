# CNVseq
A tool for detecting CNVs in CNV-seq data. By default, it identifies DUP and DEL larger than 100 kb with a mosaicism level above 0.5, and CNVs larger than 5 Mb with a mosaicism level above 0.1. Special optimization is applied for detecting big CNVs and chromosomal abnormalities, ensuring improved detection performance (recommended value: 0.3, which can reduce detection caused by some data fluctuations). The mosaicism level is adjustable.

## Install
It can be installed from github:

```
devtools::install_github('sanadamakomi/CNVseq')
library(CNVseq)
```
## Run tool

`CNVseq.R` in `exec/` is a batch tool to run calling CNVs.

Before using you need to install package `GetoptLong`:

```
installed.packages("GetoptLong")
```

Please use `Rscript` or change to mode `0755` to run this tool. See more details by command:
```
# mode name: split, cov, gender, merge, call

Rscript CNVseq.R -modeHelp modename

# see mode 'cov'

Rscript CNVseq.R -modeHelp cov

# run mode 'cov'
Rscript CNVseq.R -mode cov -bamFiles bam_path -bedFile bed_path -outDir out_dir
```
## Pipeline to call CNVs

Recommend 5 steps which can run multiple BAM file from the same pool.

### Step 1. Get `target.bed` file

Split wgs region into bins.Recommend input an `access.bed`, which is a BED file include sequence-accessible region, e.g. access.bed of [cnvkit](https://github.com/etal/cnvkit/tree/master/data). You can also learn how to create this file from your own genome, [see](https://cnvkit.readthedocs.io/en/stable/pipeline.html#access).

Recommend a bin size of 20000 for 5-8G raw data, with a median read count of 150x in each window.

```
Rscript CNVseq.R -mode split -binSize 20000 -outPath target.bed

# input access.bed

Rscript CNVseq.R -mode split -binSize 20000 -outPath target.bed -accessBed access.bed
```

### Step 2. Get `*.cnn` file
Calculating read count in BAM files. You can use multiple threads to run faster. The `mapq` value will control which reads are counted; only those with a `mapq` value greater than the input threshold will be considered.

```
# input multiple BAM files
Rscript CNVseq.R -mode cov --bamFiles bam_1,bam_2,bam_3 -bedFile target.bed -outDir .  -thread 10

# Or input multiple directories storing BAM files
Rscript CNVseq.R -mode cov -bamDirs bam_dir_1,bam_dir_2,bam_dir_3 -bedFile target.bed -outDir .  -thread 10
```

### Step 3. Get `*.cnm` file
Merge multiple read count `*.cnn` file and creating a merged data frame to output.

```
Rscript CNVseq.R -mode merge -cnnFiles cnn_1,cnn_2,cnn_2 -outPath merge.cnm

# Or input a directory path of cnn files.

Rscript CNVseq.R -mode merge -cnnDirs cnn_dir_1,cnn_dir_2,cnn_dir_3 -outPath merge.cnm
```

### Step 4. Calling CNVs

Recommend running separately when the same gender sample count is less than 10. This allows calling CNVs in autosomes using more samples, and calling CNVs in sex chromosomes using the same gender samples. After calling, you can use the `merge` mode to combine the results. And the mosaicism level is adjustable, a level 0.5 will filter results with a copy number <=1.5 or >= 2.5, while 0.3 is <=1.7 or >= 2.3.
```
Rscript CNVseq.R -mode call -cnmFile merge.cnm -testId sample -outDir . -byGender -smallCNV 0.5 -bigCNV 0.3

# Build a baseline with samples of the same gender with the input sample id

Rscript CNVseq.R -mode call -cnmFile merge.cnm -testId sample -outDir . -byGender -smallCNV 0.5 -bigCNV 0.3 -byGender

```

### Step 5. Merge result
Merge VCF, CNR, CNS file of autosome and sex chromosome result. You can use the merged VCF to do annotation, or use CNR file to plot.
```
# merge VCF
Rscript CNVseq.R -mode merge -autoFile auto.cnv.vcf -sexFile auto.cnv.vcf -outPath merge.vcf

# merge log2ratio result CNR file from step 4
Rscript CNVseq.R -mode merge -autoFile auto.cnr -sexFile auto.cnr -outPath merge.cnr

# merge segmentation result CNS file from step 4
Rscript CNVseq.R -mode merge -autoFile auto.cns -sexFile auto.cns -outPath merge.cns
```

## Update
* 2024.12.05: Create package.
