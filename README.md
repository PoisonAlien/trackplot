## Introduction

`trackplot.R` is an ultra-fast, simple, and minimal dependency R script to generate IGV style track plots (aka locus plots) and profile plots from bigWig files. It has two main functions:

* `trackplot()` - for IGV style track plots
* `profileplot()` - for density based profile plots

`trackplot.R` is a standalone R script and requires no installation. Just source it and you're good to go! See below for dependencies.

### Usage

### trackplot()

`trackplot()` is a simple, and minimal dependency R function to generate IGV style track plots (aka locus plots) from bigWig files.

 * Its fast since most of the heavy lifting is done by `bwtool`. Above example plot took less than a minute on my 5 year old [macbook Pro](https://support.apple.com/kb/sp715?locale=en_GB) 
 * Automatically queries UCSC genome browser for gene models.
 * Supports GTF and standard UCSC gene formats as well.
 * Customization: Each track can be customized for color, scale, and width.
 * Minimal dependency. Plots are generated in pure base R graphics. 

```r
download.file(url = "https://raw.githubusercontent.com/PoisonAlien/trackplot/master/trackplot.R", destfile = "trackplot.R")
source('trackplot.R') 

#Path to bigWig files
bigWigs = c("CD34.bw", "EC.bw", "LC.bw", "CD4p.bw", "CD8p.bw")

#1. Basic usage
trackplot(bigWigs = bigWigs, loci = "chr3:187,715,903-187,752,003")

#2. With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
trackplot(bigWigs = bigWigs, loci = "chr3:187,715,903-187,752,003", draw_gene_track = TRUE, build = "hg38")

#3. With GTF file as source for gene models
trackplot(bigWigs = bigWigs, loci = "chr3:187,715,903-187,752,003", draw_gene_track = TRUE, gene_model = "hg38_refseq.gtf.gz", isGTF = TRUE)

#4. Heighlight regions of interest

## Example regions to heighlight (optional)
markregions = data.frame(
    chr = c("chr3", "chr3"),
    start = c(187743255, 187735888),
    end = c(187747473, 187736777),
    name = c("Promoter-1", "Promoter-2")
  )
  
trackplot(
  bigWigs = bigWigs,
  loci = "chr3:187,715,903-187,752,003",
  draw_gene_track = TRUE,
  build = "hg38",
  mark_regions = markregions,
  custom_names = c("CD34", "EC", "LC", "CD4+", "CD8+")
)
```

<img src="example.png" /></a>

Available arguments

```r
#' Generate IGV style locus tracks with ease
#' @param bigWigs bigWig files. Default NULL. Required.
#' @param loci target region to plot. Should be of format "chr:start-end". e.g; chr3:187715903-187752003 OR chr3:187,715,903-187,752,003
#' @param binsize bin size to extract signal. Default 50 (bps).
#' @param draw_gene_track Default FALSE. If TRUE plots gene models overlapping with the queried region
#' @param query_ucsc Default FALSE. But switches to TRUE when `gene_model` is not given. Requires `mysql` installation.
#' @param build Genome build. Default `hg19`
#' @param tx transcript name to draw. Default NULL. Plots all transcripts overlapping with the queried region
#' @param gene gene name to draw. Default NULL. Plots all genes overlapping with the queried region
#' @param collapse_tx Default FALSE. Whether to collapse all transcripts belonging to same gene into a unified gene model
#' @param gene_model File with gene models. Can be a gtf file or UCSC file format. If you have read them into R as a data.frame, that works as well. Default NULL, automatically fetches gene models from UCSC server
#' @param isGTF Default FALSE. Set to TRUE if the `gene_model` is a gtf file.
#' @param groupAutoScale Default TRUE
#' @param gene_fsize Font size. Default 1
#' @param gene_track_height Default 2 
#' @param scale_track_height Default 1
#' @param col Color for tracks. Default `#2f3640`. Multiple colors can be provided for each track
#' @param show_axis Default FALSE
#' @param custom_names Default NULL and Parses from the file names.
#' @param custom_names_pos Default 0 (corresponds to left corner)
#' @param mark_regions genomic regions to highlight. A data.frame with at-least three columns containing chr, start and end positions.
#' @param mark_regions_col color for highlighted region. Default "#192A561A"
#' @param mark_regions_col_alpha Default 0.5
#' @param nthreads Default 1. Number of threads to use.
```

### profileplot()

`profileplot()` is a simple, and minimal dependency R function to generate profile-plots from bigWig files.

  * Below example for summarizing approx. 33,250 peaks for 5 bigWig files takes around 90 seconds on my 5 year old [macbook Pro](https://support.apple.com/kb/sp715?locale=en_GB). This includes generating signal matrix, summarizing, and plotting
  * Optionally, it can even query UCSC genome browser for refseq transcripts of desired assembly and summarize around TSS regions
  * Replicates can be collapsed into single value per condition

```r
#Example profile plot for a bed file with ~33,250 peaks, centered and extended 2500 bps
profileplot(
  bigWigs = bigWigs,
  bed = "CD34.narrowPeak",
  startFrom = "center",
  up = 2500,
  down = 2500,
  custom_names = c("CD34", "EC", "LC", "CD4+", "CD8+")
)
```

![](https://user-images.githubusercontent.com/8164062/100755019-05f25c80-33ec-11eb-900e-a9595d443f0f.png)

Available arguments

```r
#' Generate profile plots with ease
#' @param bigWigs bigWig files. Default NULL. Required.
#' @param bed bed file or a data.frame with first 3 column containing chromosome, star, end positions. 
#' @param binSize bin size to extract signal. Default 50 (bps). Should be >1
#' @param startFrom Default "center". Can be "center", "start" or "end"
#' @param up extend upstream by this many bps from `startFrom`. Default 2500
#' @param down extend downstream by this many bps from `startFrom`. Default 2500
#' @param ucsc_assembly If `bed` file not provided, setting `ucsc_assembly` to ref genome build will fetch transcripts from UCSC genome browser. e.g; 'hg19'
#' @param nthreads Default 4
#' @param custom_names Default NULL and Parses from the file names.
#' @param color Manual colors for each bigWig. Default NULL. 
#' @param condition Default. Condition associated with each bigWig. Lines will colord accordingly.
#' @param condition_colors Manual colors for each level in condition. Default NULL. 
#' @param collapse_replicates Default FALSE. If TRUE and when `condition` is given, collapse signals samples belonging to same condition
#' @param plot_se Default FALSE. If TRUE plots standard error shading
#' @param line_size Default 1
#' @param legend_fs Legend font size. Default 1
#' @param axis_fs Axis font size. Default 1
```

### Dependencies

`trackplot` has only two dependencies. 

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) R package - which itself has no dependency.
* [bwtool](https://github.com/CRG-Barcelona/bwtool) - a command line tool for processing bigWig files. Install and move the binary to a PATH (e.g; `/usr/local/bin`). If you have trouble compiling the tool, follow [these](https://gist.github.com/PoisonAlien/e19b482ac6146bfb03142a0de1c4fbc8) instructions. Alternatively, you can download the pre-built binary for [macOS](https://www.dropbox.com/s/kajx9ya6erzyrim/bwtool_macOS.tar.gz?dl=1) or [centOS](https://www.dropbox.com/s/77ek89jqfhcmouu/bwtool_centOS_x86_64.tar.gz?dl=1)

### Caveat

 * Windows OS is not supported
 
![](https://media.giphy.com/media/cKJjGbH7R5KKcJIR5u/giphy.gif)


### Citation

If you find the script useful consider [citing bwtool](https://academic.oup.com/bioinformatics/article/30/11/1618/282756)

*Pohl A, Beato M. bwtool: a tool for bigWig files. Bioinformatics. 2014 Jun 1;30(11):1618-9. doi: 10.1093/bioinformatics/btu056. Epub 2014 Jan 30. PMID: [24489365](https://pubmed.ncbi.nlm.nih.gov/24489365/); PMCID: PMC4029031.*

### Acknowledgements 

[Joschka Hey](https://github.com/HeyLifeHD) for all the cool suggestions :)
