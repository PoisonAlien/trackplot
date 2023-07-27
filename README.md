## trackplot - Fast and easy visualisation of bigWig files in R

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Introduction

`trackplot.R` is an ultra-fast, simple, and minimal dependency R script to generate IGV style track plots (aka locus plots) and profile plots from bigWig files. 
It has three main utilities:

* `track_extract()` -> `track_plot()` - for IGV style track plots
* `profile_extract()` -> `profile_plot()` - for density based profile plots
* `extract_signal()` -> `pca_plot()` - PCA analysis based on genomic regions of interest or around TSS sites of reference transcripts.

## Installation

`trackplot.R` is a standalone R script and requires no installation. Just source it and you're good to go! See below for [dependencies](https://github.com/PoisonAlien/trackplot#dependencies).

```r
source("https://github.com/PoisonAlien/trackplot/blob/master/R/trackplot.R?raw=true")

# OR

download.file(url = "https://raw.githubusercontent.com/PoisonAlien/trackplot/master/R/trackplot.R", destfile = "trackplot.R")
source('trackplot.R') 

# OR If you prefer to have it as package

remotes::install_github(repo = "poisonalien/trackplot")
```

## Features

  * It's significantly fast since most of the heavy lifting is done by [bwtool](https://github.com/CRG-Barcelona/bwtool). Below examples took less than a minute on my 5 year old [macbook Pro](https://support.apple.com/kb/sp715?locale=en_GB) 
  * Automatically queries UCSC genome browser for gene models, cytobands, and chromHMM tracks - making analysis reproducible.
  * Supports GTF and standard UCSC gene formats as well.
  * Lightweight and minimal dependency 
    - [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and [bwtool](https://github.com/CRG-Barcelona/bwtool) are the only requirements. 
    - Plots are generated in pure base R graphics (no ggplot2 or tidyverse packages)
  * Customization: Each plot can be customized for color, scale, height, width, etc.
  * Tracks can be summarized per condition (by mean, median, max, min)

## Usage

## trackplots

`track_extract()` and `track_plot()` are two functions to generate IGV style track plots (aka locus plots) from bigWig files. Additionally, `track_summarize` can summarize tracks by condition.
 
### Step-1: Extract signal from bigWig files 
```r
#Path to bigWig files
bigWigs = c("H1_Oct4.bw", "H1_Nanog.bw", "H1_k4me3.bw", 
            "H1_k4me1.bw", "H1_k27ac.bw", "H1_H2az.bw", "H1_Ctcf.bw")

#Region to plot
oct4_loci = "chr6:31125776-31144789"

#Extract bigWig signal
t = track_extract(bigWigs = bigWigs, loci = oct4_loci, build = "hg19")
```


### Step-2: Plot

#### Basic plot
```r
track_plot(summary_list = t)
```

![](https://github.com/PoisonAlien/trackplot/assets/8164062/b50457b0-6623-47f6-b00f-f6e0de5a4808)

#### Add cytoband and change colors for each track
```r
track_plot(summary_list = t, 
          col = c("#d35400","#d35400","#27ae60","#27ae60","#2980b9","#2980b9","#2980b9"), 
          show_ideogram = TRUE)
```

![](https://github.com/PoisonAlien/trackplot/assets/8164062/a0911998-aae8-4de1-96f5-18e278d19d80)

#### Add TF binding sites at the top (any bed files would do)
```r
oct4_nanog_peaks = c("H1_Nanog.bed","H1_Oct4.bed") #Peak files 
track_plot(summary_list = t, 
          col = c("#d35400","#d35400","#27ae60","#27ae60","#2980b9","#2980b9","#2980b9"), 
          show_ideogram = TRUE, 
          peaks = oct4_nanog_peaks)
```

![](https://github.com/PoisonAlien/trackplot/assets/8164062/2531af5e-7200-478e-aa90-4ff5f537f57a)

#### Add some chromHMM tracks to the bottom

chromHMM data should be a bed file with the 4th column containing chromatin state. See here for an [example](https://github.com/PoisonAlien/trackplot/blob/master/inst/extdata/narrowpeak/H1_chromHMM.bed) file. 

Note that the color code for each of the 15 states are as described [here](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeBroadHmm&db=hg19). 
In case if it is different for your data, you will have to define your own color codes for each state and pass it to the argument `chromHMM_cols`

```r
chromHMM_peaks = "H1_chromHMM.bed"

track_plot(summary_list = t, 
          col = c("#d35400","#d35400","#27ae60","#27ae60","#2980b9","#2980b9","#2980b9"), 
          show_ideogram = TRUE, 
          peaks = oct4_nanog_peaks, chromHMM = chromHMM_peaks)
```
![](https://github.com/PoisonAlien/trackplot/assets/8164062/5ef8d09f-1bdf-4622-9367-4245bdec63d5)

#### Add some chromHMM tracks from UCSC

UCSC has 9 cell lines for which chromHMM data is available. These can be added automatically in case if you dont have your own data.
In this case, use the argument `ucscChromHMM` with any values from TableName column of the below table.

```r
                    TableName    cell                      Description              Tissue Karyotype
1: wgEncodeBroadHmmGm12878HMM GM12878     B-lymphocyte, lymphoblastoid               blood    normal
2:  wgEncodeBroadHmmH1hescHMM H1-hESC             embryonic stem cells embryonic stem cell    normal
3:   wgEncodeBroadHmmHepg2HMM   HepG2         hepatocellular carcinoma               liver    cancer
4:   wgEncodeBroadHmmHepg2HMM    HMEC         mammary epithelial cells              breast    normal
5:    wgEncodeBroadHmmHsmmHMM    HSMM        skeletal muscle myoblasts              muscle    normal
6:   wgEncodeBroadHmmHuvecHMM   HUVEC umbilical vein endothelial cells        blood vessel    normal
```

```r
track_plot(summary_list = t, 
          col = c("#d35400","#d35400","#27ae60","#27ae60","#2980b9","#2980b9","#2980b9"), 
          show_ideogram = TRUE, 
          peaks = oct4_nanog_peaks, 
          ucscChromHMM = c("wgEncodeBroadHmmH1hescHMM", "wgEncodeBroadHmmNhlfHMM"))
```

![](https://github.com/PoisonAlien/trackplot/assets/8164062/fecf7ab1-44cb-4308-b3f4-d8ca03cdd15d)


## narrowPeaks and broadPeaks 

All of the above plots can also be generated with [narroPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) or [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) files as input. Here, 5th column containing scores are plotted as intensity. Color coding and binning of scores are as per [UCSC convention](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

`narrowPeak` is one of the output from macs2 peak caller and are easier to visualize in the absence of bigWig files.

```r
narrowPeaks = c("H1_Ctcf.bed", "H1_H2az.bed", "H1_k27ac.bed", "H1_k4me1.bed", 
"H1_k4me3.bed", "H1_Nanog.bed", "H1_Oct4.bed", "H1_Pol2.bed")

oct4_loci = "chr6:30,818,383-31,452,182" #633Kb region for example

#Use the same track_extract but, instead of `bigWigs` use the `bed` argument
narrowPeaks_track = track_extract(bed = narrowPeaks, loci = oct4_loci, build = "hg19")

#Rest plotting is same
track_plot(summary_list = narrowPeaks_track, 
          show_ideogram = TRUE, 
          peaks = oct4_nanog_peaks, 
          ucscChromHMM = c("wgEncodeBroadHmmH1hescHMM", "wgEncodeBroadHmmNhlfHMM"))

```
![image](https://github.com/PoisonAlien/trackplot/assets/8164062/fa3999fd-ab7f-4617-a43e-d3cac7f3a3b3)


## profileplots

`profile_extract()` and `profile_plot()` are functions to generate density based profile-plots from bigWig files.

  * Below example for summarizing approx. 33,250 peaks for 5 bigWig files takes around 90 seconds on my 5 year old [macbook Pro](https://support.apple.com/kb/sp715?locale=en_GB). This includes generating signal matrix, summarizing, and plotting
  * Optionally, it can even query UCSC genome browser for refseq transcripts of desired assembly and summarize around TSS regions
  * Replicates can be collapsed into single value per condition

```r
#Example profile plot for a bed file with ~33,250 peaks, centered and extended 2500 bps
profile_data = profile_extract(
  bigWigs = bigWigs,
  bed = "CD34.bed",
  startFrom = "center",
  up = 2500,
  down = 2500
)

profile_plot(profile_data)

#profile plot for refseq protein-coding genes (TSS +/2500)
profile_data = profile_extract(
  bigWigs = bigWigs,
  ucsc_assembly = "hg38",
  startFrom = "tss",
  up = 2500,
  down = 2500
)
profile_plot(profile_data)
```

![](https://user-images.githubusercontent.com/8164062/100755019-05f25c80-33ec-11eb-900e-a9595d443f0f.png)

## PCA

`pca_plot()` is a function to perform PCA analysis based on genomic regions of interest or around TSS sites of reference transcripts. Plot data and region summaries returned to the user.

```r
#PCA using UCSC protein coding reference transcripts (TSS+/- 2500 bp)
refseq_summary = extract_summary(
  bigWigs = bigWigs,
  ucsc_assembly = "hg38"
)

pca_plot(summary_list = refseq_summary)

#PCA using genomic regions of interest (BED file)
bed_summary = extract_summary(
  bigWigs = bigWigs,
  bed = "sample.bed",
  custom_names = c("CD34", "EC", "LC", "CD4+", "CD8+")
)

pca_plot(summary_list = bed_summary)
```

![](https://user-images.githubusercontent.com/8164062/101655218-a62a3000-3a41-11eb-8d20-38d046d6f042.png)

### Dependencies

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) R package - which itself has no dependency.
* [bwtool](https://github.com/CRG-Barcelona/bwtool) - a command line tool for processing bigWig files. Install and move the binary to a PATH (e.g; `/usr/local/bin`). 
Or, you could also add the path where bwtool is located to R session with the below command.

```r
#Example
Sys.setenv(PATH = paste("/Users/anand/Documents/bwtool_dir/", Sys.getenv("PATH"), sep=":"))
```

* If you have trouble compiling the tool, follow [these](https://gist.github.com/PoisonAlien/e19b482ac6146bfb03142a0de1c4fbc8) instructions. Alternatively, you can download the pre-built binary for [macOS](https://www.dropbox.com/s/kajx9ya6erzyrim/bwtool_macOS.tar.gz?dl=1) or [centOS](https://www.dropbox.com/s/77ek89jqfhcmouu/bwtool_centOS_x86_64.tar.gz?dl=1)

***PSA*** If you find the tool useful, consider starrig this repository or upvoting this [Biostars thread](https://www.biostars.org/p/475853/) so that more poeple can find it :)

### Caveat

 * Windows OS is not supported
 
![](https://media.giphy.com/media/cKJjGbH7R5KKcJIR5u/giphy.gif)


### Citation

If you find the script useful consider [citing bwtool](https://academic.oup.com/bioinformatics/article/30/11/1618/282756)

*Pohl A, Beato M. bwtool: a tool for bigWig files. Bioinformatics. 2014 Jun 1;30(11):1618-9. doi: 10.1093/bioinformatics/btu056. Epub 2014 Jan 30. PMID: [24489365](https://pubmed.ncbi.nlm.nih.gov/24489365/); PMCID: PMC4029031.*

### Acknowledgements 

[Joschka Hey](https://github.com/HeyLifeHD) for all the cool suggestions :)
