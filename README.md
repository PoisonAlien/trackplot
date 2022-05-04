## Introduction

`trackplot.R` is an ultra-fast, simple, and minimal dependency R script to generate IGV style track plots (aka locus plots) and profile plots from bigWig files. 
It has three main utilities:

* `track_extract()` -> `track_plot()` - for IGV style track plots
* `profile_extract()` -> `profile_plot()` - for density based profile plots
* `extract_signal()` -> `pca_plot()` - PCA analysis based on genomic regions of interest or around TSS sites of reference transcripts.

`trackplot.R` is a standalone R script and requires no installation. Just source it and you're good to go! See below for dependencies.

```r
source("https://github.com/PoisonAlien/trackplot/blob/master/R/trackplot.R?raw=true")

# OR

download.file(url = "https://raw.githubusercontent.com/PoisonAlien/trackplot/master/R/trackplot.R", destfile = "trackplot.R")
source('trackplot.R') 


# OR If you prefer to have it as package
remotes::install_github(repo = "poisonalien/trackplot")
```

Features:

  * It's significantly fast since most of the heavy lifting is done by `bwtool`. Below examples took less than 2 minutes on my 5 year old [macbook Pro](https://support.apple.com/kb/sp715?locale=en_GB) 
  * Automatically queries UCSC genome browser for gene models and cytobands, making analysis reproducible.
  * Supports GTF and standard UCSC gene formats as well.
  * Lightweight and minimal dependency 
    - [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and [bwtool](https://github.com/CRG-Barcelona/bwtool) are the only requirements. 
    - Plots are generated in pure base R graphics (no ggplot2 or tidyverse packages)
  * Customization: Each plot can customized for color, scale, width, etc.
  * Tracks can be summarized per condition (by mean, median, max, min)

### Usage

### trackplots

`track_extract()` and `track_plot()` are two functions to generate IGV style track plots (aka locus plots) from bigWig files. Additionally, `track_summarize` can summarize tracks by condition.
 
```r
#Path to bigWig files
bigWigs = c("CD34.bw", "EC.bw", "LC.bw", "CD4p.bw", "CD8p.bw")

#Step-1. Extract the siganl for your loci of interst
track_data = track_extract(bigWigs = bigWigs, loci = "chr3:187,715,903-187,752,003")

#Step-1a (optional). Summarize trcks by condition
track_data = track_summarize(summary_list = track_data, condition = c("A", "B", "B", "C", "D"), stat = "mean")

#Step-2. 
#Basic Plot 
track_plot(summary_list = track_data)

#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
track_plot(summary_list = track_data, draw_gene_track = TRUE, build = "hg38")

#With GTF file as source for gene models
track_plot(summary_list = track_data, draw_gene_track = TRUE, gene_model = "hg38_refseq.gtf.gz", isGTF = TRUE)

#Heighlight regions of interest

markregions = data.frame(
    chr = c("chr3", "chr3"),
    start = c(187743255, 187735888),
    end = c(187747473, 187736777),
    name = c("Promoter-1", "Promoter-2")
  )
  
track_plot(
  summary_list = track_data,
  draw_gene_track = TRUE,
  show_ideogram = TRUE,
  build = "hg38",
  regions = markregions
)
```

![](https://user-images.githubusercontent.com/8164062/101162153-3deae100-3632-11eb-8fad-66706f53ffe8.png)

### profileplot()

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

### bwpcaplot()

`bwpcaplot()` is a function to perform PCA analysis based on genomic regions of interest or around TSS sites of reference transcripts. Plot data and region summaries returned to the user.

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

```
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
