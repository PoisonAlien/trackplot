`trackplot.R` is a fast, simple, minimal dependency R script to generate IGV style track plots (aka locus plots) from bigWig files.


Example usage:

```r
source("trackplot.R")
trackplot(
  bigWigs = bws,
  loci = "chr3:187,715,903-187,752,003",
  draw_gene_track = TRUE,
  build = "hg38",
  scale_track_width = 2,
  nthreads = 4,
  mark_regions = data.frame(chr = "chr3", start = 187743255, end = 187747473),
  regions_track_width = 2,
  gene = "BCL6",
  gene_fsize = 1,
  custom_names = c("CD34", "EC", "LC", "CD4+", "CD8+")
)
```

<img src="example.png" /></a>