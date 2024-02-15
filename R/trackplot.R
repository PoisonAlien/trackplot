# This R script contains functions for bigWig visualization
#
# Source code: https://github.com/PoisonAlien/trackplot
#
# MIT License
# Copyright (c) 2020 Anand Mayakonda <anandmt3@gmail.com>
#
# Change log:
# Version: 1.5.10 [2024-02-14]
#   * Added argument `layout_ord` and `bw_ord` to `track_plot()` re-order the overall tracks and bigWig tracks
#   * Added `xlab` and `ylab` arguments to `profile_plot()`
# Version: 1.5.01 [2023-10-17]
#   * Bug fix parsing loci while parsing GTF
#   * Small updates to profile_heatmap()
# Version: 1.5.00 [2023-08-24]
#   * Added `read_coldata` to import bigwig and bed files along with metadata. This streamlines all the downstream processes
#   * Added `profile_heatmap` for plotting heatmap
#   * Added `diffpeak` for minimal differential peak analysis based on peak intensities
#   * Added `volcano_plot` for diffpeak results visualization
#   * Added `summarize_homer_annots`
#   * Support for GTF files with `track_extract`
#   * `track_extract` now accepts gene name as input.
#   * More customization to `profile_extract` `profile_plot` and `plot_pca`
#   * Nicer output with `extract_summary` 
#   * Update mysql query for UCSC. Added `ideoTblName` argument for `track_extract`. Issue: #19
# Version: 1.4.00 [2023-07-27]
#   * Updated track_plot to include chromHMM tracks and top peaks tracks
#   * Support to draw narrowPeak or boradPeak files with track_plot
#   * Support to query ucsc for chromHMM tracks
#   * Additional arguments to track_plot to adjust heights of all the tracks and margins
#   * Improved track_extract - (extracts gene models and cytobands to avoid repetitive calling ucsc genome browser)
#   * Additional arguments to pca_plot for better plotting
#   * Added example datasets
# Version: 1.3.10 [2021-10-06]
#   * Support for negative values (Issue: https://github.com/PoisonAlien/trackplot/issues/6 )
#   * Added y_min argument to track_plot. 
#   * Change the default value for collapse_tx to TRUE
# Version: 1.3.05 [2021-06-07]
#   * Summarize and groupScaleByCondition tracks by condition. Issue: #4
#   * Allow the script to install as a package.
#   * Added y_max argument for custom y-axis limits in track_plot. 
# Version: 1.3.01 [2021-04-26]
#   * Fix gtf bug. Issue: #3
# Version: 1.3.0 [2021-03-26]
#   * modularize the code base to avoid repetitive data extraction and better plotting
# Version: 1.2.0 [2020-12-09]
#   * Added bwpcaplot()
# Version: 1.1.11 [2020-12-07]
#   * Bug fixes in profileplot(): Typo for .check_dt() and startFrom estimation
# Version: 1.1.1 [2020-12-04]
#   * trackplot() now plots ideogram of target chromosome
# Version: 1.1.0 [2020-12-01]
#   * Added profileplot()
# Version: 1.0.0 [2020-11-27]
#   * Initial release

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'Prepares meta data table from bigWig files.
#'Output from this function is passed to all downstream functions.
#' @param bws path to bigWig files
#' @param sample_names sample names for each input files. Optional. Default NULL - creates one from file names.
#' @param build Reference genome build. Default hg38
#' @param input_type Default `bw`. Can be `bw` or `peak`
#' @examples
#' bigWigs = system.file("extdata", "bw", package = "trackplot") |> list.files(pattern = "\\.bw$", full.names = TRUE) 
#' cd = read_coldata(bws = bigWigs, build = "hg19")
#' beds = system.file("extdata", "narrowpeak", package = "trackplot") |> list.files(pattern = "\\.bed$", full.names = TRUE) 
#' cd_bed = read_coldata(bws = beds, input_type = "peak", build = "hg19")
#' @export

read_coldata = function(bws = NULL, sample_names = NULL, build = "hg38", input_type = "bw"){
  
  if(is.null(bws)){
    stop("Please provide paths to bigWig files")
  }
  
  input_type = match.arg(arg = input_type, choices = c("bw", "peak"))
  
  message("Checking for files..")
  bws = as.character(bws)
  lapply(bws, function(x){
    if(!file.exists(x)){
      stop(paste0(x, " does not exist!"))
    }
  })
  
  if(is.null(sample_names)){
    bw_sample_names = unlist(data.table::tstrsplit(x = basename(bws), split = "\\.", keep = 1))
  }else{
    bw_sample_names = as.character(sample_names)  
  }
  
  if(any(duplicated(bw_sample_names))){
    stop("Found duplicates. Samples names must be unique")
  }
  if(length(bw_sample_names) != length(bws)){
    stop("Please provide names for each input file")
  }
  
  coldata = data.table::data.table(bw_files = bws, bw_sample_names = bw_sample_names)
  
  attr(coldata, "refbuild") = build
  attr(coldata, "is_bw") =  input_type ==  "bw"
  message("Input type: ", input_type)
  message("Ref genome: ", build)
  message("OK!")
  
  coldata
}

#' Extract bigWig track data for the given loci
#' @param colData coldata from \code{read_coldata}
#' @param loci target region to plot. Should be of format "chr:start-end". e.g; chr3:187715903-187752003 OR chr3:187,715,903-187,752,003
#' @param gene gene name. This is mutually exclusive with \code{loci}
#' @param binsize bin size to extract signal. Default 10 (bps).
#' @param nthreads Default 1. Number of threads to use.
#' @param query_ucsc Default TRUE. Queries UCSC and extracts gene models and cytoband for the loci. Requires `mysql` installation.
#' @param gtf Use gtf file or data.frame as source for gene model. Default NULL.
#' @param build Reference genome build. Default hg38
#' @param padding Extend locus on both sides by this many bps. 
#' @param ideoTblName Table name for ideogram. Default `cytoBand`
#' @import data.table
#' @examples
#' bigWigs = system.file("extdata", "bw", package = "trackplot") |> list.files(pattern = "\\.bw$", full.names = TRUE) 
#' cd = read_coldata(bws = bigWigs, build = "hg19")
#' oct4_loci = "chr6:31125776-31144789"
#' t = track_extract(colData = cd, loci = oct4_loci, build = "hg19")
#' @export
track_extract = function(colData = NULL, loci = NULL, gene = NULL, binsize = 10, nthreads = 1, query_ucsc = TRUE, gtf = NULL, build = "hg38", padding = 0, ideoTblName = "cytoBand"){
  
  if(is.null(colData)){
    stop("Missing colData. Use read_coldata() to generate one.")
  }
  
  if(all(is.null(loci), is.null(gene))){
    stop("Please provide a loci or a gene name!")
  }
  
  input_bw = attr(colData, "is_bw")
  build = attr(colData, "refbuild")

  .check_windows()
  if(input_bw){
    .check_bwtool()  
  }
  
  .check_dt()
  
  options(warn = -1)
  op_dir = tempdir() #For now
  
  if(!is.null(gene)){
    if(!is.null(gtf)){
      etbl = .parse_gtf(gtf = gtf, genename = gene)
      cyto = NA
      start = min(unlist(lapply(etbl,function(x) attr(x, "start"))))
      end = max(unlist(lapply(etbl,function(x) attr(x, "end"))))
      chr = unique(unlist(lapply(etbl,function(x) attr(x, "chr"))))
    }else if(query_ucsc){
      message("Querying UCSC genome browser for gene model and cytoband..")
      etbl = .extract_geneModel_ucsc_bySymbol(genesymbol = gene, refBuild = build)
      chr = unique(as.character(etbl$chr)); start = min(as.numeric(etbl$start)); end = max(as.numeric(etbl$end))
      if(length(chr) > 1){
        message("Multiple chromosomes found! Using the first one ", chr[1])
        etbl = etbl[chr %in% chr[1]]
        chr = unique(as.character(etbl$chr)); start = min(as.numeric(etbl$start)); end = max(as.numeric(etbl$end))
      }
      if(!is.null(etbl)){
        etbl = .make_exon_tbl(gene_models = etbl)  
      }
      cyto = .extract_cytoband(chr = chr, refBuild = build)
      loci = paste0(chr, ":", start, "-", end)
    }else{
      cyto = etbl = NA
    }
    if(is.null(etbl)){
      stop("No transcript models found for ", gene)
    }
  }else{
    message("Parsing loci..")
    loci_p = .parse_loci(loci = loci)
    chr = loci_p$chr; start = loci_p$start; end = loci_p$end
    if(start >= end){
      stop("End must be larger than Start!")
    }
    message("    Queried region: ", chr, ":", start, "-", end, " [", end-start, " bps]")
    #Extract gene models for this region
    if(!is.null(gtf)){
      etbl = .parse_gtf(gtf = gtf, chr = chr, start = start, end = end)
      cyto = NA
    }else if(query_ucsc){
      message("Querying UCSC genome browser for gene model and cytoband..")
      etbl = .extract_geneModel_ucsc(chr, start = start, end = end, refBuild = build, txname = NULL, genename = NULL)
      if(!is.null(etbl)){
        etbl = .make_exon_tbl(gene_models = etbl)  
      }
      cyto = .extract_cytoband(chr = chr, refBuild = build, tblName = ideoTblName)
    }else{
      cyto = etbl = NA
    }
  }
  
  start = start - as.numeric(padding)  
  end = end + as.numeric(padding)  
  loci = paste0(chr, ":", start, "-", end)
  
  input_files = colData$bw_files
  custom_names = colData$bw_sample_names
  
  if(input_bw){
    windows = .gen_windows(chr = chr, start = start, end = end, window_size = binsize, op_dir = op_dir)
    track_summary = .get_summaries(bedSimple = windows, bigWigs = input_files, op_dir = op_dir, nthreads = nthreads)  
  }else{
    track_summary = .get_summaries_narrowPeaks(bigWigs = input_files, nthreads = nthreads, chr, start, end)  
  }
  
  names(track_summary) = custom_names
  
  attr(track_summary, "meta") = list(etbl = etbl, cyto = cyto, loci = loci)
  message("OK!")
  
  list(data = track_summary, colData = colData)
}

#' Summarize tracks per condition
#' @param summary_list Output from track_extract. Required.
#' @param condition a column name in \code{coldata} containing sample conditions. Default NULL.
#' @param stat can be `mean, median`, `max`, `min`. NAs are excluded. 
#' @export
track_summarize = function(summary_list = NULL, condition = NULL, stat = "mean"){
  
  if(is.null(summary_list)){
    stop("Missing input! Expecting output from track_extract()")
  }
  
  stat = match.arg(arg = stat, choices = c("mean", "median", "max", "min"))
  
  meta = attr(summary_list$data, "meta")
  loci = meta$loci
  etbl = meta$etbl
  cyto = meta$cyto
  
  coldata = summary_list$colData
  is_bw = attr(coldata, "is_bw")
  build = attr(coldata, "refbuild")
  
  if(is.null(condition)){
    stop("Please provide a column name containing sample condition!\nHere are available columns.\n", paste(colnames(coldata), collapse = " "))
  }
  
  summary_list = summary_list$data
  
  if(!is.null(condition)){
    if(!condition %in% colnames(coldata)){
      warning(paste0(condition, " does not exists in coldata. Here are available columns."))
      print(coldata)
      stop()
    }else{
      colnames(coldata)[which(colnames(coldata) == condition)] = "group_condition"
      coldata$group_condition = as.character(coldata$group_condition)
    }
    condition = as.character(coldata$group_condition)
  }
  
  names(summary_list) = condition
  
  summary_list = data.table::rbindlist(l = summary_list, use.names = TRUE, fill = TRUE, idcol = "sample_name")
  
  if(stat == "mean"){
    summary_list = summary_list[,mean(max, na.rm = TRUE), .(sample_name, chromosome, start, end)]  
  }else if (stat == "median"){
    summary_list = summary_list[,median(max, na.rm = TRUE), .(sample_name, chromosome, start, end)]
  }else if (stat == "max"){
    summary_list = summary_list[,max(max, na.rm = TRUE), .(sample_name, chromosome, start, end)]
  }else{
    summary_list = summary_list[,min(max, na.rm = TRUE), .(sample_name, chromosome, start, end)]
  }
  
  colnames(summary_list)[ncol(summary_list)] = "max" #this column name means nothing, just using it for the consistency
  summary_list = split(summary_list, summary_list$sample_name)
  attr(summary_list, "meta") = meta
  
  list(data = summary_list, colData = coldata)
}

#' Generate IGV style locus tracks with ease
#' @param summary_list Output from track_extract
#' @param draw_gene_track Default FALSE. If TRUE plots gene models overlapping with the queried region
#' @param show_ideogram Default TRUE. If TRUE plots ideogram of the target chromosome with query loci highlighted. Works only when `query_ucsc` is TRUE. 
#' @param txname transcript name to draw. Default NULL. Plots all transcripts overlapping with the queried region
#' @param genename gene name to draw. Default NULL. Plots all genes overlapping with the queried region
#' @param collapse_txs Default FALSE. Whether to collapse all transcripts belonging to same gene into a unified gene model
#' @param groupAutoScale Default TRUE
#' @param y_max custom y axis upper limits for each track. Recycled if required.
#' @param y_min custom y axis lower limits for each track. Recycled if required.
#' @param gene_fsize Font size. Default 1
#' @param col Color for tracks. Default `#2f3640`. Multiple colors can be provided for each track
#' @param show_axis Default FALSE
#' @param track_names Default NULL
#' @param track_names_pos Default 0 (corresponds to left corner)
#' @param track_names_to_left If TRUE, track names are shown to the left of the margin. Default FALSE, plots on top as a title
#' @param regions genomic regions to highlight. A data.frame with at-least three columns containing chr, start and end positions.
#' @param boxcol color for highlighted region. Default "#192A561A"
#' @param boxcolalpha Default 0.5
#' @param ucscChromHMM Name of the chromHMM table. Use .get_ucsc_hmm_tbls() to see the details.
#' @param chromHMM chromHMM data. Can be path to bed files or a list data.frames with first three columns containing chr,start,end and a 4th column containing integer coded state
#' @param chromHMM_names name for the chromHMM track
#' @param chromHMM_cols A named vector for each state (in the 4th column of chromHMM file). Default NULL
#' @param peaks bed file to be highlighted. Can be path to bed files or a list data.frames with first three columns containing chr,start,end.
#' @param peaks_track_names Provide a name for each loci bed file. Default NULL
#' @param cytoband_track_height Default 1
#' @param chromHMM_track_height Default 1
#' @param gene_track_height Default 2 
#' @param scale_track_height Default 1
#' @param peaks_track_height Default 2.
#' @param bw_track_height Default 3
#' @param left_mar Space to the left. Default 4
#' @param bw_ord Names of the tracks to be drawn in the provided order. Default NULL.
#' @param layout_ord Plot layout order. Deafult c("p", "b", "h", "g", "c") corresponding to peaks track, bigWig track, chromHmm track, gene track, cytoband track.
#' @examples
#' bigWigs = system.file("extdata", "bw", package = "trackplot") |> list.files(pattern = "\\.bw$", full.names = TRUE) 
#' cd = read_coldata(bws = bigWigs, build = "hg19")
#' oct4_loci = "chr6:31125776-31144789"
#' t = track_extract(colData = cd, loci = oct4_loci, build = "hg19")
#' trackplot::track_plot(summary_list = t)
#' @export
track_plot = function(summary_list = NULL,
                      draw_gene_track = TRUE,
                      show_ideogram = TRUE,
                      col = "gray70",
                      groupAutoScale = FALSE,
                      y_max = NULL,
                      y_min = NULL,
                      txname = NULL,
                      genename = NULL,
                      show_axis = FALSE,
                      gene_fsize = 1,
                      track_names = NULL,
                      track_names_pos = 0,
                      track_names_to_left = FALSE,
                      regions = NULL,
                      collapse_txs = TRUE,
                      boxcol = "#192A561A",
                      boxcolalpha = 0.2,
                      chromHMM = NULL,
                      chromHMM_cols = NULL,
                      chromHMM_names = NULL,
                      ucscChromHMM = NULL,
                      peaks = NULL,
                      bw_track_height = 3,
                      peaks_track_height = 2,
                      gene_track_height = 2,
                      scale_track_height = 2,
                      chromHMM_track_height = 1,
                      cytoband_track_height = 2,
                      peaks_track_names = NULL,
                      left_mar = NULL,
                      bw_ord = NULL,
                      layout_ord = c("p", "b", "h", "g", "c")
){
  
  if(is.null(summary_list)){
    stop("Missing input! Expecting output from track_extract()")
  }
  
  
  meta = attr(summary_list$data, "meta")
  loci = meta$loci
  etbl = meta$etbl
  cyto = meta$cyto
  
  coldata = summary_list$colData
  is_bw = attr(coldata, "is_bw")
  build = attr(coldata, "refbuild")
  
  loci_p = .parse_loci(loci = loci)
  chr = loci_p$chr; start = loci_p$start; end = loci_p$end
  
  # chr = summary_list$loci[1]
  # start = as.numeric(summary_list$loci[2])
  # end = as.numeric(summary_list$loci[3])
  # etbl = summary_list$etbl
  # cyto = summary_list$cyto
  # is_bw = attr(summary_list, "is_bw")
  # build = attr(summary_list, "refbuild")
  
  summary_list = summary_list$data
  
  #Change the order
  if(!is.null(bw_ord)){
    bw_ord = intersect(names(summary_list), bw_ord)
    
    if(length(bw_ord) == 0){
      stop("None of the provided bw_ord are presnt in the data! Available names:\n", paste(names(summary_list), collapse = ", "))
    }
    
    summary_list = summary_list[bw_ord]
    coldata = data.table::rbindlist(split(coldata, coldata$bw_sample_names)[bw_ord])
  }
  
  if(length(col) != length(summary_list)){
    col = rep(x = col, length(summary_list))
  }
  
  plot_regions = FALSE
  if(!is.null(regions)){
    if(is(object = regions, class2 = "data.frame")){
      regions = data.table::as.data.table(x = regions)
      colnames(regions)[1:3] = c("chromsome", "startpos", "endpos")
      regions = regions[chromsome %in% chr]
      if(nrow(regions) == 0){
        warning("None of the regions are within the requested chromosme: ", chr)
        plot_regions = TRUE
      }else{
        plot_regions = TRUE  
      }
    }else{
      stop("'mark_regions' must be a data.frame with first 3 columns containing : chr, start, end")
    }
  }
  
  if(!is.null(track_names)){
    names(summary_list) = track_names
  }
  
  groupScaleByCondition = FALSE #For furture
  if(groupScaleByCondition){
    plot_height = unlist(lapply(summary_list, function(x) max(x$max, na.rm = TRUE)))
    plot_height_min = unlist(lapply(summary_list, function(x) min(x$max, na.rm = TRUE)))
    plot_height = data.table::data.table(plot_height, plot_height_min, col, names(summary_list))
    plot_height$og_ord = 1:nrow(plot_height)
    plot_height = plot_height[order(col)]
    plot_height_max = plot_height[,.(.N, max(plot_height)), .(col)]
    plot_height_min = plot_height[,.(.N, max(plot_height_min)), .(col)]
    plot_height$max = rep(plot_height_max$V2, plot_height_max$N)
    plot_height$min = rep(plot_height_min$V2, plot_height_min$N)
    plot_height_min = plot_height[order(og_ord)][,min]
    plot_height = plot_height[order(og_ord)][,max]
  }else if(groupAutoScale){
    plot_height = max(unlist(lapply(summary_list, function(x) max(x$max, na.rm = TRUE))), na.rm = TRUE)
    plot_height_min = min(unlist(lapply(summary_list, function(x) min(x$max, na.rm = TRUE))), na.rm = TRUE)
    plot_height = rep(plot_height, length(summary_list))
    plot_height_min = rep(plot_height_min, length(summary_list))
  }else{
    plot_height = unlist(lapply(summary_list, function(x) max(x$max, na.rm = TRUE)))
    plot_height_min = unlist(lapply(summary_list, function(x) min(x$max, na.rm = TRUE)))
  }
  
  if(!is.null(y_max)){
    #If custom ylims are provided
    if(length(y_max) != length(summary_list)){
      y_max = rep(y_max, length(summary_list))
    }
    plot_height = y_max
    
  }else{
    plot_height = round(plot_height, digits = 2)  
  }
  
  if(!is.null(y_min)){
    #If custom ylims are provided
    if(length(y_min) != length(summary_list)){
      y_min = rep(y_min, length(summary_list))
    }
    plot_height_min = y_min
    
  }else{
    plot_height_min = round(plot_height_min, digits = 2)  
  }
  
  ntracks = length(summary_list)
  
  lo = .make_layout(ntracks = ntracks, ntracks_h = bw_track_height, cytoband = show_ideogram, cytoband_h = cytoband_track_height, genemodel = draw_gene_track, 
               genemodel_h = gene_track_height, chrHMM = any(!is.null(ucscChromHMM), !is.null(chromHMM)), chrHMM_h = chromHMM_track_height, loci = !is.null(peaks), 
               loci_h = peaks_track_height, scale_track_height = scale_track_height, lord = layout_ord)
  
  query = data.table::data.table(chr = chr, start = start, end = end)
  data.table::setkey(x = query, chr, start, end)
  
  if(is.null(left_mar)){
    left_mar = ifelse(test = show_axis, yes = 4, no = 2)
  }
  
  #Draw top peaks
  if(!is.null(peaks)){
    
    if(is.list(peaks)){
      peaks_data = lapply(peaks, function(l){
        colnames(l)[1:3] = c("chr", "start", "end")
        data.table::setDT(l, key = c("chr", "start", "end"))
        l
      })
    }else{
      peaks_data = lapply(peaks, function(l){
        l = data.table::fread(file = l)
        colnames(l)[1:3] = c("chr", "start", "end")
        data.table::setDT(l, key = c("chr", "start", "end"))
        l
      })
    }
    
    if(is.null(peaks_track_names)){
      names(peaks_data) = paste0("Bed", 1:length(peaks_data))  
    }else{
      names(peaks_data) = peaks_track_names
    }
    
    if(show_axis){
      par(mar = c(0.25, left_mar, 0.25, 1))
    }else{
      par(mar = c(0.25, left_mar, 0.25, 1))  
    }
    
    plot(NA, xlim = c(start, end), ylim = c(0, length(peaks)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
    
    for(idx in seq_along(peaks_data)){
      l_idx = peaks_data[[idx]]
      l_idx = data.table::foverlaps(x = query, y = l_idx, type = "any", nomatch = NULL)[,.(chr, start, end)]
      rect(xleft = start, ybottom = idx - 0.49, xright = end, ytop = idx - 0.51, col = "gray90", border = NA)
      if(nrow(l_idx) > 0){
        rect(xleft = l_idx$start, ybottom = idx - 0.9, xright = l_idx$end, ytop = idx - 0.1, col = "#34495e", border = NA)
      }
      text(x = start, y = idx - 0.5, labels = names(peaks_data)[idx], adj = 1.2, xpd = TRUE)
    }
  }
  
  #Draw bigWig signals
  if(is_bw){
    for(idx in 1:length(summary_list)){
      x = summary_list[[idx]]
      if(show_axis){
        par(mar = c(0.5, left_mar, 2, 1))
      }else{
        par(mar = c(0.5, left_mar, 2, 1))  
      }
      
      plot(NA, xlim = c(start, end), ylim = c(plot_height_min[idx], plot_height[idx]), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      #If there is no signal, just add the track names and go to next
      if(nrow(x) == 0){
        if(track_names_to_left){
          text(x = start, y = 0.5, labels = names(summary_list)[idx], adj = 1, cex = gene_fsize, xpd = TRUE)
          #mtext(text = names(summary_list)[idx], side = 2, line = -2, outer = TRUE, xpd = TRUE, las = 2, adj = 0)
          #title(main = , adj = track_names_pos, font.main = 3)  
        }else{
          title(main = names(summary_list)[idx], adj = track_names_pos, font.main = 3)
        }
        next
      }
      rect(xleft = x$start, ybottom = 0, xright = x$end, ytop = x$max, col = col[idx], border = col[idx])
      if(show_axis){
        axis(side = 2, at = c(plot_height_min[idx], plot_height[idx]), las = 2)  
      }else{
        text(x = start, y = plot_height[idx], labels = paste0("[", plot_height_min[idx], "-", plot_height[idx], "]"), adj = 0, xpd = TRUE)
      }
      #plot(NA, xlim = c(start, end), ylim = c(0, nrow(regions)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      
      if(plot_regions){
        # boxcol = "#192a56"
        boxcol = grDevices::adjustcolor(boxcol, alpha.f = boxcolalpha)
        if(nrow(regions) > 0){
          for(i in 1:nrow(regions)){
            if(idx == length(summary_list)){
              #If its a last plot, draw rectangle till 0
              rect(xleft = regions[i, startpos], ybottom = 0, xright = regions[i, endpos], ytop = plot_height[idx]+10, col = boxcol, border = NA, xpd = TRUE)
            }else if (idx == 1){
              if(ncol(regions) > 3){
                text(x = mean(c(regions[i, startpos], regions[i, endpos])), y = plot_height[idx]+(plot_height[idx]*0.1), labels = regions[i, 4], adj = 0.5, xpd = TRUE, font = 1, cex = 1.2)
              }else{
                text(x = mean(c(regions[i, startpos], regions[i, endpos])), y = plot_height[idx]+(plot_height[idx]*0.1), labels = paste0(regions[i, startpos], "-", regions[i, endpos]), adj = 0.5, xpd = TRUE, font = 1, cex = 1.2)
              }
              rect(xleft = regions[i, startpos], ybottom = -10, xright = regions[i, endpos], ytop = plot_height[idx], col = boxcol, border = NA, xpd = TRUE)
            }else{
              rect(xleft = regions[i, startpos], ybottom = -10, xright = regions[i, endpos], ytop = plot_height[idx]+10, col = boxcol, border = NA, xpd = TRUE)  
            }
          }
        }
      }
      
      if(track_names_to_left){
        text(x = start, y = (plot_height_min[idx] + plot_height[idx])/2, labels = names(summary_list)[idx], adj = 1.1, cex = gene_fsize, xpd = TRUE)
      }else{
        title(main = names(summary_list)[idx], adj = track_names_pos, font.main = 3)  
      }
    }
  }else{
    for(idx in 1:length(summary_list)){
      x = summary_list[[idx]]
      
      if(show_axis){
        par(mar = c(0, left_mar, 0, 1))
      }else{
        par(mar = c(0.5, left_mar, 1, 1))  
      }
      
      plot(NA, xlim = c(start, end), ylim = c(0, 1), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      #If there is no signal, just add the track names and go to next
      if(nrow(x) == 0){
        if(track_names_to_left){
          text(x = start, y = 0.5, labels = names(summary_list)[idx], adj = 1, cex = gene_fsize, xpd = TRUE)
          #mtext(text = names(summary_list)[idx], side = 2, line = -2, outer = TRUE, xpd = TRUE, las = 2, adj = 0)
          #title(main = , adj = track_names_pos, font.main = 3)  
        }else{
          title(main = names(summary_list)[idx], adj = track_names_pos, font.main = 3)
        }
        next
      }
      cols = cut(x$max, breaks = c(0, 166, 277, 389, 500, 612, 723, 834, 945, max(x$max)), labels = c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", 
                                                                                                      "#525252", "#252525", "#000000"))
      rect(xleft = x$start, ybottom = 0.01, xright = x$end, ytop = 0.99, col = as.character(cols), border = NA)
      
      if(plot_regions){
        # boxcol = "#192a56"
        boxcol = grDevices::adjustcolor(boxcol, alpha.f = boxcolalpha)
        if(nrow(regions) > 0){
          for(i in 1:nrow(regions)){
            if(idx == length(summary_list)){
              #If its a last plot, draw rectangle till 0
              rect(xleft = regions[i, startpos], ybottom = 0, xright = regions[i, endpos], ytop = plot_height[idx]+10, col = boxcol, border = NA, xpd = TRUE)
            }else if (idx == 1){
              if(ncol(regions) > 3){
                text(x = mean(c(regions[i, startpos], regions[i, endpos])), y = plot_height[idx]+(plot_height[idx]*0.1), labels = regions[i, 4], adj = 0.5, xpd = TRUE, font = 1, cex = 1.2)
              }else{
                text(x = mean(c(regions[i, startpos], regions[i, endpos])), y = plot_height[idx]+(plot_height[idx]*0.1), labels = paste0(regions[i, startpos], "-", regions[i, endpos]), adj = 0.5, xpd = TRUE, font = 1, cex = 1.2)
              }
              rect(xleft = regions[i, startpos], ybottom = -10, xright = regions[i, endpos], ytop = plot_height[idx], col = boxcol, border = NA, xpd = TRUE)
            }else{
              rect(xleft = regions[i, startpos], ybottom = -10, xright = regions[i, endpos], ytop = plot_height[idx]+10, col = boxcol, border = NA, xpd = TRUE)  
            }
          }
        }
      }
      
      if(track_names_to_left){
        text(x = start, y = 0.5, labels = names(summary_list)[idx], adj = 1, cex = gene_fsize, xpd = TRUE)
        #mtext(text = names(summary_list)[idx], side = 2, line = -2, outer = TRUE, xpd = TRUE, las = 2, adj = 0)
        #title(main = , adj = track_names_pos, font.main = 3)  
      }else{
        title(main = names(summary_list)[idx], adj = track_names_pos, font.main = 3)
      }
    }
    
  }
  
  #Draw chrom HMM tracks
  plotHMM = FALSE
  if(!is.null(chromHMM)){
    if(is.list(chromHMM)){
      chromHMM = lapply(chromHMM, function(l){
        colnames(l)[1:4] = c("chr", "start", "end", "name")
        data.table::setDT(l, key = c("chr", "start", "end"))
        l
      })
    }else{
      chromHMM = lapply(chromHMM, function(l){
        l = data.table::fread(file =  l)
        colnames(l)[1:4] = c("chr", "start", "end", "name")
        data.table::setDT(l, key = c("chr", "start", "end"))
      })
    }
    
    hmmdata = lapply(chromHMM, function(hmm){
      .load_chromHMM(chr = chr, start = start, end = end, ucsc = hmm)
      #.extract_chromHmm_ucsc()
    })
    
    if(is.null(chromHMM_names)){
      names(hmmdata) = paste0("chromHMM_", 1:length(hmmdata))
    }else{
      names(hmmdata) = chromHMM_names
    }
    
    plotHMM = TRUE
  }else if(!is.null(ucscChromHMM)){
    hmmdata = lapply(ucscChromHMM, function(hmmtbl){
      .extract_chromHmm_ucsc(chr = chr, start = start, end = end, refBuild = build, tbl = hmmtbl)
    })
    names(hmmdata) = ucscChromHMM
    plotHMM = TRUE
    #return(hmmdata)
  }
  
  if(plotHMM){
    if(show_axis){
      par(mar = c(0.1, left_mar, 0, 1))
    }else{
      par(mar = c(0.1, left_mar, 0, 1))  
    }
    
    if(is.null(chromHMM_cols)){
      chromHMM_cols = .get_ucsc_hmm_states_cols()
    }
    
    .plot_ucsc_chrHmm(d = hmmdata, start = start, end = end, hmm_cols = chromHMM_cols)
  }
  
  #Draw gene models
  if(draw_gene_track){
    
    #etbl = .make_exon_tbl(gene_models = etbl, txname = txname, genename = genename)
    
    if(!is.null(etbl)){
      
      if(!is.null(genename)){
        if(length(etbl[unlist((lapply(etbl, attr, "gene"))) %in% genename]) == 0){
          message("Note: Could not find any of the requested gene names! Available genes are:")
          print(unique(unlist((lapply(etbl, attr, "gene")))))
        }else{
          etbl = etbl[unlist((lapply(etbl, attr, "gene"))) %in% genename]  
        }
      }
      
      if(collapse_txs){
        etbl = .collapse_tx(etbl)
      }
      
      if(show_axis){
        par(mar = c(0.25, left_mar, 0, 1))
      }else{
        par(mar = c(0.25, left_mar, 0, 1))  
      }
      
      plot(NA, xlim = c(start, end), ylim = c(0, length(etbl)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      exon_col = "#192a56"
      for(tx_id in 1:length(etbl)){
        txtbl = etbl[[tx_id]]
        segments(x0 = attr(txtbl, "start"), y0 = tx_id-0.45, x1 = attr(txtbl, "end"), y1 = tx_id-0.45, col = exon_col, lwd = 1)
        name_at = min(c(txtbl[[1]], txtbl[[2]]))
        if(is.na(attr(txtbl, "tx"))){
          text(x = name_at, y = tx_id-0.45, labels = paste0(attr(txtbl, "gene")), adj = 1, cex = gene_fsize, xpd = TRUE, pos = 2) #x = start for outer margin
        }else{
          text(x = name_at, y = tx_id-0.45, labels = paste0(attr(txtbl, "tx"), " [", attr(txtbl, "gene"), "]"), cex = gene_fsize, adj = 0, xpd = TRUE, pos = 2)  
        }
        
        rect(xleft = txtbl[[1]], ybottom = tx_id-0.75, xright = txtbl[[2]], ytop = tx_id-0.25, col = exon_col, border = NA)
        if(attr(txtbl, "strand") == "+"){
          dirat = pretty(x = c(min(txtbl[[1]]), max(txtbl[[2]])))
          dirat[1] = min(txtbl[[1]]) #Avoid drawing arrows outside gene length
          dirat[length(dirat)] = max(txtbl[[2]])
          points(x = dirat, y = rep(tx_id-0.45, length(dirat)), pch = ">", col = exon_col)
        }else{
          dirat = pretty(x = c(min(txtbl[[1]]), max(txtbl[[2]])))
          dirat[1] = min(txtbl[[1]]) #Avoid drawing arrows outside gene length
          dirat[length(dirat)] = max(txtbl[[2]])
          points(x = dirat, y = rep(tx_id-0.45, length(dirat)), pch = "<", col = exon_col)
        }
      }
    }
  }
  
  #Draw scale
  if(show_axis){
    par(mar = c(0, left_mar, 0, 1))  
  }else{
    par(mar = c(0, left_mar, 0, 1))
  }
  lab_at = pretty(c(start, end))
  lab_at_lab = ifelse(test = lab_at > 1e6, yes = paste0(lab_at/1e6, "M"), no = ifelse(lab_at > 100000, yes = paste0(lab_at/1e5, "K"), no = lab_at))
  plot(NA, xlim = c(start, end), ylim = c(0, 1), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  rect(xleft = start, ybottom = 0.5, xright = end, ytop = 0.5, lty = 2, xpd = TRUE)
  rect(xleft = lab_at, ybottom = 0.45, xright = lab_at, ytop = 0.5, xpd = TRUE)
  text(x = lab_at, y = 0.2, labels = lab_at_lab, xpd = FALSE)
  #axis(side = 1, at = lab_at, lty = 2, line = -3)
  text(x = end, y = 0.9, labels = paste0(chr, ":", start, "-", end), adj = 1, xpd = TRUE)
  
  #Draw ideogram
  if(is.list(cyto)){
    if(show_ideogram){
      par(mar = c(0.2, 1, 0, 1))
      plot(NA, xlim = c(0, max(cyto$end)), ylim = c(0, 1), axes = FALSE, frame.plot = FALSE, xlab = NA, ylab = NA)
      rect(xleft = cyto$start, ybottom = 0.1, xright = cyto$end, ytop = 0.6, col = cyto$color, border = "#34495e")
      rect(xleft = start, ybottom = 0, xright = end, ytop = 0.7, col = "#d35400", lwd = 2, border = "#d35400")
      text(x = 0, y = 0.5, labels = chr, adj = 1.2, font = 2, xpd = TRUE)
    }
  }
  
}

# profileplot is an ultra-fast, simple, and minimal dependency R script to generate profile-plots from bigWig files
#' Generate bigWig signal matrix for given genomic regions or ucsc refseq transcripts 
#' @param colData from \code{read_coldata}
#' @param bed bed file or a data.frame with first 3 column containing chromosome, star, end positions. 
#' @param binSize bin size to extract signal. Default 50 (bps). Should be >1
#' @param startFrom Default "start". For bed files this can be "start", "center" or "end". For `ucsc_assembly` this can only be "start" or "end"
#' @param up extend upstream by this many bps from `startFrom`. Default 2500
#' @param down extend downstream by this many bps from `startFrom`. Default 2500
#' @param ucsc_assembly If `bed` file not provided, setting `ucsc_assembly` to TRUE will fetch transcripts from UCSC genome browser.
#' @param pc_genes Use only protein coding genes when `ucsc_assembly` is used. Default TRUE
#' @param nthreads Default 4
#' @seealso \code{\link{profile_summarize}} \code{\link{profile_plot}} \code{\link{profile_heatmap}}
#' @export

profile_extract = function(colData = NULL, bed = NULL, ucsc_assembly = TRUE, startFrom = "start", binSize = 50,
                           up = 2500, down = 2500, pc_genes = TRUE, nthreads = 4){
  .check_windows()
  .check_bwtool(warn = FALSE)
  .check_dt()
  
  if(is.null(colData)){
    stop("Missing colData. Use read_coldata() to generate one.")
  }
  
  bigWigs = colData$bw_files
  custom_names = colData$bw_sample_names
  
  op_dir = tempdir() #For now
  
  if(is.null(bed)){
    if(ucsc_assembly){
      ucsc_assembly = attr(colData, "refbuild")
      message("No bed file was given. Defaulting to ucsc refseq..")
      startFrom = match.arg(arg = startFrom, choices = c("start", "end"))
      bed = .make_genome_bed(refBuild = ucsc_assembly, up = as.numeric(up), down = as.numeric(down), tss = startFrom, op_dir = op_dir, pc_genes = pc_genes, for_profile = TRUE)
      bed_annot = bed[[2]]
      bed = bed[[1]]
    }else{
      stop("Please provide either a BED file or set ucsc_assembly to TRUE")
    }
  }else{
    startFrom = match.arg(arg = startFrom, choices = c("start", "end", "center"))
    bed = .make_bed(bed = bed, op_dir = op_dir, up = as.numeric(up), down = as.numeric(up), tss = startFrom, for_profile = TRUE)
    bed_annot = NA
  }
  
  message("Extracting signals..")
  mats = parallel::mclapply(bigWigs, function(x){
    .bwt_mats(bw = x, binSize = binSize, bed = bed, size = paste0(up, ":", down), startFrom = startFrom, op_dir = op_dir)
  }, mc.cores = nthreads)
  
  mats = as.character(unlist(x = mats))
  sig_list = lapply(mats, data.table::fread)
  
  if(!is.null(custom_names)){
    names(sig_list) = custom_names
  }else{
    names(sig_list) = gsub(pattern = "*\\.matrix$", replacement = "", x = basename(path = mats))
  }
  
  attr(sig_list, "args") = c(up, down)
  
  #Remove intermediate files
  lapply(mats, function(x) system(command = paste0("rm ", x), intern = TRUE))
  
  list(data = sig_list, colData = colData)
}


#' Summarize data for profile plots
#' @param sig_list Output generated from `profile_extract`
#' @param stat Default `mean`. Can be `mean`, `median`
#' @param condition column name with conditions in `colData`. If provided summarizes signals from samples belonging to same group or condition
#' @seealso \code{\link{profile_extract}} \code{\link{profile_plot}} \code{\link{profile_heatmap}}
#' @export
profile_summarize = function(sig_list = NULL, stat = "mean", condition = NULL){
  
  if(is.null(sig_list)){
    stop("Missing input! Use profile_extract() to generate one.")
  }
  
  stat = match.arg(arg = stat, choices = c("mean", 'median'))
  
  colData = sig_list$colData
  collapse_replicates = FALSE
  if(!is.null(condition)){
    if(!condition %in% colnames(colData)){
      stop(condition , " not found in colData!\nAvailable columns: ", print(paste(colnames(colData), collapse = " ")))
    }
    collapse_by_idx = which(colnames(colData) == condition)
    condition = unlist(colData[,collapse_by_idx, with = FALSE], use.names = FALSE)
    collapse_replicates = TRUE
  }
  
  message("Summarizing..")
  sig_summary = .summarizeMats(mats = sig_list$data, group = condition, collapse_reps = collapse_replicates, summarizeBy = stat)
  attr(sig_summary, "args") = attr(sig_list$data, "args")
  list(data = sig_summary, colData = colData)
}


#' Draw a profile plot
#' @param sig_list Output generated from profile_summarize
#' @param color Manual colors for each bigWig. Default NULL. 
#' @param line_size Default 1
#' @param legend_fs Legend font size. Default 1
#' @param axis_fs Axis font size. Default 1
#' @param xlab x axis label. Default NA
#' @param ylab y axis label. Default NA
#' @export
profile_plot = function(sig_list = NULL, color = NULL, line_size = 1, legend_fs = 1, axis_fs = 1, xlab = NA, ylab = NA){
  
  if(is.null(sig_list)){
    stop("Missing input! Expecting output from profile_summarize()")
  }
  
  args = attr(sig_list$data, "args")
  up = as.numeric(args[1])
  down = as.numeric(args[2])
  
  sig_summary = sig_list$data
  
  if(is.null(color)){
    color = c("#2f4f4f", "#8b4513", "#228b22", "#00008b", "#ff0000", "#ffd700", "#7fff00", "#00ffff", "#ff00ff", "#6495ed", "#ffe4b5", "#ff69b4") #hcl.colors(n = 10, palette = "Dark 2")
    color = color[1:length(sig_summary)]
  }
  
  y_max = max(unlist(lapply(sig_summary, max, na.rm = TRUE)))
  y_min = min(unlist(lapply(sig_summary, min, na.rm = TRUE)))
  ylabs = pretty(c(y_min, y_max), n = 5)
  
  
  x_max = max(unlist(lapply(sig_summary, length)))
  xlabs = c(up, 0, down)
  xticks = xticks = c(0,
                      as.integer(length(sig_summary[[1]])/sum(as.numeric(xlabs[1]), as.numeric(xlabs[3])) * as.numeric(xlabs[1])),
                      length(sig_summary[[1]]))
  
  #line_size = 1
  par(mar = c(4, 4, 2, 1))
  plot(NA, xlim = c(0, x_max), ylim = c(min(ylabs), max(ylabs)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  abline(h = ylabs, v = pretty(xticks), col = "gray90", lty = 2)
  
  lapply(1:length(sig_summary), function(idx){
    points(sig_summary[[idx]], type = 'l', lwd = line_size, col = color[idx])
  })
  axis(side = 1, at = xticks, labels = xlabs, cex.axis = axis_fs)
  axis(side = 2, at = ylabs, las = 2, cex.axis = axis_fs)
  
  legend(x = "topright", legend = names(sig_summary), col = color, bty = "n", lty = 1, lwd = 1.2, cex = legend_fs, xpd = TRUE)
  
  mtext(text = xlab, side = 1, line = 2.5, cex = 1)
  mtext(text = ylab, side = 2, line = 2.5, cex = 1)
  
  invisible(list(mean_signal = sig_summary, color_codes = color, xticks = xticks, xlabs = xlabs))
}

#' Draw a heatmap
#' @details This function takes output from extract_matrix and draws a heatmap
#' @param mat_list Input matrix list generated by \code{\link{profile_extract}}
#' @param sortBy Sort matrix by.., Can be mean, median. Default mean.
#' @param col_pal Color palette to use. Default Blues. Use hcl.pals(type = "sequential") to see available palettes
#' @param revpal Reverse color palette? Default FALSE.
#' @param sample_names Manually specify sample names.
#' @param title_size size of title. Default 0.8
#' @param top_profile Boolean. Whether to draw top profile plot.
#' @param top_profile_h Default 2.
#' @param zmins Manually specify min scores to include in heatmap
#' @param zmaxs Manually specify max scores to include in heatmap
#' @param scale Whether to row scale the matrix. Default FALSE
#' @param file_name Default NULL. If provided saves plot as a png.
#' @param hm_width Width of the plot. Default 1024
#' @param hm_height Height of the plot. Default 600
#' @param mat_order Default NULL. Sample order in which the heatmaps are drawn.
#' @export
profile_heatmap = function(mat_list, sortBy = "mean", col_pal = "Blues", revpal = FALSE, sample_names = NULL,
                        title_size = 1, top_profile = FALSE, top_profile_h = 2, zmins = NULL, zmaxs = NULL,
                        scale = FALSE, file_name = NULL, hm_width = 1024, hm_height = 600, mat_order = NULL){
  
  
  if(!sortBy %in% c("mean", "median")){
    stop("sortBy can only be mean, median")
  }
  
  col_pal = match.arg(arg = col_pal, choices = hcl.pals(type = "sequential"))
  hmcols = rev(hcl.colors(n = 255, palette = col_pal))
  if(revpal){
    hmcols = rev(hmcols)
  }
  
  cdata = mat_list$colData
  size = attr(mat_list$data, "args")
  mat_list = mat_list$data
  
  if(!is.null(mat_order)){
    if(length(mat_order) != length(mat_list)){
      stop("Length of mat_order should be the same as number of bw files! [", length(mat_list), "]")
    }
    mat_list = mat_list[mat_order]
  }
  
  mat_list = .order_matrix(mats = mat_list, sortBy = sortBy)
  
  if(top_profile){
    profile_dat = .summarizeMats(mats = mat_list, summarizeBy = "mean")
    yl = c(min(unlist(x = profile_dat), na.rm = TRUE),
           max(unlist(x = profile_dat), na.rm = TRUE))
    yl = round(x = yl, digits = 2)
    
  }
  
  nsamps = length(mat_list)
  
  if(!is.null(zmins)){
    if(length(zmins) != length(mat_list)){
      warning("zmins are recycled")
      zmins = rep(zmins, length(mat_list))
    }
  }else{
    zmins = unlist(lapply(profile_dat, min, na.rm = TRUE))
  }
  
  if(!is.null(zmaxs)){
    if(length(zmaxs) != length(mat_list)){
      warning("zmaxs are recycled")
      zmaxs = rep(zmaxs, length(mat_list))
    }
  }else{
    zmaxs = unlist(lapply(profile_dat, max, na.rm = TRUE))
  }
  
  if(!is.null(sample_names)){
    if(length(sample_names) != length(mat_list)){
      stop("Number of sample names should be equal to number of samples in the matrix")
    }else{
      sample_names = cdata$bw_sample_names
      names(mat_list) = sample_names
    }
  }
  
  xlabs = c(size[1], 0, size[2])
  xticks = c(0,
             1/(sum(as.integer(xlabs[1]), as.integer(xlabs[3]))/as.integer(xlabs[1])),
             1)
  
  if(!is.null(file_name)){
    png(filename = paste0(file_name, ".png"), res = 100, height = 1024, width = 600)  
  }
  
  if(top_profile){
    lo_mat = matrix(data = 1:(nsamps*2), nrow = 2, ncol = nsamps)
    lo = layout(mat = lo_mat, heights = c(top_profile_h, 9))
  }else{
    lo_mat = matrix(data = 1:nsamps, nrow = 1)
    lo = layout(mat = lo_mat)
  }
  
  for(i in 1:nsamps){
    
    if(top_profile){
      .plot_profile_mini(plot_dat = profile_dat, index = i, ylims = c(zmins[i], zmaxs[i]))
    }
    
    par(mar = c(3,2,1.5,1), las=1, tcl=-.25)
    hm.dat = mat_list[[i]]
    data.table::setDF(x = hm.dat)
    
    if(is.null(zmins)){
      #colMin = apply(hm.dat, 2, min, na.rm = TRUE)
      zmin = round(min(hm.dat, na.rm = TRUE), digits = 2)
      #zmin = 0
    }else{
      zmin = zmins[i]
    }
    
    if(scale){
      hm.dat = scale(x = t(hm.dat), scale = TRUE)
      #hm.dat =  apply(hm.dat, 2, rev)
    }else{
      hm.dat =  t(apply(hm.dat, 2, rev))
    }
    
    if(is.null(zmaxs)){
      zmax = max(rowMeans(hm.dat, na.rm = TRUE)) #mean(apply(hm.dat, 2, max, na.rm = TRUE))
      #colMax = which(x = colMax == max(colMax), arr.ind = TRUE)[1]
      #zmax = round(max(boxplot.stats(unlist(hm.dat[colMax,]))$stats), 2)
    }else{
      zmax = zmaxs[i]
    }
    
    #hmcols = grDevices::colorRampPalette(hmcols)(255)
    hm.dat[hm.dat >= zmax] = zmax
    #return(hm.dat)
    image(hm.dat, axes = FALSE, col = hmcols, useRaster = TRUE,
          zlim = c(zmin, zmax), xlim = c(-0.2, 1))
    
    #Add legend
    image(x = c(-0.1, -0.05), y = seq(0, 1, length.out = length(hmcols)-1),
          z = matrix(data = 1:(length(hmcols)-1), nrow = 1), add = TRUE, col = hmcols)
    axis(side = 2, at = seq(0, 1, length.out = 5),
         labels = NA,
         line = -0.6, font.axis = 2, yaxp  = c(1.1, 1.2, 3), lwd = 1)
    mtext(text = round(seq(zmin, zmax, length.out = 5), digits = 2), side = 2,
          line = -0.1, at = seq(0, 1, length.out = 5), font = 1)
    
    title(main = names(mat_list)[i], cex.main = title_size, font.main = 2)
    axis(side = 1, at = xticks,
         labels = NA, lty = 1, lwd = 1,
         font.axis = 2, cex.axis = 1, line = 0.7, tick = FALSE)
    mtext(text = c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]), side = 1, line = 1, at = xticks, font = 1)
    
    rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, border = "black", lwd = 1)
  }
  
  if(!is.null(file_name)){
    dev.off()
  }
  
}

# bwpcaplot function to perform PCA analysis based on genomic regions of interest or around TSS sites.
#' Extract area under the curve for every peak from from given bigWig files.
#' @param colData bigWig files. Default NULL. Required.
#' @param bed bed file or a data.frame with first 3 column containing chromosome, star, end positions. 
#' @param binSize bin size to extract signal. Default 50 (bps). Should be >1
#' @param ucsc_assembly If `bed` file not provided, setting `ucsc_assembly` to TRUE will fetch transcripts from UCSC genome browser.
#' @param startFrom Default "start". For bed files this can be "start", "center" or "end". For `ucsc_assembly` this can only be "start" or "end"
#' @param pc_genes Use only protein coding genes when using `ucsc_assembly`. Default TRUE
#' @param up extend upstream by this many bps from `startFrom`. Default 2500
#' @param down extend downstream by this many bps from `startFrom`. Default 2500
#' @param nthreads Default 4
#' @export
extract_summary = function(colData, bed = NULL, ucsc_assembly = TRUE, startFrom = "start", binSize = 50,
                           up = 2500, down = 2500, pc_genes = TRUE, nthreads = 4){
  
  if(is.null(colData)){
    stop("Missing colData. Use read_coldata() to generate one.")
  }
  
  .check_windows()
  .check_bwtool()
  .check_dt()
  
  bigWigs = colData$bw_files
  custom_names = colData$bw_sample_names
  
  op_dir = tempdir() #For now
  
  if(is.null(bed)){
    if(ucsc_assembly){
      ucsc_assembly = attr(colData, "refbuild")
      message("No bed file was given. Defaulting to ucsc refseq..")
      startFrom = match.arg(arg = startFrom, choices = c("start", "end"))
      bed = .make_genome_bed(refBuild = ucsc_assembly, up = as.numeric(up), down = as.numeric(down), tss = startFrom, op_dir = op_dir, pc_genes = pc_genes, for_profile = FALSE)
      bed_annot = bed[[2]]
      bed = bed[[1]]
    }else{
      stop("Please provide either a BED file or an ucsc_assembly name")
    }
  }else{
    startFrom = match.arg(arg = startFrom, choices = c("start", "end", "center"))
    bed = .make_bed(bed = bed, op_dir = op_dir, up = as.numeric(up), down = as.numeric(up), tss = startFrom)
    bed_annot = NA
  }
  
  op_dir = tempdir() #For now
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("Extracting summaries..")
  summaries = parallel::mclapply(seq_along(bigWigs), FUN = function(idx){
    bw = bigWigs[idx]
    bn = custom_names[idx]
    cmd = paste("bwtool summary -with-sum -keep-bed -header", bed, bw, paste0(op_dir, bn, ".summary"))
    system(command = cmd, intern = TRUE)
    paste0(op_dir, bn, ".summary")
  }, mc.cores = nthreads)

  merge_anno = FALSE
  
  if(!is.na(bed_annot)){
    bed_annot = data.table::fread(file = bed_annot)
    bed_annot[, id := paste0(V1, ":", V2, "-", V3)]
    colnames(bed_annot) = c("chromosome", "start", "end", "tx", "gene", "id")
    merge_anno = TRUE
  }
  
  summary_list = lapply(summaries, function(x){
    x = data.table::fread(x)
    colnames(x)[1] = 'chromosome'
    x = x[,.(chromosome, start, end, size, sum)]
    #x[,id := paste0(chromosome, ":", start, "-", end)]
    x$sum
  })
  summary_list = data.frame(summary_list)
  colnames(summary_list) = custom_names
  data.table::setDT(summary_list)
  
  if(merge_anno){
    summary_list = cbind(bed_annot, summary_list)
  }else{
    bed = data.table::fread(file = bed)
    colnames(bed) = c("chromosome", "start", "end")
    summary_list = cbind(bed, summary_list)
  }
  
  list(data = summary_list, colData = colData)
}

#' Differential Peak Analysis
#' @details Takes output from \code{extract_summary} and performs differential peak analysis with Limma
#'
#' @param summary_list Output from \code{\link{extract_summary}}
#' @param condition a column name in \code{coldata} containing sample conditions. Default NULL.
#' @param log2 log2 convert data prior to testing. Default TRUE
#' @param num Numerator condition. Default NULL
#' @param den Denominator condition. Default NULL
#' @export
#'
diffpeak = function(summary_list = NULL, condition = NULL, log2 = TRUE, num = NULL, den = NULL){
  
  if (is.null(condition)){
    stop("Define a condition for the differential peak analysis\n")
  }
  
  sum_tbl = as.data.frame(summary_list$data)
  sum_tbl$rid = paste0("rid_", 1:nrow(sum_tbl))
  rownames(sum_tbl) = sum_tbl$rid
  coldata = as.data.frame(summary_list$colData)
  
  if(condition %in% colnames(coldata)){
    condition <- coldata[, condition]
  }else {
    print(coldata)
    stop("The condition argument must be a column name in coldata. See above for list of valid colnames")
  }
  
  exprs = sum_tbl[,coldata$bw_sample_names]
  if(log2){
    exprs = log2(x = exprs + 0.1)
  }
  #sum_tbl = sum_tbl[complete.cases(sum_tbl),, drop = FALSE]
  #sds_idx = order(apply(sum_tbl, 1, sd, na.rm = TRUE), decreasing = TRUE, na.last = TRUE)
  #sum_tbl = sum_tbl[sds_idx,,drop = FALSE]
  
  
  design <- model.matrix(~0 + as.factor(condition))
  colnames(design) <- as.character(levels(as.factor(condition)))
  
  if(is.null(num) & is.null(den)){
    contrast <- vector()
    for (a in 1:length(levels(as.factor(condition)))) {
      for (b in 1:length(levels(as.factor(condition)))) {
        if (a != b) {
          if (a < b) {
            contrast[length(contrast) + 1] <- paste(levels(as.factor(condition))[a],
                                                    levels(as.factor(condition))[b], sep = "-")
          }
        }
      }
    }
    message("Argument num and den are missing. Pefroming diffpeak analysis for below contrast:")
    print(contrast[1])
    num = unlist(data.table::tstrsplit(x = contrast[1], split = "-"))[1]
    den = unlist(data.table::tstrsplit(x = contrast[1], split = "-"))[2]
  }else{
    if(is.null(num) | is.null(den)){
      stop("Num and Den must be provided")
    }else{
      contrast = paste0(num, "-", den)
    }
  }
  
  fit <- limma::lmFit(object = exprs, design = design)
  
  cnt <- paste(colnames(design)[1], colnames(design)[2], sep = "-")
  cMat <- limma::makeContrasts(contrasts = contrast, levels = design)
  fit2 <- limma::contrasts.fit(fit, cMat)
  efit <- limma::eBayes(fit2)
  
  tt = limma::topTable(fit = efit, coef = 1, number = "all")
  tt = merge(sum_tbl, tt, by = "row.names")
  data.table::setDT(x = tt)
  tt$Row.names = tt$rid = NULL
  
  tt = tt[order(P.Value, decreasing = FALSE)]
  attr(tt, "contrast") = contrast
  tt
}

#' Volcano plot for limma results from differential peak Analysis
#' @details Takes output from \code{diffpeak} and draws volcano plot
#' @param res Output from \code{\link{diffpeak}}
#' @param fdr FDR threshold. Default 0.1
#' @param upcol color for up-regulated. Default "#d35400"
#' @param downcol color for up-regulated. Default "#1abc9c"
#' @param alpha Default 0.6
#' @param size Point size. Default 0.8
#' @export
#'
volcano_plot = function(res = NULL, fdr = 0.1, upcol = "#d35400", downcol = "#1abc9c", alpha = 0.6, size = 0.8){
  
  if(is.null(res)){
    stop("Missing input. Expecting output diffpeak()")
  }
  
  contrast = attr(res, "contrast")
  
  ylims = max(-log10(res$P.Value), na.rm = TRUE)
  xlims = range(res$logFC)
  
  res_sig = res[adj.P.Val < fdr]
  res_nonsig = res[!adj.P.Val < fdr]
  
  par(mar = c(3.5, 3.5, 3, 1))
  plot(NA, xlim = xlims, ylim = c(0, ylims), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  grid()
  points(res_nonsig$logFC, y = -log10(res_nonsig$P.Value), pch = 19, col = adjustcolor("gray", alpha.f = alpha), cex = size)
  if(nrow(res_sig)>0){
    points(res_sig$logFC, y = -log10(res_sig$P.Value), pch = 19, col = ifelse(test = res_sig$logFC > 0, adjustcolor(col = upcol, alpha.f = alpha), adjustcolor(col = downcol, alpha.f = alpha)), cex = size)  
  }
  axis(side = 1, at = pretty(xlims))
  axis(side = 2, at = pretty(c(0, ylims)), las = 2)
  mtext(text = "Log Fold Change", side = 1, line = 2, col = "#34495e", font = 2)
  mtext(text = "-log10(P-Value)", side = 2, line = 2, col = "#34495e", font = 2)
  title(main = contrast, col.main = "#34495e")
  
  #Small hack to align legend properly
  leg_pos = which(abs(xlims) == max(abs(xlims)))
  leg_pos = ifelse(test = leg_pos == 1, yes = "bottomleft", no = "bottomright")
  
  legend(x = leg_pos, legend = c(paste0("Down [", nrow(res_sig[logFC < 0]), "]"), paste0("Up [", nrow(res_sig[logFC > 0]), "]")), 
         col = c(downcol, upcol), pch = 19, title = paste0("Diff. peaks [FDR<", fdr, "]"), adj = 0, bty = "n")
}

#' Draw a PCA plot
#' @param summary_list output from extract_summary
#' @param top Top most variable peaks to consider for PCA. Default 1000
#' @param log2 log transform data? Default FALSE. IF TRUE, adds a small positive value and log2 converts.
#' @param xpc Default PC1
#' @param ypc Default PC2
#' @param color_by a column name in \code{coldata} to color by. Default NULL.
#' @param pch_by a column name in \code{coldata} to pch by. Default NULL. 
#' @param color Manual colors for each level in `color_by` Default NULL. 
#' @param show_cree If TRUE draws a cree plot. Default TRUE
#' @param size Point size. Default 1
#' @param lab_size Font size for labels. Default 1
#' @param legend_size Default 1
#' @param legendpos Default topright
#' @param legendpos2 Default bottomright
#' @export
pca_plot = function(summary_list = NULL, top = 1000, log2 = FALSE, xpc = "PC1", ypc = "PC2", color_by = NULL, pch_by = NULL, color = NULL, 
                     show_cree = TRUE, lab_size = 1, size = 1, legend_size = 1, legendpos = "topright", legendpos2 = "bottomright"){
  
  
  if(is.null(summary_list)){
    stop("Missing input! Expecting output from extract_summary()")
  }
  
  sum_tbl = as.data.frame(summary_list$data)
  colData = summary_list$colData
  
  sum_tbl = sum_tbl[,colData$bw_sample_names]
  if(log2){
    sum_tbl = log2(x = sum_tbl + 0.1)
  }
  #sum_tbl = sum_tbl[complete.cases(sum_tbl),, drop = FALSE]
  sds_idx = order(apply(sum_tbl, 1, sd, na.rm = TRUE), decreasing = TRUE, na.last = TRUE)
  sum_tbl = sum_tbl[sds_idx,,drop = FALSE]
  
  
  if(is.null(color)){
    condition_colors = c("#2f4f4f", "#8b4513", "#228b22", "#00008b", "#ff0000", "#ffd700", "#7fff00", "#00ffff", "#ff00ff", "#6495ed", "#ffe4b5", "#ff69b4")
    # condition_colors = c("#A6CEE3FF", "#1F78B4FF", "#B2DF8AFF", "#33A02CFF", "#FB9A99FF", 
    #                      "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF", "#CAB2D6FF", "#6A3D9AFF", 
    #                      "#FFFF99FF", "#9E0142FF", "#D53E4FFF", "#F46D43FF", "#000000FF", 
    #                      "#EE82EEFF", "#4169E1FF", "#7B7060FF", "#535C68FF")
  }else{
    condition_colors = color
  }
  
  condition = color_by
  pch = pch_by
  
  if(!is.null(condition)){
    group_df = cbind(sample = colData$bw_sample_names, condition = colData[,which(colnames(colData) == condition), with = FALSE])
    colnames(group_df) = c("sample", "condition")
    condition_colors = condition_colors[1:nrow(group_df[,.N,condition])]
    names(condition_colors) = group_df[,.N,condition][,condition]
    group_df$color = condition_colors[group_df$condition]
    condition_leg = TRUE
  }else{
    group_df = data.table::data.table(sample = colData$bw_sample_names)
    group_df$color = "black"
    condition_leg = FALSE
  }
  
  
  if(!is.null(pch)){
    pch_vec = 15:19
    pch_df = cbind(sample = colData$bw_sample_names, pch = colData[,which(colnames(colData) == pch), with = FALSE])
    colnames(pch_df) = c("sample", "pch")
    pchs = pch_vec[1:nrow(pch_df[,.N,pch])]
    names(pchs) = pch_df[,.N,pch][,pch]
    pch_df$pch = pchs[pch_df$pch]
    pch_leg = TRUE
  }else{
    pch_df = cbind(sample = colData$bw_sample_names, pch = 19)
    colnames(pch_df) = c("sample", "pch")
    pchs = NA
    pch_leg = FALSE
  }
  
  group_df = merge(group_df, pch_df, by = "sample", all = TRUE)
  
  if(nrow(sum_tbl) < top){
    pca = prcomp(t(sum_tbl))
  }else{
    pca = prcomp(t(sum_tbl[1:top,]))
  }
  
  pca_dat = as.data.frame(pca$x)
  pca_var_explained = pca$sdev^2/sum(pca$sdev^2)
  names(pca_var_explained) = paste0("PC", 1:length(pca_var_explained))
  data.table::setDT(x = pca_dat, keep.rownames = "sample")
  pca_dat = merge(pca_dat, group_df, by = 'sample')
  data.table::setDF(x = pca_dat)
  attr(pca_dat, "percentVar") <- round(pca_var_explained, digits = 3)
  #print(head(pca_dat))
  
  if(show_cree){
    lo = layout(mat = matrix(data = c(1, 2), ncol = 2))
  }
  
  grid_cols = "gray90"
  par(mar = c(3, 4, 2, 1))
  plot(NA, axes = FALSE, xlab = NA, ylab = NA, cex = 1.2, xlim = range(pretty(pca_dat[, xpc])), ylim = range(pretty(pca_dat[, ypc])))
  abline(h = pretty(pca_dat[, xpc]), v = pretty(pca_dat[, ypc]), col = grid_cols, lty = 2, lwd = 0.1)
  abline(h = 0, v = 0, col = grid_cols, lty = 2, lwd = 0.8)
  points(x = pca_dat[, xpc], y = pca_dat[, ypc], col = pca_dat$color, bg = pca_dat$color, pch = as.integer(pca_dat$pch), cex = size)
  axis(side = 1, at = pretty(pca_dat[, xpc]), cex.axis = 0.8)
  axis(side = 2, at = pretty(pca_dat[, ypc]), las = 2, cex.axis = 0.8)
  mtext(text = paste0(xpc, " [", round(pca_var_explained[xpc], digits = 2), "]"), side = 1, line = 2, cex = 0.8)
  mtext(text = paste0(ypc, " [", round(pca_var_explained[ypc], digits = 2), "]"), side = 2, line = 2, cex = 0.8)
  text(x = pca_dat[, xpc], y = pca_dat[, ypc], labels = pca_dat$sample, pos = 3, col = pca_dat$color, xpd = TRUE, cex = lab_size)
  #title(main = NA, sub = paste0("N = ", top, " peaks"), adj = 0, outer = TRUE)
  
  data.table::setDT(x = pca_dat)
  
  if(condition_leg){
    legend(x = legendpos, legend = names(condition_colors), col = condition_colors, 
           pch = 19, cex = legend_size, xpd = TRUE, ncol = 2, bty = "n")
  }
  
  if(pch_leg){
    legend(x = legendpos2, legend = names(pchs), col = "black", 
           pch = as.integer(pchs), cex = legend_size, xpd = TRUE, ncol = 1, bty = "n")
  }
  
  if(show_cree){
    par(mar = c(3, 4, 2, 4))
    b = barplot(height = pca_var_explained, names.arg = NA, col = "#2c3e50", ylim = c(0, 1), las = 2, axes = FALSE)
    points(x = b, y = cumsum(pca_var_explained), type = 'l', lty = 2, lwd = 1.2, xpd = TRUE, col = "#c0392b")
    points(x = b, y = cumsum(pca_var_explained), type = 'p', pch = 19, xpd = TRUE, col = "#c0392b")
    mtext(text = paste0("PC", 1:length(pca_var_explained)), side = 1, at = b, las = 2, line = 0.5, cex = 0.8)
    axis(side = 2, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
    mtext(text = "var. explained", side = 2, line = 2.5)
    axis(side = 4, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
    mtext(text = "cumulative var. explained", side = 4, line = 2.5)
  }
  
  invisible(x = pca_dat)
}


#' Parse peak annotations generated by homer annotatePeaks.pl
#' @details summarizes peak annotations generated with homer annotatePeaks.pl, generates a pie chart of peak distributions.
#'
#' @param anno Raw annotations generated by homer `annotatePeaks.pl`. Can be more than one file.
#' @param sample_names Sample names correspoding to each input file. Default parses from input file.
#' @param legend_font_size font size for legend. Default 1.
#' @param label_size font size for labels. Default 0.8.
#' @export

summarize_homer_annots = function(anno, sample_names = NULL, legend_font_size = 1, label_size = 0.8){
  
  homer = lapply(anno, function(homer.anno){
    homer.anno = data.table::fread(homer.anno)
    colnames(homer.anno)[1] = 'peakid'
    homer.anno$anno = sapply(strsplit(x = as.character(homer.anno$Annotation), split = ' (', fixed = T), '[[', 1)
    
    homer.anno = homer.anno[,.(Chr, Start, End, Strand, anno, `Gene Name`, `Gene Type`,`Distance to TSS`, `Nearest PromoterID`, `Nearest Ensembl`)]
    homer.anno = homer.anno[order(anno, `Gene Name`)]
    colnames(homer.anno) = c('Chr', 'Start', 'End', 'Strand', 'Annotation', 'Hugo_Symbol',
                             'Biotype', 'Distance_to_TSS', 'Nearest_PromoterID', 'Nearest_Ens')
    homer.anno$Annotation = gsub(pattern = "3' UTR", replacement = '3pUTR', x = homer.anno$Annotation)
    homer.anno$Annotation = gsub(pattern = "5' UTR", replacement = '5pUTR', x = homer.anno$Annotation)
    homer.anno
  })
  
  if(is.null(sample_names)){
    names(homer) = unlist(data.table::tstrsplit(x = basename(anno), split = "\\.", keep = 1))
  }else{
    names(homer) = sample_names
  }
  
  homer.anno.stats = lapply(homer, function(h){
    homer.anno.stats = h[,.N,Annotation][,fract := N/sum(N)][order(N, decreasing = TRUE)]
    homer.anno.stats[,leg := paste0(Annotation, " [", N, "]")]
    homer.anno.stats
  })
  npeaks = unlist(lapply(homer, nrow))
  
  homer.anno.stats = data.table::rbindlist(l = homer.anno.stats, idcol = "Sample")
  homer.anno.stats = data.table::dcast(data = homer.anno.stats, Annotation ~ Sample, value.var = "fract", fill = 0)
  data.table::setDF(x = homer.anno.stats, rownames = homer.anno.stats$Annotation)
  homer.anno.stats =homer.anno.stats[,-1]
  
  pie.cols = c('3pUTR' = '#E7298A', '5pUTR' = '#D95F02', 'Intergenic' = '#BEBADA',
               'TTS' = '#FB8072', 'exon' = '#80B1D3', 'intron' = '#FDB462',
               'non-coding' = '#FFFFB3',
               'NA' = 'gray70', 'promoter-TSS' = '#1B9E77')
  
  homer.anno.stats = homer.anno.stats[names(pie.cols)[names(pie.cols) %in% rownames(homer.anno.stats)],,]
  
  pie.col = pie.cols[rownames(homer.anno.stats)]
  
  #layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(4, 1.25))
  par(mar = c(2, 4, 5, 3))
  b = barplot(height = as.matrix(homer.anno.stats), col = pie.col, horiz = TRUE,
              las = 2, axes = FALSE, names.arg = rep(NA, ncol(homer.anno.stats)), border = pie.cols)
  axis(side = 1, at = seq(0, 1, 0.25), font = 2, lwd = 2)
  mtext(text = npeaks[colnames(x = homer.anno.stats)], side = 4, at = b, las = 2, font = 4, adj = -0.2)
  mtext(text = colnames(homer.anno.stats), side = 2, at = b, las = 2, adj = 1, font = 2)
  
  #plot.new()
  #par(mar = c(1, 0, 1, 1))
  .add_legend("topleft", legend = names(pie.col), col = pie.col,  bty = "n", border=NA,
             xpd = TRUE, text.font = 2, pch = 15, xjust = 0, yjust = 0,
             cex = 0.8, y.intersp = 1.5, x.intersp = 1,
             pt.cex = 1.2 * 0.8, ncol = 3)
  
  homer
}

# if(plot){
#   homer.anno.stats = homer.anno[,.N,Annotation][,fract := round(N/sum(N), digits = 2)][order(N, decreasing = TRUE)]
#   homer.anno.stats[,leg := paste0(Annotation, " [", N, "]")]
#
#   pie.cols = c('3pUTR' = '#E7298A', '5pUTR' = '#D95F02', 'Intergenic' = '#BEBADA',
#                'TTS' = '#FB8072', 'exon' = '#80B1D3', 'intron' = '#FDB462', 'non-coding' = '#FFFFB3',
#                'NA' = 'gray70', 'promoter-TSS' = '#1B9E77')
#
#   par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=4,xpd=NA, mar=c(0,0,1,3))
#   pie(x = homer.anno.stats$fract, col = pie.cols[homer.anno.stats$Annotation], labels = homer.anno.stats$fract,
#       border = "white", radius = 0.8, init.angle = 0, font = 4, cex = label_size)
#   symbols(0,0,circles=.3, inches=FALSE, col="white", bg="white", lty=0, add=TRUE)
#   add_legend(x = "topright", bty = "n", legend = homer.anno.stats$leg,
#              col = pie.cols[homer.anno.stats$Annotation],
#              cex = legend_font_size, pch = 15, text.font = 4)
# }

#------------------------------------------------------------------------------------------------------------------------------------
#                                                   Undocumented Accessory functions
#------------------------------------------------------------------------------------------------------------------------------------

# Add legends outside the margin
.add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

.make_layout = function(ntracks, ntracks_h = 3, cytoband = TRUE, cytoband_h = 1, genemodel = TRUE, genemodel_h = 1, chrHMM = TRUE, chrHMM_h = 1, loci = TRUE, loci_h = 2, scale_track_height = 1, lord = NULL){
  
  #track heights (peaks, bigwig tracks, chromHMM, gene model, cytoband, scale)
  lo_h_ord =  list("p" = loci_h, "b" = rep(ntracks_h, ntracks), "h" = chrHMM_h, "g" = genemodel_h, "c" = cytoband_h, "s" = scale_track_height)
  mat_ord = c("p" = 1, "b" = 2, "h" = 3, "g" = 4, "c" = 5, "s" = 6)
  
  case = NULL
  #print(paste(cytoband, genemodel, chrHMM, loci, sep = ", "))
  
  if(cytoband & genemodel & chrHMM & loci){
    #print("all")
    case = 1
    lo_h_ord = lo_h_ord[c("p", "b", "h", "g", "s", "c")] #c(loci_h, rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height, cytoband_h)
  }else if(cytoband & genemodel & chrHMM == FALSE & loci){
    #print("no hmm")
    case = 2
    lo_h_ord = lo_h_ord[c("p", "b", "g", "s", "c")] #c(loci_h, rep(ntracks_h, ntracks), genemodel_h, scale_track_height, cytoband_h)
  }else if(cytoband & genemodel & chrHMM == FALSE & loci == FALSE){
    #print("no hmm and no loci")
    case = 3
    lo_h_ord = lo_h_ord[c("b", "g", "s", "c")] #c(rep(ntracks_h, ntracks), genemodel_h, scale_track_height, cytoband_h)
  }else if(cytoband & genemodel == FALSE & chrHMM == FALSE & loci == FALSE){
    #print("no hmm and no loci and no genemodel")
    case = 4
    lo_h_ord = lo_h_ord[c("b", "s", "c")] #c(rep(ntracks_h, ntracks), scale_track_height, cytoband_h)
  }else if(cytoband == FALSE & genemodel == FALSE & chrHMM == FALSE & loci == FALSE){
    #print("no hmm and no loci and no genemodel and no cytoband")
    case = 5
    lo_h_ord = lo_h_ord[c("b", "s")] #c(rep(ntracks_h, ntracks), scale_track_height)
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == TRUE & loci == FALSE){
    #print("no loci and no cytoband")
    case = 6
    lo_h_ord = lo_h_ord[c("b", "h", "g", "s")] #c(rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height)
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == TRUE & loci == TRUE){
    #print("no cytoband")
    case = 7
    lo_h_ord = lo_h_ord[c("p", "b", "h", "g", "s")] #c(loci_h, rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height)
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == FALSE & loci == TRUE){
    #print("no cytoband no chrHMM")
    case = 8
    lo_h_ord = lo_h_ord[c("p", "b", "g", "s")] #c(loci_h, rep(ntracks_h, ntracks), genemodel_h, scale_track_height)
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == FALSE & loci == FALSE){
    #print("no cytoband no chrHMM no loci")
    case = 9
    lo_h_ord = lo_h_ord[c("b", "g", "s")] #c(rep(ntracks_h, ntracks), genemodel_h, scale_track_height)
  }else if(cytoband == TRUE & genemodel == TRUE & chrHMM == TRUE & loci == FALSE){
    #print("no loci")
    case = 10
    lo_h_ord = lo_h_ord[c("b", "h", "g", "s", "c")] #c(rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height, cytoband_h)
  }else if(cytoband == TRUE & genemodel == FALSE & chrHMM == FALSE & loci == TRUE){
    #print("no genemodel no chrhmm")
    case = 11
    lo_h_ord = lo_h_ord[c("p", "b", "s", "c")] #c(loci_h, rep(ntracks_h, ntracks), scale_track_height, cytoband_h)
  }else{
    #print("Something is wrong!")
  }
  
  lord = c(lord, "s")
  lord = c(intersect(lord, names(lo_h_ord)), setdiff(names(lo_h_ord), lord))
  lo_heights = unlist(lo_h_ord[lord], use.names = FALSE)
  lo = lo_h_ord[lord]
  
  #Order in which plots are drawn by default
  plot_ord = data.frame(row.names = c("p", "b", "h", "g", "s", "c"), name = c("p", "b", "h", "g", "s", "c"), ord = 1:6)
  plot_ord = plot_ord[names(lo),,drop = FALSE]
  plot_ord$ord_req = 1:nrow(plot_ord) #Required order
  plot_ord = plot_ord[order(plot_ord$ord),]
  
  ### Re-organize the layout by user specification
  dt = data.table::data.table(name = plot_ord$name, ord = plot_ord$ord, ord_req = plot_ord$ord_req)
  dt$n_tracks = ifelse(test = dt$name == 'b', yes = ntracks, no = 1)
  dt_s = split(dt, dt$n_tracks)
  if(length(dt_s) > 1){
    dt_mult = dt_s[[2]]
    dt_sing = dt_s[[1]]
    dt_sing[,n_tracks := 1]
    dt_sing$ord_req2 = dt_sing$ord_req
    dt_mult = data.table::data.table(name = rep(dt_mult$name, dt_mult$n_tracks), ord = dt_mult$ord, ord_req = dt_mult$ord_req, ord_req2 = seq(dt_mult$ord_req, dt_mult$ord_req + dt_mult$n_tracks -1 , 1), n_tracks = dt_mult$n_tracks)
    
    dt = data.table::rbindlist(l = list(dt_sing, dt_mult), use.names = TRUE, fill = TRUE)
    dt = dt[order(ord_req)]
    
    dt = split(dt, dt$ord) |> data.table::rbindlist()
    dt[, ord_req4 := 1:nrow(dt)]
    dt = dt[order(ord_req)]
    data = dt$ord_req4
  }else{
    dt = data.table::data.table(name = rep(dt$name, dt$n_tracks), ord = dt$ord, ord_req = dt$ord_req, ord_req2 = seq(dt$ord_req, dt$ord_req + dt$n_tracks -1 , 1), n_tracks = dt$n_tracks)
    dt = split(dt, dt$ord) |> data.table::rbindlist()
    dt[, ord_req4 := 1:nrow(dt)]
    dt = dt[order(ord_req)]
    data = dt$ord_req
  }

  if(case == 1){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 2){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 3){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 4){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 5){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 6){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 7){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 8){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 9){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 10){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }else if(case == 11){
    lo = layout(mat = matrix(data = data), heights = lo_heights)
  }
  
  lo_h_ord[lord]
}

.gen_windows = function(chr = NA, start, end, window_size = 50, op_dir = getwd()){
  #chr = "chr19"; start = 15348301; end = 15391262; window_size = 50; op_dir = getwd()
  message(paste0("Generating windows ", "[", window_size, " bp window size]"))
  
  window_dat = data.table::data.table()
  #temp = start;
  while(start <= end){
    window_dat = data.table::rbindlist(l = list(window_dat, data.table::data.table(start, end = start + window_size)), fill = TRUE)
    start = start + window_size
  }
  window_dat$chr = chr
  window_dat = window_dat[, .(chr, start, end)]
  
  op_dir = paste0(op_dir, "/")
  
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  temp_op_bed = tempfile(pattern = "trackr", tmpdir = op_dir, fileext = ".bed")
  data.table::fwrite(x = window_dat, file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  temp_op_bed
}


.get_summaries = function(bedSimple, bigWigs, op_dir = getwd(), nthreads = 1){
  #bedSimple = temp_op_bed; bigWigs = list.files(path = "./", pattern = "bw"); op_dir = getwd(); nthreads = 1
  op_dir = paste0(op_dir, "/")
  
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  message(paste0("Extracting signals"))
  
  summaries = parallel::mclapply(bigWigs, FUN = function(bw){
    bn = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "", x = basename(bw))
    message(paste0("    Processing ", bn, " .."))
    cmd = paste("bwtool summary -with-sum -keep-bed -header", bedSimple, bw, paste0(op_dir, bn, ".summary"))
    system(command = cmd, intern = TRUE)
    paste0(op_dir, bn, ".summary")
  }, mc.cores = nthreads)
  
  summary_list = lapply(summaries, function(x){
    x = data.table::fread(x)
    colnames(x)[1] = 'chromosome'
    x = x[,.(chromosome, start, end, size, max)]
    if(all(is.na(x[,max]))){
      message("No signal! Possible cause: chromosome name mismatch between bigWigs and queried loci.")
      x[, max := 0]
    }
    x
  })
  
  
  #Remove intermediate files
  lapply(summaries, function(x) system(command = paste0("rm ", x), intern = TRUE))
  system(command = paste0("rm ", bedSimple), intern = TRUE)
  
  names(summary_list) = gsub(pattern = "*\\.summary$", replacement = "", x = basename(path = unlist(summaries)))
  summary_list
}


.extract_geneModel = function(ucsc_tbl = NULL, chr = NULL, start = NULL, end = NULL, txname = txname, genename = genename){
  
  message("Parsing UCSC file..")
  if(is(object = ucsc_tbl, class2 = "data.frame")){
    ucsc = data.table::as.data.table(ucsc_tbl)
  }else if(file.exists(ucsc_tbl)){
    ucsc = data.table::fread(file = ucsc_tbl)  
  }else{
    stop("Gene model must be a file or a data.frame")
  }
  
  query = data.table::data.table(chr, start, end)
  query[,chr := as.character(chr)]
  query[,start := as.numeric(as.character(start))]
  query[,end := as.numeric(as.character(end))]
  data.table::setkey(x = query, chr, start, end)
  
  colnames(ucsc)[c(3, 5, 6)] = c("chr", "start", "end")
  data.table::setkey(x = ucsc, chr, start, end)
  
  gene_models = data.table::foverlaps(x = query, y = ucsc, type = "any", nomatch = NULL)
  
  if(nrow(gene_models) == 0){
    message("No features found within the requested loci! If you are not sure why..\n    1.Make sure there are no discripancies in chromosome names i.e, chr prefixes\n")
    return(NULL)  
  }else{
    if(!is.null(txname)){
      gene_models = gene_models[name %in% txname]
    }
    
    if(!is.null(genename)){
      gene_models = gene_models[name2 %in% genename]
    }
    
    exon_tbls = lapply(seq_along(along.with = 1:nrow(gene_models)), function(idx){
      exon_start = as.numeric(unlist(data.table::tstrsplit(x = gene_models[idx, exonStarts], split = ",")))
      exon_end = as.numeric(unlist(data.table::tstrsplit(x = gene_models[idx, exonEnds], split = ",")))
      exon_tbl = data.frame(start = exon_start, end = exon_end)
      attributes(exon_tbl) = list(start = gene_models[idx, start], end = gene_models[idx, end], strand = gene_models[idx, strand], tx = gene_models[idx, name], gene = gene_models[idx, name2])
      exon_tbl
    })
    
    return(exon_tbls)
  }
}

.extract_cytoband = function(chr = NULL, refBuild = "hg19", tblName = "cytoBand"){
  
  if(!grepl(pattern = "^chr", x = chr)){
    message("Adding chr prefix to target chromosome for UCSC query..")
    chr = paste0("chr", chr)
  }
  
  message(paste0("Extracting cytobands from UCSC:\n", "    chromosome: ", chr, "\n", "    build: ", refBuild))
  mydb = RMySQL::dbConnect(RMySQL::MySQL(), user = "genome", db = refBuild, host = "genome-mysql.soe.ucsc.edu")
  cmd = paste0("select chrom, chromStart, chromEnd, name, gieStain from ", tblName, " WHERE chrom = '", chr, "'")
  cyto = DBI::dbGetQuery(mydb, cmd)
  
  colnames(cyto) = c("chr", "start", "end", "band", "stain")
  data.table::setDT(cyto)
  data.table::setkey(x = cyto, chr, start, end)
  
  #Color codes from https://github.com/jianhong/trackViewer (Thank you..)
  ### gieStain #############################
  # #FFFFFF - gneg    - Giemsa negative bands
  # #999999 - gpos25  - Giemsa positive bands
  # #666666 - gpos50  - Giemsa positive bands
  # #333333 - gpos75  - Giemsa positive bands
  # #000000 - gpos100 - Giemsa positive bands
  # #660033 - acen    - centromeric regions
  # #660099 - gvar    - variable length heterochromatic regions
  # #6600cc - stalk   - tightly constricted regions on the short arms of
  #                     the acrocentric chromosomes
  colorSheme = c(
    "gneg"    = "#FFFFFF",
    "acen"    = "#660033",
    "gvar"    = "#660099",
    "stalk"   = "#6600CC"
  )
  gposCols <- sapply(1:100, function(i){
    i <- as.hexmode(round(256-i*2.56, digits = 0))
    i <- toupper(as.character(i))
    if(nchar(i)==1) i <- paste0("0", i)
    return(paste0("#", i, i, i))
  })
  names(gposCols) <- paste0("gpos", 1:100)
  colorSheme <- c(gposCols, colorSheme)
  cyto$color = colorSheme[cyto[,stain]]
  cyto
  
  cyto
}


.load_chromHMM = function(chr, start, end, ucsc){
  
  if(nrow(ucsc) == 0){
    message("No features found within the requested loci!")
    return(NULL)
  }
  colnames(ucsc)[1:4] = c("chr", "start", "end", "name")
  if(!grepl(pattern = "^chr", x = chr)){
    ucsc[, chr := gsub(pattern = "^chr", replacement = "", x = chr)]
  }
  data.table::setkey(x = ucsc, chr, start, end)
  
  query = data.table::data.table(chr = chr, start = start, end = end)
  data.table::setkey(x = query, chr, start, end)
  
  data.table::foverlaps(x = query, y = ucsc, type = "any", nomatch = NULL)[,.(chr, start, end, name)]
}

.extract_chromHmm_ucsc = function(chr, start, end, refBuild = "hg38", tbl){
  
  if(!grepl(pattern = "^chr", x = chr)){
    message("Adding chr prefix to target chromosome for UCSC query..")
    tar_chr = paste0("chr", chr)
  }else{
    tar_chr = chr
  }
  
  ucsc_tbls = .get_ucsc_hmm_tbls()
  tbl = match.arg(arg = tbl, choices = ucsc_tbls$TableName)
  
  message(paste0("Extracting chromHMM from UCSC:\n", "    chromosome: ", tar_chr, "\n", "    build: ", refBuild))
  mydb = RMySQL::dbConnect(RMySQL::MySQL(), user = "genome", db = refBuild, host = "genome-mysql.soe.ucsc.edu")
  cmd = paste0("select chrom, chromStart, chromEnd, name from ", tbl, " WHERE chrom = '", chr, "'")
  ucsc = DBI::dbGetQuery(mydb, cmd)
  
  if(nrow(ucsc) == 0){
    message("No features found within the requested loci!")
    return(NULL)
  }
  data.table::setDT(x = ucsc)
  colnames(ucsc) = c("chr", "start", "end", "name")
  if(!grepl(pattern = "^chr", x = chr)){
    ucsc[, chr := gsub(pattern = "^chr", replacement = "", x = chr)]
  }
  data.table::setkey(x = ucsc, chr, start, end)
  
  query = data.table::data.table(chr = chr, start = start, end = end)
  data.table::setkey(x = query, chr, start, end)
  
  data.table::foverlaps(x = query, y = ucsc, type = "any", nomatch = NULL)[,.(chr, start, end, name)]
}

.extract_geneModel_ucsc_bySymbol = function(genesymbol, refBuild){
  
  message(paste0("Extracting gene models from UCSC:\n", "    Gene: ", genesymbol, "\n", "    build: ", refBuild))
  mydb = RMySQL::dbConnect(RMySQL::MySQL(), user = "genome", db = refBuild, host = "genome-mysql.soe.ucsc.edu")
  cmd = paste0("select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds from refGene WHERE name2 = '", genesymbol, "'")
  ucsc = DBI::dbGetQuery(mydb, cmd)
  data.table::setDT(ucsc)
  
  if(nrow(ucsc) == 0){
    message("No features found within the requested loci!")
    return(NULL)
  }
  
  colnames(ucsc) = c("chr", "start", "end", "strand", "name", "name2", "exonStarts", "exonEnds")
  data.table::setkey(x = ucsc, chr, start, end)
  ucsc
}

.extract_geneModel_ucsc = function(chr, start = NULL, end = NULL, refBuild = "hg19", txname = NULL, genename = NULL){
  
  if(!grepl(pattern = "^chr", x = chr)){
    message("Adding chr prefix to target chromosome for UCSC query..")
    tar_chr = paste0("chr", chr)
  }else{
    tar_chr = chr
  }
  
  message(paste0("Extracting gene models from UCSC:\n", "    chromosome: ", tar_chr, "\n", "    build: ", refBuild))
  mydb = RMySQL::dbConnect(RMySQL::MySQL(), user = "genome", db = refBuild, host = "genome-mysql.soe.ucsc.edu")
  cmd = paste0("select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds from refGene WHERE chrom = '", tar_chr, "'")
  ucsc = DBI::dbGetQuery(mydb, cmd)
  data.table::setDT(ucsc)
  
  if(nrow(ucsc) == 0){
    message("No features found within the requested loci!")
    return(NULL)
  }
  colnames(ucsc) = c("chr", "start", "end", "strand", "name", "name2", "exonStarts", "exonEnds")
  if(!grepl(pattern = "^chr", x = chr)){
    ucsc[, chr := gsub(pattern = "^chr", replacement = "", x = chr)]
  }
  data.table::setkey(x = ucsc, chr, start, end)
  
  query = data.table::data.table(chr = chr, start = start, end = end)
  data.table::setkey(x = query, chr, start, end)
  
  gene_models = data.table::foverlaps(x = query, y = ucsc, type = "any", nomatch = NULL)
  
  if(nrow(gene_models) == 0){
    message("No features found within the requested loci!")
    return(NULL) 
  }else{
    return(gene_models)
  }
}

.make_exon_tbl = function(gene_models, txname = NULL, genename = NULL){
  if(!is.null(txname)){
    gene_models = gene_models[name %in% txname]
  }
  
  if(nrow(gene_models) == 0){
    message("    Requested transcript ", txname, " does not exist within the queried region!\n    Skipping gene track plotting..")
    return(NULL)
  }
  
  if(!is.null(genename)){
    gene_models = gene_models[name2 %in% genename]
  }
  
  if(nrow(gene_models) == 0){
    message("    Requested gene ", genename, " does not exist within the queried region!\n    Skipping gene track plotting..")
    return(NULL)
  }
  
  exon_tbls = lapply(seq_along(along.with = 1:nrow(gene_models)), function(idx){
    exon_start = as.numeric(unlist(data.table::tstrsplit(x = gene_models[idx, exonStarts], split = ",")))
    exon_end = as.numeric(unlist(data.table::tstrsplit(x = gene_models[idx, exonEnds], split = ",")))
    exon_tbl = data.frame(start = exon_start, end = exon_end)
    attributes(exon_tbl) = list(start = gene_models[idx, start], end = gene_models[idx, end], strand = gene_models[idx, strand], tx = gene_models[idx, name], gene = gene_models[idx, name2])
    exon_tbl
  })
  
  return(exon_tbls)
}

.get_ucsc_hmm_states_cols = function(){
  states = c("red", "red4", "purple", "orange", "orange", "yellow", "yellow", 
             "blue", "darkgreen", "darkgreen", "lightgreen", "gray", "gray90", 
             "gray90", "gray90")
  names(states) = 1:15
  states
}

.plot_ucsc_chrHmm = function(d, start = NULL, end = NULL, hmm_cols = NULL){
  #hmm_cols = .get_ucsc_hmm_states_cols()
  
  if(is.null(start)){
    start = min(data.table::rbindlist(l = d, use.names = TRUE, fill = TRUE)[,start], na.rm = TRUE)  
  }
  
  if(is.null(end)){
    end = max(data.table::rbindlist(l = d, use.names = TRUE, fill = TRUE)[,end], na.rm = TRUE)
  }
  
  plot(NA, ylim = c(0, length(d)), xlim = c(start, end), axes = FALSE, xlab = NA, ylab = NA)
  lapply(seq_along(d), function(i){
    di = d[[i]]
    di$state = unlist(data.table::tstrsplit(x = di$name, split = "_", keep = 1))
    rect(xleft = di$start, ybottom = i-0.9, xright = di$end, ytop = i-0.1, col = hmm_cols[di$state], border = NA)
    di_name = gsub(pattern = "wgEncodeBroadHmm|HMM", replacement = "", x = names(d)[i])
    text(x = start, y = i - 0.5, labels = di_name, adj = 1.2, xpd = TRUE)
    #mtext(text = di_name, side = 2, line = 1, outer = TRUE, cex = 1)
  })
  
}

.collapse_tx = function(exon_tbls){
  message("Collapsing transcripts..")
  tx_tbl = lapply(exon_tbls, function(x){
    xdt = data.table::data.table(start = x[[1]], end = x[[2]])
    xdt$tx = attr(x = x, which = "tx")
    xdt$gene = attr(x = x, which = "gene")
    xdt$strand = attr(x = x, which = "strand")
    xdt$tx_start = attr(x = x, which = "start")
    xdt$tx_end = attr(x = x, which = "end")
    xdt
  })
  tx_tbl = data.table::rbindlist(l = tx_tbl)
  tx_tbl[,id := paste0(start, ":", end)]
  tx_tbl = tx_tbl[!duplicated(id)]
  
  exon_tbls = lapply(split(tx_tbl, as.factor(as.character(tx_tbl$gene))), function(x){
    x = x[order(start)]
    exon_start = as.numeric(x[,start])
    exon_end = as.numeric(x[,end])
    gene_tbl = data.frame(start = exon_start, end = exon_end)
    attributes(gene_tbl) = list(start = min(x[, tx_start]), end = max(x[, tx_end]), strand = unique(x[, strand]), tx = NA, gene = unique(x[, gene]))
    gene_tbl
  })
  
  exon_tbls
}

.parse_gtf = function(gtf = NULL, chr, start = NULL, end = NULL, refBuild = "hg19", txname = NULL, genename = NULL){
  message("Parsing gtf file..")
  if(is(object = gtf, class2 = "data.frame")){
    gtf = data.table::as.data.table(gtf)
  }else if(file.exists(gtf)){
    gtf = data.table::fread(file = gtf)  
  }else{
    stop("Gene model must be a file or a data.frame")
  }
  
  colnames(gtf) = c("chr", "source", "feature", "start", "end", "ph", "strand", "ph2", "info")
  gtf[,chr := as.character(chr)]
  gtf[,start :=as.numeric(as.character(start))]
  gtf[,end := as.numeric(as.character(end))]
  data.table::setkey(x = gtf, chr, start, end)
  
  if(!is.null(genename)){
    gene_models = gtf[info %like% genename]
    if(nrow(gene_models) == 0){
      message("No features found for the gene ", genename)
      return(NULL)  
    }
  }else{
    query = data.table::data.table(chr, start, end)
    data.table::setkey(x = query, chr, start, end)
    gene_models = data.table::foverlaps(x = query, y = gtf, type = "any", nomatch = NULL)
    if(nrow(gene_models) == 0){
      message("No features found within the requested loci!")
      return(NULL)  
    }
  }
  
  
  gene_models_exon = gene_models[feature %in% c("exon", "transcript")]
  #gene_models_rest = gene_models[!feature %in% "exon"]
  
  feature_ids = data.table::tstrsplit(x = gene_models_exon$info, split = "; ")
  feature_id_names = lapply(feature_ids, function(x){
    #x = x[1:50] #sample rows
    x = unique(unlist(data.table::tstrsplit(x = x, split = " ", keep = 1)))
    x = x[complete.cases(x)]
    x[1]
  })
  names(feature_ids) = feature_id_names
  req_fields = c("gene_id", "transcript_id")
  req_fields = req_fields[req_fields %in% unlist(feature_id_names)]
  feature_ids = feature_ids[req_fields]
  
  feature_ids = sapply(feature_ids, function(x){
    gsub(pattern = "\"|;", replacement = "", x = unlist(data.table::tstrsplit(x = x, split = " ", keep = 2)))
  })
  feature_ids = data.frame(feature_ids)
  colnames(feature_ids) = c("name2", "name")
  gene_models_exon = cbind(gene_models_exon, feature_ids)
  #gene_models = data.table::rbindlist(list(gene_models_rest, gene_models_exon), use.names = TRUE, fill = TRUE)
  gene_models = gene_models_exon[order(as.numeric(as.character(start)))]
  #gene_models[,.(chr,start,end,strand, tx, gene)]
  
  if(nrow(gene_models) == 0){
    warning("No features found within the requested loci!")
    return(NULL)
  }else{
    if(!is.null(txname)){
      gene_models = gene_models[name %in% txname]
    }
    
    if(!is.null(genename)){
      gene_models = gene_models[name2 %in% genename]
    }
    
    if(nrow(gene_models) == 0){
      warning("Requested gene or transcript could not be found within the requested loci!")
      return(NULL)
    }
    
    gene_models = split(gene_models, as.factor(as.character(gene_models$name)))
    
    exon_tbls = lapply(seq_along(along.with = 1:length(gene_models)), function(idx){
      x = gene_models[[idx]]
      exon_start = as.numeric(as.character(x[feature %in% "exon"][, start]))
      exon_end = as.numeric(as.character(x[feature %in% "exon"][, end]))
      exon_tbl = data.frame(start = exon_start, end = exon_end)
      attributes(exon_tbl) = list(start = min(x[,start], na.rm = TRUE), end = max(x[,end], na.rm = TRUE), strand = unique(x[,strand]), tx = unique(x[, name]), gene = unique(x[, name2]), chr = unique(x[,chr]))
      exon_tbl
    })
  }
  exon_tbls
}

.make_genome_bed = function(refBuild = "hg19", tss = "start", up = 2500, down = 2500, op_dir = tempdir(), pc_genes = FALSE, for_profile = TRUE){
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  tss = match.arg(arg = tss, choices = c("start", "end"))
  
  temp_op_bed = tempfile(pattern = "profileplot_ucsc", tmpdir = op_dir, fileext = ".bed")
  
  message(paste0("Extracting gene models from UCSC:\n", "    build: ", refBuild))
  mydb = RMySQL::dbConnect(RMySQL::MySQL(), user = "genome", db = refBuild, host = "genome-mysql.soe.ucsc.edu")
  cmd = paste0("select chrom, txStart, txEnd, strand, name, name2 from refGene")
  ucsc = DBI::dbGetQuery(mydb, cmd)
  data.table::setDT(ucsc)
  colnames(ucsc) = c("chr", "start", "end", "strand", "tx_id", "gene_id")
  
  main_contigs = paste0("chr", c(1:22, "X", "Y"))
  ucsc = ucsc[chr %in% main_contigs]
  
  if(pc_genes){
    ucsc = ucsc[tx_id %like% "^NM"]
  }
  
  message("Fetched ", nrow(ucsc), " transcripts from ", nrow(ucsc[,.N,.(chr)]), " contigs")
  
  if(for_profile){
    #If it is only for profile plot where tss or tes are extended, we donyt extend manually here. Instead invert tss for negative strand txs. Let bwtool do the hard work
    ucsc_minus = ucsc[strand %in% "-"]
    colnames(ucsc_minus) = c("chr", "end", "start", "strand", "tx_id", "gene_id")
    ucsc_plus = ucsc[strand %in% "+"]
    ucsc_bed = data.table::rbindlist(l = list(ucsc_plus, ucsc_minus), use.names = TRUE, fill = TRUE)[,.(chr, start, end, tx_id, gene_id)]
  }else{
    #If for summary, extrend the regions
    ucsc_minus = ucsc[strand %in% "-"]
    if(nrow(ucsc_minus) > 0){
      if(tss == "start"){
        ucsc_minus[, bed_start := end-up]
        ucsc_minus[, bed_end := end+down]
      }else{
        ucsc_minus[, bed_start := start-up]
        ucsc_minus[, bed_end := start+down]
      }
    }
    
    ucsc_plus = ucsc[strand %in% "+"]
    if(nrow(ucsc_plus) > 0){
      if(tss == "start"){
        ucsc_plus[, bed_start := start-up]
        ucsc_plus[, bed_end := start+down]
      }else{
        ucsc_plus[, bed_start := end-up]
        ucsc_plus[, bed_end := end+down]
      }
    }
    
    ucsc_bed = data.table::rbindlist(
      l = list(ucsc_plus[, .(chr, bed_start, bed_end, tx_id, gene_id)], ucsc_minus[, .(chr, bed_start, bed_end, tx_id, gene_id)]),
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  colnames(ucsc_bed) = c("chr", "start", "end", "tx", "gene")
  data.table::setkey(x = ucsc_bed, chr, start, end)
  
  data.table::fwrite(x = ucsc_bed[,1:3], file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  data.table::fwrite(x = ucsc_bed, file = paste0(temp_op_bed, "2"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  list(temp_op_bed, paste0(temp_op_bed, "2"))
}

.make_bed = function(bed, op_dir = tempdir(), up = 2500, down = 2500, tss = "center", for_profile = FALSE){
  #bwtool tool requires only three columns
  
  tss = match.arg(arg = tss, choices = c("start", "end", "center"))
  
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  temp_op_bed = tempfile(pattern = "profileplot", tmpdir = op_dir, fileext = ".bed")
  
  if(is.data.frame(bed)){
    bed = data.table::as.data.table(x = bed)
    #data.table::setDT(x = bed)
    colnames(bed)[1:3] = c("chr", "start", "end")
    bed[, chr := as.character(chr)]
    bed[, start := as.numeric(as.character(start))]
    bed[, end := as.numeric(as.character(end))]
  }else if(file.exists(bed)){
    bed = data.table::fread(file = bed, select = list(character = 1, numeric = c(2, 3)), col.names = c("chr", "start", "end"))
    bed = bed[,.(chr, start, end)]
  }
  
  if(!for_profile){
    if(tss == "center"){
      bed[, focal_point := as.integer(apply(bed[,2:3], 1, mean))]
      bed[, bed_start := focal_point-up]
      bed[, bed_end := focal_point+down]
    }else if(tss == "start"){
      bed[, bed_start := start-up]
      bed[, bed_end := start+down]
    }else{
      bed[, bed_start := end-up]
      bed[, bed_end := end+down]
    }
    bed = bed[,.(chr, bed_start, bed_end)]
    data.table::setkey(x = bed, chr, bed_start, bed_end)
  }
  
  data.table::fwrite(x = bed, file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(temp_op_bed)
}

.bwt_mats = function(bw, binSize, bed, size, startFrom, op_dir){
  
  bn = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "",
            x = basename(bw), ignore.case = TRUE)
  message(paste0("Processing ", bn, ".."))
  
  bw = gsub(pattern = " ", replacement = "\\ ", x = bw, fixed = TRUE) #Change spaces with \ for unix style paths
  
  if(startFrom == "start"){
    cmd = paste0("bwtool matrix -starts -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))  
  }else if(startFrom == "end"){
    cmd = paste0("bwtool matrix -ends -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))
  }else{
    cmd = paste0("bwtool matrix -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))
  }
  print(cmd)

  system(command = cmd, intern = TRUE)
  paste0(op_dir, "/", bn, ".matrix")
}


.summarizeMats = function(mats = NULL, summarizeBy = 'mean', group = NULL, collapse_reps = FALSE){
  
  if(!is.null(group)){
    # gdf = data.table::data.table(sample = names(mats), condition = group)
    # gdf = gdf[order(group, sample)]
    group_u = unique(group)
    if(collapse_reps){
      summarizedMats = lapply(group_u, function(g){
        x = apply(data.table::rbindlist(l = mats[which(group == g)], fill = TRUE, use.names = TRUE), 2, summarizeBy, na.rm = TRUE)
        x
      })
      names(summarizedMats) = group_u
    }else{
      summarizedMats = lapply(mats[1:(length(mats))], function(x){
        if(!is.null(dim(x))){
          x = apply(x, 2, summarizeBy, na.rm = TRUE)
        }
        x
      })
    }
  }else{
    summarizedMats = lapply(mats[1:(length(mats))], function(x){
      if(!is.null(dim(x))){
        x = apply(x, 2, summarizeBy, na.rm = TRUE)
      }
      x
    })
  }
  
  #summarizedMats$param = mats$param
  summarizedMats
}


# estimate tandard deviation for CI
.estimateCI = function(mats = NULL, group = NULL, collapse_reps = FALSE){
  
  if(!is.null(group)){
    if(collapse_reps){
      group_u = unique(group)
      ciMats = lapply(group_u, function(g){
        x = apply(data.table::rbindlist(l = mats[which(group == g)], fill = TRUE, use.names = TRUE), 2, function(y){
          sd(y, na.rm = TRUE)/sqrt(length(y))
        })
        x
      })
      names(ciMats) = group_u
    }else{
      ciMats = lapply(mats[1:(length(mats))], function(x){
        if(!is.null(dim(x))){
          x = apply(x, 2, function(y){
            sd(y, na.rm = TRUE)/sqrt(length(y))
          })
        }
        x
      })
    }
  }else{
    ciMats = lapply(mats[1:(length(mats))], function(x){
      if(!is.null(dim(x))){
        x = apply(x, 2, function(y){
          sd(y, na.rm = TRUE)/sqrt(length(y))
        })
      }
      x
    })
  }
  
  ciMats
}

.order_by_sds = function(mat, keep_sd = FALSE){
  mat_sd = apply(as.matrix(mat), 1, sd, na.rm = TRUE) #order it based on SD
  mat_sd = sort(mat_sd, decreasing = TRUE)
  mat = mat[names(mat_sd),, drop = FALSE]
  mat
}

.extract_summary = function(bw, binSize, bed, op_dir){
  bn = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "",
            x = basename(bw), ignore.case = TRUE)
  message(paste0("    Processing ", bn, ".."))
  
  cmd = paste("bwtool summary -with-sum -keep-bed -header", bedSimple, bw, paste0(op_dir, bn, ".summary"))
  paste0(op_dir, "/", bn, ".summary")
}

.loci2df = function(loci){
  chr = as.character(unlist(data.table::tstrsplit(x = loci, spli = ":", keep = 1)))
  start = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-", keep = 1))
  start = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(start))))
  end = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-", keep = 2))
  end = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(end))))
  data.table::data.table(chr, start, end)
}

.check_bwtool = function(warn = FALSE){
  check = as.character(Sys.which(names = 'bwtool'))[1]
  if(check != ""){
    if(warn){
      message("Checking for bwtool installation")
      message(paste0("    All good! Found bwtool at: ", check))
    }else{
      return(invisible(0))
    }
  }else{
    stop("Could not locate bwtool. Download it from here: https://github.com/CRG-Barcelona/bwtool/releases")
  }
}

.check_mysql = function(warn = FALSE){
  check = as.character(Sys.which(names = 'mysql'))[1]
  if(check != ""){
    if(warn){
      message("Checking for mysql installation")
      message(paste0("    All good! Found mysql at: ", check))
    }else{
      return(invisible(0))
    }
  }else{
    stop("Could not locate mysql.\nInstall:\n apt install mysql-server [Debian]\n yum install mysql-server [centOS]\n brew install mysql [macOS]\n conda install -c anaconda mysql [conda]")
  }
}

.check_dt = function(){
  if(!requireNamespace("data.table", quietly = TRUE)){
    message("Could not find data.table library. Attempting to install..")
    install.packages("data.table")
  }
  suppressPackageStartupMessages(expr = library("data.table", quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE))
}

.check_windows = function(){
  if(Sys.info()["sysname"] == "Windows"){
    stop("Windows is not supported :(")
  }
}

.get_ucsc_hmm_tbls = function(){
  hmm = structure(
    list(
      TableName = c(
        "wgEncodeBroadHmmGm12878HMM",
        "wgEncodeBroadHmmH1hescHMM",
        "wgEncodeBroadHmmHepg2HMM",
        "wgEncodeBroadHmmHepg2HMM",
        "wgEncodeBroadHmmHsmmHMM",
        "wgEncodeBroadHmmHuvecHMM",
        "wgEncodeBroadHmmK562HMM",
        "wgEncodeBroadHmmNhekHMM",
        "wgEncodeBroadHmmNhlfHMM"
      ),
      cell = c(
        "GM12878",
        "H1-hESC",
        "HepG2",
        "HMEC",
        "HSMM",
        "HUVEC",
        "K562",
        "NHEK",
        "NHLF"
      ),
      Tier = c(1L,
               1L, 2L, 3L, 3L, 2L, 1L, 3L, 3L),
      Description = c(
        "B-lymphocyte, lymphoblastoid",
        "embryonic stem cells",
        "hepatocellular carcinoma",
        "mammary epithelial cells",
        "skeletal muscle myoblasts",
        "umbilical vein endothelial cells",
        "leukemia",
        "epidermal keratinocytes",
        "lung fibroblasts"
      ),
      Lineage = c(
        "mesoderm",
        "inner cell mass",
        "endoderm",
        "ectoderm",
        "mesoderm",
        "mesoderm",
        "mesoderm",
        "ectoderm",
        "endoderm"
      ),
      Tissue = c(
        "blood",
        "embryonic stem cell",
        "liver",
        "breast",
        "muscle",
        "blood vessel",
        "blood",
        "skin",
        "lung"
      ),
      Karyotype = c(
        "normal",
        "normal",
        "cancer",
        "normal",
        "normal",
        "normal",
        "cancer",
        "normal",
        "normal"
      ),
      Sex = c("F",
              "M", "M", "U", "U", "U", "F", "U", "U")
    ),
    row.names = c(NA,-9L),
    class = c("data.table", "data.frame")
  )
  
  hmm
}

.get_summaries_narrowPeaks = function(bigWigs, nthreads = 1, chr = NA, start = NA, end = NA){
  #bedSimple = temp_op_bed; bigWigs = list.files(path = "./", pattern = "bw"); op_dir = getwd(); nthreads = 1
  query = data.table::data.table(chromosome = chr, start = start, end = end, key = c("chromosome", "start", "end"))
  
  message(paste0("Extracting signals"))
  
  summaries = parallel::mclapply(bigWigs, FUN = function(bw){
    bn = unlist(data.table::tstrsplit(x =  basename(bw), split = "\\.", keep = 1))
    message(paste0("    Processing ", bn, " .."))
    bw = data.table::fread(file = bw)
    colnames(bw)[c(1:3, 5)] = c("chromosome", "start", "end", "max")
    data.table::setkeyv(x = bw, cols = c("chromosome", "start", "end"))
    bw = data.table::foverlaps(x = query, y = bw, type = "any", nomatch = NULL)
    bw[,.(chromosome, start, end, max)]
  }, mc.cores = nthreads)
  
  names(summaries) = unlist(data.table::tstrsplit(x =  basename(bigWigs), split = "\\.", keep = 1))
  summaries
}

.parse_loci = function(loci){
  chr = as.character(unlist(data.table::tstrsplit(x = loci, split = ":", keep = 1)))
  start = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-"))[1]
  start = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(start))))
  end = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-"))[2]
  end = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(end))))
  list(chr = chr, start = start, end = end)
}

#Plot profile mini version
.plot_profile_mini = function(plot_dat, index = 1, ylims = c(0, 2)){
  #size = as.character(plot_dat$param["size"])
  
  y = plot_dat[[index]]
  x = (1:(length(y)))/length(y)
  par(mar = c(1,2,1,1))
  plot(NA, axes = FALSE, xlim = c(-0.2, 1), ylim = ylims, xlab = NA, ylab = NA)
  #grid(col = "gray90")
  points(x, y, type = "l", lwd = 1.2)
  axis(side = 2, at = ylims, lwd = 1, line = -0.6, labels = NA)
  mtext(side = 2, at = ylims, lwd = 1, line = -0.05, text = round(ylims, digits = 2), cex = 0.8, font = 1, las = 2)
}

# Order matrices
.order_matrix = function(mats, sortBy = NULLs, k = NULL){
  
  mats_avg = lapply(mats, function(x){
    x = as.matrix(as.data.frame(x = x))
    apply(x, 1, mean, na.rm = TRUE)
  })
  
  mats_avg = as.data.frame(mats_avg)
  cluster_row_cut = NULL
  
  if(sortBy == "mean"){
    mats_avg$oall_avg = rowMeans(mats_avg, na.rm = TRUE)
    row_idx = order(mats_avg$oall_avg, decreasing = TRUE)
  }else if(sortBy == "median"){
    mats_avg$oall_avg = apply(mats_avg, 1, median, na.rm = TRUE)
    row_idx = order(mats_avg$oall_avg, decreasing = TRUE)
  }else if(sortBy == "hclust"){
    set.seed(seed = 1024)
    hc = hclust(d = dist(mats_avg))
    mats_avg$hc_order = hc$order
    if(!is.null(k)){
      mats_avg$cluster = cutree(tree = hc, k = k)
      mats_avg$row_idx = 1:nrow(mats_avg)
      mats_avg_spl = split(mats_avg, f = as.factor(as.character(mats_avg$cluster)))
      cluster_row_cut = unlist(lapply(mats_avg_spl, nrow))
      print(cluster_row_cut)
      mats_avg_spl = lapply(mats_avg_spl, function(x){
        xhc = hclust(d = dist(x[,1:(ncol(x)-3)]))
        x$row_idx2 = xhc$order
        x = x[order(x$row_idx2, decreasing = FALSE),, drop = FALSE]
        x
      })
      mats_avg_spl = data.table::rbindlist(l = mats_avg_spl, fill = TRUE, use.names = TRUE)
      #cluster_row_cut = cumsum(xhm[,.N, cluster][,N])
      row_idx = mats_avg_spl$row_idx
    }else{
      row_idx = mats_avg$hc_order
    }
  }
  
  mats = lapply(mats, function(x){
    x = as.data.frame(x = x)
    x = x[row_idx,,drop = FALSE]
    x
  })
  
  mats
}
#------------------------------------------------------------------------------------------------------------------------------------