# This R script contains functions for bigWig visualization
#
# track_extract() -> track_plot() : to generate IGV style track plots (aka locus plots) from bigWig files.
# profile_extract() -> profile_plot() : to generate profile-plots from bigWig files.
# extract_signal() -> pca_plot() : to perform PCA analysis based on genomic regions of interest or around TSS sites.
#
# Source code: https://github.com/PoisonAlien/trackplot
#
# MIT License
# Copyright (c) 2020 Anand Mayakonda <anandmt3@gmail.com>
#
# Changelog:
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

#' Extract bigWig track data for the given loci
#' @param bigWigs bigWig files. Default NULL. Required.
#' @param bed narrowPeak or broadPeak files. This argument is mutually exclusive with bigWigs.
#' @param loci target region to plot. Should be of format "chr:start-end". e.g; chr3:187715903-187752003 OR chr3:187,715,903-187,752,003
#' @param binsize bin size to extract signal. Default 50 (bps).
#' @param nthreads Default 1. Number of threads to use.
#' @param custom_names Default NULL and Parses from the file names.
#' @param query_ucsc Default TRUE. Queries UCSC and extracts gene models and cytoband for the loci. Requires `mysql` installation.
#' @param build Reference genome build. Default hg38
#' @import data.table
#' @examples
#' bigWigs = system.file("extdata", "bw", package = "trackplot") |> list.files(pattern = "\\.bw$", full.names = TRUE) 
#' oct4_loci = "chr6:31125776-31144789"
#' t = track_extract(bigWigs = bigWigs, loci = oct4_loci, build = "hg19")
#' @export
track_extract = function(bigWigs = NULL, bed = NULL, loci = NULL, binsize = 50, custom_names = NULL, nthreads = 1, query_ucsc = TRUE, build = "hg38"){
  
  if(all(is.null(bed), is.null(bigWigs))){
    stop("Please provide bigWigs or bed [narrowPeak or broadPeak] files")
  }
  
  if(all(!is.null(bed), !is.null(bigWigs))){
    stop("Please provide either bigWigs or bed [narrowPeak or broadPeak] files")
  }
  
  #Check if inputs are bigWig
  input_bw = ifelse(test = is.null(bed), yes = TRUE, no = FALSE)
  print(input_bw)
  
  .check_windows()
  if(input_bw){
    .check_bwtool()  
  }
  
  .check_dt()
  
  options(warn = -1)
  op_dir = tempdir() #For now
  
  message("Parsing loci..")
  if(is.null(loci)){
    stop("Missing loci. Provide a target region to plot.\n  Should be of format \"chr:start-end\". e.g; chr3:187715903-187752003 OR chr3:187,715,903-187,752,003")
  }else{
    chr = as.character(unlist(data.table::tstrsplit(x = loci, split = ":", keep = 1)))
    start = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-"))[1]
    start = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(start))))
    end = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, split = ":", keep = 2)), split = "-"))[2]
    end = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(end))))
  }
  
  if(start >= end){
    stop("End must be larger than Start!")
  }
  message("    Queried region: ", chr, ":", start, "-", end, " [", end-start, " bps]")
  
  message("Checking for files..")
  
  if(input_bw){
    input_files = bigWigs
  }else{
    input_files = bed
  }
  
  for(i in 1:length(input_files)){
    if(!file.exists(as.character(input_files)[i])){
      stop(paste0(as.character(input_files)[i], " does not exist!"))
    }
  }
  
  if(!is.null(custom_names)){
    if(length(custom_names) != length(input_files)){
      stop("Please provide names for all input files")
    }
  }
  
  if(input_bw){
    windows = .gen_windows(chr = chr, start = start, end = end, window_size = binsize, op_dir = op_dir)
    track_summary = .get_summaries(bedSimple = windows, bigWigs = bigWigs, op_dir = op_dir, nthreads = nthreads)  
  }else{
    track_summary = .get_summaries_narrowPeaks(bigWigs = bed, nthreads = nthreads, chr, start, end)  
  }
  
  
  if(!is.null(custom_names)){
    names(track_summary) = custom_names  
  }
  track_summary$loci = c(chr, start, end)
  
  #Extract gene models for this region
  if(query_ucsc){
    message("Querying UCSC genome browser for gene model and cytoband..")
    etbl = .extract_geneModel_ucsc(chr, start = start, end = end, refBuild = build, txname = NULL, genename = NULL)
    cyto = .extract_cytoband(chr = chr, refBuild = build)
  }else{
    cyto = etbl = NA
    
  }
  track_summary$etbl = etbl
  track_summary$cyto = cyto
  attr(track_summary, "is_bw") = input_bw
  attr(track_summary, "refbuild") = build
  
  track_summary
}

#' Summarize tracks per condition
#' @param summary_list Output from track_extract. Required.
#' @param condition Default NULL. Provide condition for each bigWig file 
#' @param stat can be `mean, median`, `max`, `min`. NAs are excluded. 
#' @export
track_summarize = function(summary_list = NULL, condition = NULL, stat = "mean"){
  
  if(is.null(summary_list)){
    stop("Missing input! Expecting output from track_extract()")
  }
  
  if(is.null(condition)){
    stop("Missing condition! Provide condition for each bigWig file.")
  }
  
  stat = match.arg(arg = stat, choices = c("mean", "median", "max", "min"))
  
  loci = summary_list$loci
  etbl = summary_list$etbl
  cyto = summary_list$cyto
  is_bw = attr(summary_list, "is_bw")
  build = attr(summary_list, "refbuild")
  
  summary_list$loci = NULL
  summary_list$etbl = NULL
  summary_list$cyto = NULL
  
  if(length(condition) != length(summary_list)){
    stop("Incorrect conditions! Provide condition for each bigWig file.")
  }
  
  names(summary_list) = condition
  
  summary_list = data.table::rbindlist(l = summary_list, use.names = TRUE, fill = TRUE, idcol = "sample_name")
  print(summary_list)
  
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
  summary_list$loci = loci
  summary_list$etbl = etbl
  summary_list$cyto = cyto
  attr(summary_list, "is_bw") = is_bw
  attr(summary_list, "refbuild") = build
  summary_list
}

#' Generate IGV style locus tracks with ease
#' @param summary_list Output from track_extract
#' @param draw_gene_track Default FALSE. If TRUE plots gene models overlapping with the queried region
#' @param show_ideogram Default TRUE. If TRUE plots ideogram of the target chromosome with query loci highlighted. Works only when `query_ucsc` is TRUE. 
#' @param txname transcript name to draw. Default NULL. Plots all transcripts overlapping with the queried region
#' @param genename gene name to draw. Default NULL. Plots all genes overlapping with the queried region
#' @param collapse_txs Default FALSE. Whether to collapse all transcripts belonging to same gene into a unified gene model
#' @param gene_model File with gene models. Can be a gtf file or UCSC file format. If you have read them into R as a data.frame, that works as well. Default NULL, automatically fetches gene models from UCSC server
#' @param isGTF Default FALSE. Set to TRUE if the `gene_model` is a gtf file.
#' @param groupAutoScale Default TRUE
#' @param groupScaleByCondition Scale tracks by condition
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
#' @examples
#' bigWigs = system.file("extdata", "bw", package = "trackplot") |> list.files(pattern = "\\.bw$", full.names = TRUE) |> rev()
#' oct4_loci = "chr6:31125776-31144789"
#' t = track_extract(bigWigs = bigWigs, loci = oct4_loci, build = "hg19")
#' trackplot::track_plot(summary_list = t)
#' @export
track_plot = function(summary_list = NULL,
                      draw_gene_track = TRUE,
                      show_ideogram = FALSE,
                      col = "gray70",
                      groupAutoScale = FALSE,
                      groupScaleByCondition = FALSE,
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
                      gene_model = NULL,
                      isGTF = FALSE,
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
                      left_mar = NULL
){
  
  if(is.null(summary_list)){
    stop("Missing input! Expecting output from track_extract()")
  }
  
  chr = summary_list$loci[1]
  start = as.numeric(summary_list$loci[2])
  end = as.numeric(summary_list$loci[3])
  etbl = summary_list$etbl
  cyto = summary_list$cyto
  is_bw = attr(summary_list, "is_bw")
  build = attr(summary_list, "refbuild")
  
  summary_list$loci = NULL
  summary_list$etbl = NULL
  summary_list$cyto = NULL
  
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
               loci_h = peaks_track_height, scale_track_height = scale_track_height)
  
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
    
    lapply(seq_along(peaks_data), function(idx){
      l_idx = peaks_data[[idx]]
      l_idx = data.table::foverlaps(x = query, y = l_idx, type = "any", nomatch = NULL)[,.(chr, start, end)]
      rect(xleft = start, ybottom = idx - 0.49, xright = end, ytop = idx - 0.51, col = "gray90", border = NA)
      rect(xleft = l_idx$start, ybottom = idx - 0.9, xright = l_idx$end, ytop = idx - 0.1, col = "#34495e", border = NA)
      text(x = start, y = idx - 0.5, labels = names(peaks_data)[idx], adj = 1.2, xpd = TRUE)
    })
  }
  
  #Draw bigWig signals
  if(is_bw){
    lapply(1:length(summary_list), function(idx){
      x = summary_list[[idx]]
      if(show_axis){
        par(mar = c(0.5, left_mar, 2, 1))
      }else{
        par(mar = c(0.5, left_mar, 2, 1))  
      }
      
      plot(NA, xlim = c(start, end), ylim = c(plot_height_min[idx], plot_height[idx]), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
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
      
    })
    
  }else{
    lapply(1:length(summary_list), function(idx){
      x = summary_list[[idx]]
      if(show_axis){
        par(mar = c(0, left_mar, 0, 1))
      }else{
        par(mar = c(0.5, left_mar, 1, 1))  
      }
      
      plot(NA, xlim = c(start, end), ylim = c(0, 1), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
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
      
    })
    
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
    
    if(isGTF){
      etbl = .parse_gtf(gtf = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)
    }else{
      etbl = .make_exon_tbl(gene_models = etbl, txname = txname, genename = genename)
    }
    
    if(!is.null(etbl)){
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
  if(show_ideogram){
    par(mar = c(0.2, 1, 0, 1))
    plot(NA, xlim = c(0, max(cyto$end)), ylim = c(0, 1), axes = FALSE, frame.plot = FALSE, xlab = NA, ylab = NA)
    rect(xleft = cyto$start, ybottom = 0.1, xright = cyto$end, ytop = 0.6, col = cyto$color, border = "#34495e")
    rect(xleft = start, ybottom = 0, xright = end, ytop = 0.7, col = "#d35400", lwd = 2, border = "#d35400")
    text(x = 0, y = 0.5, labels = chr, adj = 1.2, font = 2, xpd = TRUE)
  }
}

# profileplot is an ultra-fast, simple, and minimal dependency R script to generate profile-plots from bigWig files
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
#' @export

profile_extract = function(bigWigs = NULL, bed = NULL, binSize = 50, startFrom = "center",
                           up = 2500, down = 2500, ucsc_assembly = NULL, nthreads = 4, 
                           custom_names = NULL){
  .check_windows()
  .check_bwtool()
  .check_dt()
  
  if(is.null(bigWigs)){
    stop("Provide at-least one bigWig file")
  }
  
  message("Checking for files..")
  for(i in 1:length(bigWigs)){
    if(!file.exists(as.character(bigWigs)[i])){
      stop(paste0(as.character(bigWigs)[i], " does not exist!"))
    }
  }
  
  if(!is.null(custom_names)){
    if(length(custom_names) != length(bigWigs)){
      stop("Please provide names for all bigWigs")
    }
  }
  
  op_dir = tempdir() #For now
  
  
  if(is.null(bed)){
    if(is.null(ucsc_assembly)){
      stop("Please provide either a BED file or ucsc_assembly name")
    }else{
      if(length(ucsc_assembly) > 1){
        warning("Multiple assemblies provided. Using first one..")
        ucsc_assembly = ucsc_assembly[1]
      }
      bed = .make_genome_bed(refBuild = ucsc_assembly, up = as.numeric(up), down = as.numeric(down), tss = startFrom, op_dir = op_dir, pc_genes = FALSE)
    }
  }else{
    startFrom = match.arg(arg = startFrom, choices = c("start", "end", "center"))
    bed = .make_bed(bed = bed, op_dir = op_dir, up = as.numeric(up), down = as.numeric(up), tss = startFrom)
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
  #return(sig_list)
  
  sig_list$args = c(up, down)
  
  #Remove intermediate files
  lapply(mats, function(x) system(command = paste0("rm ", x), intern = TRUE))
  
  sig_list
}


#' Draw a profile plot
#' @param sig_list Output generated from profile_extract
#' @param color Manual colors for each bigWig. Default NULL. 
#' @param condition Default. Condition associated with each bigWig. Lines will colord accordingly.
#' @param condition_colors Manual colors for each level in condition. Default NULL. 
#' @param collapse_replicates Default FALSE. If TRUE and when `condition` is given, collapse signals samples belonging to same condition
#' @param plot_se Default FALSE. If TRUE plots standard error shading
#' @param line_size Default 1
#' @param legend_fs Legend font size. Default 1
#' @param axis_fs Axis font size. Default 1
#' @export

profile_plot = function(sig_list = NULL, color = NULL, condition = NULL, condition_colors = NULL, 
                        collapse_replicates = FALSE, plot_se = FALSE, line_size = 1, legend_fs = 1, axis_fs = 1){
  
  
  if(is.null(sig_list)){
    stop("Missing input! Expecting output from profile_extract()")
  }
  
  up = as.numeric(sig_list$args[1])
  down = as.numeric(sig_list$args[2])
  sig_list$args = NULL
  
  message("Summarizing..")
  sig_summary = .summarizeMats(mats = sig_list, group = condition, collapse_reps = collapse_replicates)
  
  if(is.null(color)){
    color = c("#A6CEE3FF", "#1F78B4FF", "#B2DF8AFF", "#33A02CFF", "#FB9A99FF", 
              "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF", "#CAB2D6FF", "#6A3D9AFF", 
              "#FFFF99FF", "#9E0142FF", "#D53E4FFF", "#F46D43FF", "#000000FF", 
              "#EE82EEFF", "#4169E1FF", "#7B7060FF", "#535C68FF")
    color = color[1:length(sig_list)]
  }
  
  
  sig_se_summary = NULL
  if(plot_se){
    sig_se_summary = .estimateCI(mats = sig_list, group = condition, collapse_reps = collapse_replicates)
  }
  
  if(!is.null(condition)){
    group_df = data.table::data.table(sample = names(sig_summary), condition = condition)
    if(is.null(condition_colors)){
      condition_colors = c("#A6CEE3FF", "#1F78B4FF", "#B2DF8AFF", "#33A02CFF", "#FB9A99FF", 
                           "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF", "#CAB2D6FF", "#6A3D9AFF", 
                           "#FFFF99FF", "#9E0142FF", "#D53E4FFF", "#F46D43FF", "#000000FF", 
                           "#EE82EEFF", "#4169E1FF", "#7B7060FF", "#535C68FF")[1:nrow(group_df[,.N,condition])]
    }
    names(condition_colors) = group_df[,.N,condition][,condition]
    color = condition_colors[group_df$condition]
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
  par(mar = c(3, 3, 2, 1))
  plot(NA, xlim = c(0, x_max), ylim = c(min(ylabs), max(ylabs)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  abline(h = ylabs, v = pretty(xticks), col = "gray90", lty = 2)
  
  lapply(1:length(sig_summary), function(idx){
    points(sig_summary[[idx]], type = 'l', lwd = line_size, col = color[idx])
    if(plot_se){
      polygon(x = c(1:length(sig_summary[[idx]]), rev(1:length(sig_summary[[idx]]))),
              y = c(sig_summary[[idx]]-sig_se_summary[[idx]], rev(sig_summary[[idx]]+sig_se_summary[[idx]])),
              col = grDevices::adjustcolor(col = color[idx], #color[names(mat_list)[[i]]],
                                           alpha.f = 0.4), border = NA)
    }
  })
  axis(side = 1, at = xticks, labels = xlabs, cex.axis = axis_fs)
  axis(side = 2, at = ylabs, las = 2, cex.axis = axis_fs)
  
  if(!is.null(condition)){
    legend(x = "topright", legend = unique(names(color)), col = unique(color), bty = "n", lty = 1, lwd = 1.2, cex = legend_fs, xpd = TRUE)
  }else{
    legend(x = "topright", legend = names(sig_summary), col = color, bty = "n", lty = 1, lwd = 1.2, cex = legend_fs, xpd = TRUE)
  }
  
  invisible(list(mean_signal = sig_summary, std_err = sig_se_summary, color_codes = color, xticks = xticks, xlabs = xlabs))
  
}

# bwpcaplot function to perform PCA analysis based on genomic regions of interest or around TSS sites.

#' Extract area under the curve for every peak from from given bigWig files.
#' @param bigWigs bigWig files. Default NULL. Required.
#' @param bed bed file or a data.frame with first 3 column containing chromosome, star, end positions. 
#' @param binSize bin size to extract signal. Default 50 (bps). Should be >1
#' @param ucsc_assembly If `bed` file not provided, setting `ucsc_assembly` to ref genome build will fetch transcripts from UCSC genome browser. e.g; 'hg19'
#' @param ucsc_startFrom Default "center". Can be "center", "start" or "end"
#' @param ucsc_up extend upstream by this many bps from `startFrom`. Default 2500
#' @param ucsc_down extend downstream by this many bps from `startFrom`. Default 2500
#' @param nthreads Default 4
#' @param custom_names Default NULL and Parses from the file names.
#' @export
extract_summary = function(bigWigs = NULL, bed = NULL, binSize = 50,
                           nthreads = 4, ucsc_assembly = NULL, ucsc_startFrom = "start", ucsc_up = 2500, ucsc_down = 2500, custom_names = NULL){
  
  .check_windows()
  .check_bwtool()
  .check_dt()
  
  if(is.null(bigWigs)){
    stop("Missing bigWig files")
  }
  
  message("Checking for files..")
  for(i in 1:length(bigWigs)){
    if(!file.exists(as.character(bigWigs)[i])){
      stop(paste0(as.character(bigWigs)[i], " does not exist!"))
    }
  }
  
  if(!is.null(custom_names)){
    if(length(custom_names) != length(bigWigs)){
      stop("Please provide names for all bigWigs")
    }
  }else{
    custom_names = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "", x = basename(bigWigs))
  }
  
  op_dir = tempdir() #For now
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  
  if(!is.null(ucsc_assembly)){
    #If assembly is requested, fetch protein coding transcripts (+/- 2500bp.
    temp_op_bed = .make_genome_bed(refBuild = ucsc_assembly, tss = "start", pc_genes = TRUE, up = ucsc_up, down = ucsc_down)
  }else{
    temp_op_bed = tempfile(pattern = "pcaplot", tmpdir = op_dir, fileext = ".bed")
    if(is.data.frame(bed)){
      bed = data.table::as.data.table(x = bed)
      #data.table::setDT(x = bed)
      colnames(bed)[1:3] = c("chr", "start", "end")
      bed[, chr := as.character(chr)]
      bed[, start := as.numeric(as.character(start))]
      bed[, end := as.numeric(as.character(end))]
    }else if(file.exists(bed)){
      bed = data.table::fread(file = bed, select = list(character = 1, numeric = c(2, 3)), col.names = c("chr", "start", "end"))
    }
    data.table::fwrite(x = bed[,.(chr, start, end)], file = temp_op_bed, sep = "\t", col.names = FALSE)
  }
  
  message("Extracting BED summaries..")
  summaries = parallel::mclapply(bigWigs, FUN = function(bw){
    bn = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "", x = basename(bw))
    message(paste0("    Processing ", bn, " .."))
    cmd = paste("bwtool summary -with-sum -keep-bed -header", temp_op_bed, bw, paste0(op_dir, bn, ".summary"))
    system(command = cmd, intern = TRUE)
    paste0(op_dir, bn, ".summary")
  }, mc.cores = nthreads)
  
  summary_list = lapply(summaries, function(x){
    x = data.table::fread(x)
    colnames(x)[1] = 'chromosome'
    x = x[,.(chromosome, start, end, size, sum)]
    x[,id := paste0(chromosome, ":", start, "-", end)]
    x
  })
  
  names(summary_list) = custom_names
  
  summary_list
  
}

#' Draw a PCA plot
#' @param summary_list output from extract_summary
#' @param top Top most variable peaks to consider for PCA. Default 1000
#' @param xpc Default PC1
#' @param ypc Default PC2
#' @param color Manual colors for each bigWig. Default NULL. 
#' @param condition Default. Condition associated with each bigWig. Lines will colord accordingly.
#' @param condition_colors Manual colors for each level in condition. Default NULL. 
#' @param condition_pch Manual plotting character for each level in condition. Default 21. 
#' @param show_cree If TRUE draws a cree plot. Default TRUE
#' @param size Point size. Default 1
#' @param lab_size Font size for labels. Default 1
#' @param legend_size Default 1
#' @param log2 log transform data? Default FALSE. IF TRUE, adds a small positive value and log2 converts.
#' @param legendpos Default topright
#' @export
pca_plot = function(summary_list = NULL, top = 1000, xpc = "PC1", ypc = "PC2", size = 1, color = NULL, condition = NULL, condition_colors = NULL, 
                    condition_pch = NULL, show_cree = TRUE, lab_size = 1, log2 = FALSE, legendpos = "topright", legend_size = 1){
  
  
  if(is.null(summary_list)){
    stop("Missing input! Expecting output from extract_summary()")
  }
  
  sum_tbl = data.table::rbindlist(l = summary_list, idcol = "sample", use.names = TRUE, fill = TRUE)
  sum_tbl = data.table::dcast(data = sum_tbl, formula =  id ~ sample, value.var = "sum", fun.aggregate = max)
  data.table::setDF(x = sum_tbl, rownames = sum_tbl$id)
  sum_tbl$id = NULL
  sum_tbl
  
  if(is.null(color)){
    color = c("#A6CEE3FF", "#1F78B4FF", "#B2DF8AFF", "#33A02CFF", "#FB9A99FF", 
              "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF", "#CAB2D6FF", "#6A3D9AFF", 
              "#FFFF99FF", "#9E0142FF", "#D53E4FFF", "#F46D43FF", "#000000FF", 
              "#EE82EEFF", "#4169E1FF", "#7B7060FF", "#535C68FF")
    color = color[1:length(summary_list)]
  }
  
  if(!is.null(condition)){
    group_df = data.table::data.table(sample = names(summary_list), condition = condition)
    if(is.null(condition_colors)){
      condition_colors = c("#A6CEE3FF", "#1F78B4FF", "#B2DF8AFF", "#33A02CFF", "#FB9A99FF", 
                           "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF", "#CAB2D6FF", "#6A3D9AFF", 
                           "#FFFF99FF", "#9E0142FF", "#D53E4FFF", "#F46D43FF", "#000000FF", 
                           "#EE82EEFF", "#4169E1FF", "#7B7060FF", "#535C68FF")[1:nrow(group_df[,.N,condition])]
    }
    names(condition_colors) = group_df[,.N,condition][,condition]
    group_df$color = condition_colors[group_df$condition]
  
    if(!is.null(condition_pch)){
      names(condition_pch) = group_df[,.N,condition][,condition]
      group_df$pch = condition_pch[group_df$condition]
    }else{
      group_df$pch = 19
    }
    
  }else{
    group_df = data.table::data.table(sample = names(summary_list), color = color)
    group_df$pch = 21
  }
  
  if(log2){
    sum_tbl = log2(x = sum_tbl + 0.1)
  }
  sum_tbl = sum_tbl[complete.cases(sum_tbl),, drop = FALSE]
  
  if(nrow(sum_tbl) < top){
    pca = prcomp(t(.order_by_sds(mat = sum_tbl)))
  }else{
    pca = prcomp(t(.order_by_sds(mat = sum_tbl)[1:top,]))  
  }
  
  pca_dat = as.data.frame(pca$x)
  pca_var_explained = pca$sdev^2/sum(pca$sdev^2)
  names(pca_var_explained) = paste0("PC", 1:length(pca_var_explained))
  data.table::setDT(x = pca_dat, keep.rownames = "sample")
  data.table::setDF(x = pca_dat)
  pca_dat = merge(pca_dat, group_df, by = 'sample')
  attr(pca_dat, "percentVar") <- round(pca_var_explained, digits = 3)
  #print(head(pca_dat))
  if(show_cree){
    lo = layout(mat = matrix(data = c(1, 2), ncol = 2))
  }
  
  grid_cols = "gray90"
  par(mar = c(3, 4, 1, 1))
  plot(NA, axes = FALSE, xlab = NA, ylab = NA, cex = 1.2, xlim = range(pretty(pca_dat[, xpc])), ylim = range(pretty(pca_dat[, ypc])))
  abline(h = pretty(pca_dat[, xpc]), v = pretty(pca_dat[, ypc]), col = grid_cols, lty = 2, lwd = 0.1)
  abline(h = 0, v = 0, col = grid_cols, lty = 2, lwd = 0.8)
  points(x = pca_dat[, xpc], y = pca_dat[, ypc], col = pca_dat$color, bg = pca_dat$color, pch = pca_dat$pch, cex = size)
  axis(side = 1, at = pretty(pca_dat[, xpc]), cex.axis = 0.8)
  axis(side = 2, at = pretty(pca_dat[, ypc]), las = 2, cex.axis = 0.8)
  mtext(text = paste0(xpc, " [", round(pca_var_explained[xpc], digits = 2), "]"), side = 1, line = 2, cex = 0.8)
  mtext(text = paste0(ypc, " [", round(pca_var_explained[ypc], digits = 2), "]"), side = 2, line = 2, cex = 0.8)
  text(x = pca_dat[, xpc], y = pca_dat[, ypc], labels = pca_dat$sample, pos = 3, col = pca_dat$color, xpd = TRUE, cex = lab_size)
  #title(main = NA, sub = paste0("N = ", top, " peaks"), adj = 0, outer = TRUE)
  
  data.table::setDT(x = pca_dat)
  if(!is.null(condition)){
    leg_dat = pca_dat[,.N,.(condition, color, pch)]
    legend(x = legendpos, legend = leg_dat$condition, col = leg_dat$color, 
           bty = "n", pch = leg_dat$pch, cex = legend_size, xpd = TRUE, ncol = 2)
  }
  
  if(show_cree){
    par(mar = c(3, 4, 1, 4))
    b = barplot(height = pca_var_explained, names.arg = NA, col = "#2c3e50", ylim = c(0, 1), las = 2, axes = FALSE)
    points(x = b, y = cumsum(pca_var_explained), type = 'l', lty = 2, lwd = 1.2, xpd = TRUE, col = "#c0392b")
    points(x = b, y = cumsum(pca_var_explained), type = 'p', pch = 19, xpd = TRUE, col = "#c0392b")
    mtext(text = paste0("PC", 1:length(pca_var_explained)), side = 1, at = b, las = 2, line = 0.5, cex = 0.8)
    axis(side = 2, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
    mtext(text = "var. explained", side = 2, line = 2.5)
    axis(side = 4, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
    mtext(text = "cumulative var. explained", side = 4, line = 2.5)
  }
  
  invisible(x = list(pca_data = pca_dat, bed_summary = sum_tbl))
}

#------------------------------------------------------------------------------------------------------------------------------------
#                                                   Undocumented Accessory functions
#------------------------------------------------------------------------------------------------------------------------------------

.make_layout = function(ntracks, ntracks_h = 3, cytoband = TRUE, cytoband_h = 1, genemodel = TRUE, genemodel_h = 1, chrHMM = TRUE, chrHMM_h = 1, loci = TRUE, loci_h = 2, scale_track_height = 1){
  
  case = NULL
  #print(paste(cytoband, genemodel, chrHMM, loci, sep = ", "))
  
  if(cytoband & genemodel & chrHMM & loci){
    #print("all")
    case = 1
  }else if(cytoband & genemodel & chrHMM == FALSE & loci){
    #print("no hmm")
    case = 2
  }else if(cytoband & genemodel & chrHMM == FALSE & loci == FALSE){
    #print("no hmm and no loci")
    case = 3
  }else if(cytoband & genemodel == FALSE & chrHMM == FALSE & loci == FALSE){
    #print("no hmm and no loci and no genemodel")
    case = 4
  }else if(cytoband == FALSE & genemodel == FALSE & chrHMM == FALSE & loci == FALSE){
    #print("no hmm and no loci and no genemodel and no cytoband")
    case = 5
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == TRUE & loci == FALSE){
    #print("no loci and no cytoband")
    case = 6
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == TRUE & loci == TRUE){
    #print("no cytoband")
    case = 7
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == FALSE & loci == TRUE){
    #print("no cytoband no chrHMM")
    case = 8
  }else if(cytoband == FALSE & genemodel == TRUE & chrHMM == FALSE & loci == FALSE){
    #print("no cytoband no chrHMM no loci")
    case = 9
  }else if(cytoband == TRUE & genemodel == TRUE & chrHMM == TRUE & loci == FALSE){
    #print("no loci")
    case = 10
  }else if(cytoband == TRUE & genemodel == FALSE & chrHMM == FALSE & loci == TRUE){
    #print("no genemodel no chrhmm")
    case = 11
  }else{
    #print("Something is wrong!")
  }
  
  if(case == 1){
    lo = layout(mat = matrix(data = seq_len(ntracks+5)), heights = c(loci_h, rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height, cytoband_h))
  }else if(case == 2){
    lo = layout(mat = matrix(data = seq_len(ntracks+4)), heights = c(loci_h, rep(ntracks_h, ntracks), genemodel_h, scale_track_height, cytoband_h))
  }else if(case == 3){
    lo = layout(mat = matrix(data = seq_len(ntracks+3)), heights = c(rep(ntracks_h, ntracks), genemodel_h, scale_track_height, cytoband_h))
  }else if(case == 4){
    lo = layout(mat = matrix(data = seq_len(ntracks+2)), heights = c(rep(ntracks_h, ntracks), scale_track_height, cytoband_h))
  }else if(case == 5){
    lo = layout(mat = matrix(data = seq_len(ntracks+1)), heights = c(rep(ntracks_h, ntracks), scale_track_height))
  }else if(case == 6){
    lo = layout(mat = matrix(data = seq_len(ntracks+3)), heights = c(rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height))
  }else if(case == 7){
    lo = layout(mat = matrix(data = seq_len(ntracks+4)), heights = c(loci_h, rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height))
  }else if(case == 8){
    lo = layout(mat = matrix(data = seq_len(ntracks+3)), heights = c(loci_h, rep(ntracks_h, ntracks), genemodel_h, scale_track_height))
  }else if(case == 9){
    lo = layout(mat = matrix(data = seq_len(ntracks+2)), heights = c(rep(ntracks_h, ntracks), genemodel_h, scale_track_height))
  }else if(case == 10){
    lo = layout(mat = matrix(data = seq_len(ntracks+4)), heights = c(rep(ntracks_h, ntracks), chrHMM_h, genemodel_h, scale_track_height, cytoband_h))
  }else if(case == 11){
    lo = layout(mat = matrix(data = seq_len(ntracks+3)), heights = c(loci_h, rep(ntracks_h, ntracks), scale_track_height, cytoband_h))
  }
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

.extract_cytoband = function(chr = NULL, refBuild = "hg19"){
  
  if(!grepl(pattern = "^chr", x = chr)){
    message("Adding chr prefix to target chromosome for UCSC query..")
    chr = paste0("chr", chr)
  }
  
  cmd = paste0(
    "mysql --user genome --host genome-mysql.cse.ucsc.edu -NAD ",
    refBuild,
    " -e 'select chrom, chromStart, chromEnd, name, gieStain from cytoBand WHERE chrom =\"",
    chr,
    "\"'"
  )
  message(paste0("Extracting cytobands from UCSC:\n", "    chromosome: ", chr, "\n", "    build: ", refBuild, "\n    query: ", cmd))
  
  cyto = data.table::fread(cmd = cmd, colClasses = c("character", "numeric", "numeric", "character", "character"))
  colnames(cyto) = c("chr", "start", "end", "band", "stain")
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
  
  .check_mysql()
  
  cmd = paste0("mysql --user genome --host genome-mysql.cse.ucsc.edu -NAD ",  refBuild,  " -e 'select chrom, chromStart, chromEnd, name from ", tbl, " WHERE chrom =\"", tar_chr, "\"'")
  message(paste0("Extracting chromHMM from UCSC:\n", "    chromosome: ", tar_chr, "\n", "    build: ", refBuild, "\n    query: ", cmd))
  #system(command = cmd)
  ucsc = data.table::fread(cmd = cmd)
  if(nrow(ucsc) == 0){
    message("No features found within the requested loci!")
    return(NULL)
  }
  colnames(ucsc) = c("chr", "start", "end", "name")
  if(!grepl(pattern = "^chr", x = chr)){
    ucsc[, chr := gsub(pattern = "^chr", replacement = "", x = chr)]
  }
  data.table::setkey(x = ucsc, chr, start, end)
  
  query = data.table::data.table(chr = chr, start = start, end = end)
  data.table::setkey(x = query, chr, start, end)
  
  data.table::foverlaps(x = query, y = ucsc, type = "any", nomatch = NULL)[,.(chr, start, end, name)]
}

.extract_geneModel_ucsc = function(chr, start = NULL, end = NULL, refBuild = "hg19", txname = NULL, genename = NULL){
  .check_mysql()
  op_file = tempfile(pattern = "ucsc", fileext = ".tsv")
  
  if(!grepl(pattern = "^chr", x = chr)){
    message("Adding chr prefix to target chromosome for UCSC query..")
    tar_chr = paste0("chr", chr)
  }else{
    tar_chr = chr
  }
  
  cmd = paste0("mysql --user genome --host genome-mysql.cse.ucsc.edu -NAD ",  refBuild,  " -e 'select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds from refGene WHERE chrom =\"", tar_chr, "\"'")
  message(paste0("Extracting gene models from UCSC:\n", "    chromosome: ", tar_chr, "\n", "    build: ", refBuild, "\n    query: ", cmd))
  #system(command = cmd)
  ucsc = data.table::fread(cmd = cmd)
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
    message("No features found within the requested loci! If you are not sure why..\n    1.Make sure there are no discripancies in chromosome names i.e, chr prefixes\n")
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
  tx_tbl = tx_tbl[,id := paste0(start, ":", end)][!duplicated(id), ,.(gene)]
  
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
  
  query = data.table::data.table(chr, start, end)
  data.table::setkey(x = query, chr, start, end)
  
  gene_models = data.table::foverlaps(x = query, y = gtf, type = "any", nomatch = NULL)
  
  if(nrow(gene_models) == 0){
    message("No features found within the requested loci! If you are not sure why..\n    1.Make sure there are no discripancies in chromosome names i.e, chr prefixes\n")
    return(NULL)  
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
  colnames(feature_ids) = c("gene", "tx")
  gene_models_exon = cbind(gene_models_exon, feature_ids)
  #gene_models = data.table::rbindlist(list(gene_models_rest, gene_models_exon), use.names = TRUE, fill = TRUE)
  gene_models = gene_models_exon[order(as.numeric(as.character(start)))]
  
  if(nrow(gene_models) == 0){
    warning("No features found within the requested loci!")
    return(NULL)
  }else{
    if(!is.null(txname)){
      gene_models = gene_models[tx %in% txname]
    }
    
    if(!is.null(genename)){
      gene_models = gene_models[gene %in% genename]
    }
    
    if(nrow(gene_models) == 0){
      warning("Requested gene or transcript could not be found within the requested loci!")
      return(NULL)
    }
    
    gene_models = split(gene_models, as.factor(as.character(gene_models$tx)))
    
    exon_tbls = lapply(seq_along(along.with = 1:length(gene_models)), function(idx){
      x = gene_models[[idx]]
      exon_start = as.numeric(as.character(x[feature %in% "exon"][, start]))
      exon_end = as.numeric(as.character(x[feature %in% "exon"][, end]))
      exon_tbl = data.frame(start = exon_start, end = exon_end)
      attributes(exon_tbl) = list(start = min(x[,start], na.rm = TRUE), end = max(x[,end], na.rm = TRUE), strand = unique(x[,strand]), tx = unique(x[, tx]), gene = unique(x[, gene]))
      exon_tbl
    })
  }
  exon_tbls
}

.make_genome_bed = function(refBuild = "hg19", tss = "start", up = 2500, down = 2500, op_dir = tempdir(), pc_genes = FALSE){
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  tss = match.arg(arg = tss, choices = c("start", "end"))
  
  temp_op_bed = tempfile(pattern = "profileplot_ucsc", tmpdir = op_dir, fileext = ".bed")
  
  cmd = paste0("mysql --user genome --host genome-mysql.cse.ucsc.edu -NAD ",  refBuild,  " -e 'select chrom, txStart, txEnd, strand, name, name2 from refGene'")
  message(paste0("Extracting gene models from UCSC:\n", "    build: ", refBuild, "\n    query: ", cmd))
  #system(command = cmd)
  ucsc = data.table::fread(cmd = cmd)
  colnames(ucsc) = c("chr", "start", "end", "strand", "tx_id", "gene_id")
  
  main_contigs = paste0("chr", c(1:22, "X", "Y"))
  ucsc = ucsc[chr %in% main_contigs]
  
  if(pc_genes){
    ucsc = ucsc[tx_id %like% "^NM"]
  }
  
  message("Fetched ", nrow(ucsc), " transcripts from ", nrow(ucsc[,.N,.(chr)]), " contigs")
  
  if(tss == "end"){
    colnames(ucsc)[2:3] = c("end", "start")
  }
  
  ucsc_minus = ucsc[strand %in% "-"]
  if(nrow(ucsc_minus) > 0){
    ucsc_minus[, bed_start := end-up]
    ucsc_minus[, bed_end := end+down]
  }
  
  ucsc_plus = ucsc[strand %in% "+"]
  if(nrow(ucsc_plus) > 0){
    ucsc_plus[, bed_start := start-up]
    ucsc_plus[, bed_end := start+down]
  }
  
  ucsc_bed = data.table::rbindlist(l = list(ucsc_plus[, .(chr, bed_start, bed_end)], ucsc_minus[, .(chr, bed_start, bed_end)]), use.names = TRUE, fill = TRUE)
  colnames(ucsc_bed) = c("chr", "start", "end")
  data.table::setkey(x = ucsc_bed, chr, start, end)
  
  data.table::fwrite(x = ucsc_bed[,1:3], file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  temp_op_bed
}

.make_bed = function(bed, op_dir = tempdir(), up = 2500, down = 2500, tss = "center"){
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
  }
  
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
  
  data.table::fwrite(x = bed, file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(temp_op_bed)
}

.bwt_mats = function(bw, binSize, bed, size, startFrom, op_dir){
  
  bn = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "",
            x = basename(bw), ignore.case = TRUE)
  message(paste0("Processing ", bn, ".."))
  
  bw = gsub(pattern = " ", replacement = "\\ ", x = bw, fixed = TRUE) #Change spaces with \ for unix style paths
  
  cmd = paste0("bwtool matrix -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))
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

.check_bwtool = function(warn = TRUE){
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


#------------------------------------------------------------------------------------------------------------------------------------