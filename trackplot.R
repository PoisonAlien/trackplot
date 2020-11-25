#' Generate IGV style locus tracks with ease
#' @param bigWigs bigWig files. Default NULL. Required.
#' @param loci target region to plot. Should be of format "chr:start-end". e.g; chr3:187715903-187752003
#' @param chr chromosome of interest. Default NULL. This argument is mutually exclusive with `loci`
#' @param start start position. Default NULL. This argument is mutually exclusive with `loci`
#' @param end end interest. Default NULL. This argument is mutually exclusive with `loci`
#' @param nthreads Default 1.
#' @param binsize bin size to extract signal. Default 50 (bpbs).
#' @param draw_gene_track Default FALSE. If TRUE plots gene models overlapping with the queried region
#' @param gene_model File with gene models. Can be a gtf file or UCSC file format. Default NULL, automatically fetches gene models from UCSC server
#' @param isGTF Default FALSE. Set to TRUE if the `gene_model` is a gtf file.
#' @param tx transcript name to draw. Default NULL. Plots all transcripts overlapping with the queried region
#' @param gene gene name to draw. Default NULL. Plots all genes overlapping with the queried region
#' @param gene_fsize Font size. Default 0.8
#' @param gene_track_width Default 1 
#' @param scale_track_width Default 1 
#' @param query_ucsc Default FALSE. But switches to TREI when `gene_model` is FALSE
#' @param build Genome build. Default hg19
#' @param col Color for tracks. Default `#2f3640`. Multiple colors can be provided for each track
#' @param groupAutoScale Default TRUE
#' @param show_axis Default FALSE
#' @param custom_names Default NULL and Parses from the file names.
#' @param mark_regions genomic regions to highlight. A data.frame with three columns containing chr, start and end positions.
#' @param regions_track_width Default 1
#' 
trackplot = function(bigWigs = NULL,
                  loci = NULL,
                  chr = NULL,
                  start = NULL,
                  end = NULL,
                  nthreads = 1,
                  binsize = 50,
                  draw_gene_track = FALSE,
                  gene_model = NULL,
                  isGTF = FALSE,
                  tx = NULL,
                  gene = NULL,
                  gene_fsize = 0.8,
                  gene_track_width = 2,
                  scale_track_width = 1,
                  query_ucsc = FALSE,
                  build = "hg19",
                  col = "#2f3640",
                  groupAutoScale = TRUE,
                  show_axis = FALSE,
                  custom_names = NULL,
                  mark_regions = NULL,
                  regions_track_width = 1
){
  
  message("Parsing loci..")
  if(!is.null(loci)){
    chr = as.character(unlist(data.table::tstrsplit(x = loci, spli = ":", keep = 1)))
    start = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, spli = ":", keep = 2)), split = "-"))[1]
    start = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(start))))
    end = unlist(data.table::tstrsplit(x = unlist(data.table::tstrsplit(x = loci, spli = ":", keep = 2)), split = "-"))[2]
    end = as.numeric(as.character(gsub(pattern = ",", replacement = "", x = as.character(end))))
  }else{
    if(any(is.null(chr), is.null(start), is.null(end))){
      stop("Provide a target region to plot. Use argument loci or chr, start, end")
    }
    chr = as.character(chr)
    start = as.numeric(as.character(start))
    end = as.numeric(as.character(end))
  }
  
  if(start >= end){
    stop("End must be larger than Start!")
  }
  message("    Queried region: ", chr, ":", start, "-", end, " [", end-start, " bps]")
  
  if(is.null(bigWigs)){
    stop("Provide at-least one bigWig file")
  }
  message("Checking for files..")
  for(i in 1:length(bigWigs)){
    if(!file.exists(as.character(bigWigs)[i])){
      stop(paste0(as.character(bigWigs)[i], " does not exist!"))
    }
  }
  
  
  .check_bwtool()
  
  if(draw_gene_track){
    .check_mysql()
  }
  
  op_dir = tempdir() #For now
  windows = .gen_windows(chr = chr, start = start, end = end, window_size = binsize, op_dir = op_dir)
  signals = .get_summaries(bedSimple = windows, bigWigs = bigWigs, op_dir = op_dir, nthreads = nthreads)
  .plot_track(
    summary_list = signals,
    chr = chr,
    start = start,
    end = end,
    draw_gene_track = draw_gene_track,
    gene_model = gene_model,
    gtf = isGTF,
    query_ucsc = query_ucsc,
    col = col,
    autoscale = groupAutoScale,
    txname = tx,
    genename = gene,
    gene_width = gene_track_width,
    scale_width = scale_track_width,
    plot_axis = show_axis,
    u_fount = gene_fsize,
    track_names = custom_names,
    regions = mark_regions,
    region_width = regions_track_width, 
    build = build
  )
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

.check_mysql = function(warn = TRUE){
  check = as.character(Sys.which(names = 'mysql'))[1]
  if(check != ""){
    if(warn){
      message("Checking for mysql installation")
      message(paste0("    All good! Found mysql at: ", check))
    }else{
      return(invisible(0))
    }
  }else{
    stop("Could not locate mysql.\nInstall:\n apt install mysql-server\n yum install mysql-server\n conda install -c anaconda mysql")
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
  window_dat[, chr := chr]
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
    x
  })
  
  #Remove intermediate files
  lapply(summaries, function(x) system(command = paste0("rm ", x), intern = TRUE))
  system(command = paste0("rm ", bedSimple), intern = TRUE)

  names(summary_list) = gsub(pattern = "*\\.summary$", replacement = "", x = basename(path = unlist(summaries)))
  summary_list
}

.plot_track = function(summary_list, chr, start, end, draw_gene_track = FALSE, gene_model = NULL, gtf = FALSE, query_ucsc = NULL, build = "hg19", col = "gray70", autoscale = TRUE, txname = NULL, genename = NULL, gene_width = 2, scale_width = 1, plot_axis = TRUE, u_fount = 0.6, track_names = NULL, regions = NULL, region_width = 1){
  
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
    if(length(track_names) != length(summary_list)){
      track_names = rep(x = track_names, length(summary_list))
    }
    names(summary_list) = track_names
  }
  
  if(autoscale){
    plot_height = max(unlist(lapply(summary_list, function(x) max(x$max, na.rm = TRUE))), na.rm = TRUE)
    plot_height = rep(plot_height, length(summary_list))
  }else{
    plot_height = unlist(lapply(summary_list, function(x) max(x$max, na.rm = TRUE)))
  }
  
  plot_height = round(plot_height, digits = 2)
  
  ntracks = length(summary_list)
  if(draw_gene_track){
    if(is.null(gene_model)){
      message("Missing gene model. Trying to query UCSC genome browser..")
      etbl = .extract_geneModel_ucsc(chr = chr, start = start, end = end, refBuild = build, txname = txname, genename = genename)
      if(plot_regions){
        lo = layout(mat = matrix(data = seq_len(ntracks+3)), heights = c(region_width, rep(3, ntracks), gene_width, scale_width))
      }else{
        lo = layout(mat = matrix(data = seq_len(ntracks+2)), heights = c(rep(3, ntracks), gene_width, scale_width))  
      }
    }else{
      if(gtf){
        etbl = .parse_gtf(gtf = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)  
      }else{
        etbl = .extract_geneModel(ucsc_tbl = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)  
      }
      
      if(plot_regions){
        lo = layout(mat = matrix(data = seq_len(ntracks+3)), heights = c(region_width, rep(3, ntracks), gene_width, scale_width))
      }else{
        lo = layout(mat = matrix(data = seq_len(ntracks+2)), heights = c(rep(3, ntracks), gene_width, scale_width))  
      }
    }
  }else{
    if(plot_regions){
      lo = layout(mat = matrix(data = seq_len(ntracks+1)), heights = c(rep(3, ntracks), scale_width))
    }else{
      lo = layout(mat = matrix(data = seq_len(ntracks+1)), heights = c(rep(3, ntracks), scale_width))  
    }
  }
  
  
  #Draw BED regions from `mark_regions`
  if(plot_regions){
    if(plot_axis){
      par(mar = c(0.5, 4, 2, 1))
    }else{
      par(mar = c(0.5, 1, 2, 1))  
    }
    if(nrow(regions) > 0){
      plot(NA, xlim = c(start, end), ylim = c(0, nrow(regions)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      for(i in 1:nrow(regions)){
        rect(xleft = regions[i, startpos], ybottom = i-0.75, xright = regions[i, endpos], ytop = i-0.25, col = "#192a56", border = "#192a56")
        text(x = regions[i, endpos], y = i-0.5, labels = paste0(regions[i, chromsome], ":", regions[i, startpos], "-", regions[i, endpos]), adj = -0.05, xpd = TRUE)
      }
    }else{
      plot(NA, xlim = c(start, end), ylim = c(0, 1), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      message("No overlapping regions!")
    }
    
  }
  
  #Draw bigWig signals
  lapply(1:length(summary_list), function(idx){
    x = summary_list[[idx]]
    if(plot_axis){
      par(mar = c(0.5, 4, 2, 1))
    }else{
      par(mar = c(0.5, 1, 2, 1))  
    }
    
    plot(NA, xlim = c(start, end), ylim = c(0, plot_height[idx]), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
    rect(xleft = x$start, ybottom = 0, xright = x$end, ytop = x$max, col = col[idx], border = col[idx])
    if(plot_axis){
      axis(side = 2, at = c(0, plot_height[idx]), las = 2)  
    }else{
      text(x = start, y = plot_height[idx], labels = paste0("[0-", plot_height[idx], "]"), adj = 0, xpd = TRUE)
    }
    title(main = names(summary_list)[idx], adj = 0, font.main = 3)
  })
  
  #Draw bigWig signals
  if(draw_gene_track){
    if(!is.null(etbl)){
      if(plot_axis){
        par(mar = c(0.5, 4, 0, 0.5))
      }else{
        par(mar = c(0.5, 1, 0, 0.5))  
      }
      
      plot(NA, xlim = c(start, end), ylim = c(0, length(etbl)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      exon_col = "#192a56"
      for(tx_id in 1:length(etbl)){
        txtbl = etbl[[tx_id]]
        segments(x0 = attr(txtbl, "start"), y0 = tx_id-0.45, x1 = attr(txtbl, "end"), y1 = tx_id-0.45, col = exon_col, lwd = 1)
        text(x = start, y = tx_id-0.45, labels = paste0(attr(txtbl, "tx"), " [", attr(txtbl, "gene"), "]"), adj = 0, cex = u_fount)
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
  
  if(plot_axis){
    par(mar = c(3, 4, 0, 1))  
  }else{
    par(mar = c(3, 1, 0, 1))
  }
  lab_at = pretty(c(start, end))
  plot(NA, xlim = c(start, end), ylim = c(0, 1), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  axis(side = 1, at = lab_at, lty = 2, line = 1)
  text(x = end, y = 0.1, labels = paste0(chr, ":", start, "-", end), adj = 1, xpd = TRUE)
}

.extract_geneModel = function(ucsc_tbl = NULL, chr = NULL, start = NULL, end = NULL, txname = txname, genename = genename){
  ucsc = data.table::fread(file = ucsc_tbl)
  query = data.table::data.table(chr, start, end)
  data.table::setkey(x = query, chr, start, end)
  
  colnames(ucsc)[c(3, 5, 6)] = c("chr", "start", "end")
  data.table::setkey(x = ucsc, chr, start, end)
  
  gene_models = data.table::foverlaps(x = query, y = ucsc, type = "any")
  
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

.extract_geneModel_ucsc = function(chr, start = NULL, end = NULL, refBuild = "hg19", txname = NULL, genename = NULL){
  op_file = tempfile(pattern = "ucsc", fileext = ".tsv")
  cmd = paste0("mysql --user genome --host genome-mysql.cse.ucsc.edu -NAD ",  refBuild,  " -e 'select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds from refGene WHERE chrom =\"", chr, "\"'")
  message(paste0("Extracting gene models from UCSC:\n", " chromosome: ", chr, "\n", " build: ", refBuild, "\n query: ", cmd))
  #system(command = cmd)
  ucsc = data.table::fread(cmd = cmd)
  colnames(ucsc) = c("chr", "start", "end", "strand", "name", "name2", "exonStarts", "exonEnds")
  data.table::setkey(x = ucsc, chr, start, end)
  
  query = data.table::data.table(chr, start, end)
  data.table::setkey(x = query, chr, start, end)
  
  gene_models = data.table::foverlaps(x = query, y = ucsc, type = "any")
  
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

.parse_gtf = function(gtf = NULL, chr, start = NULL, end = NULL, refBuild = "hg19", txname = NULL, genename = NULL){
  message("Parsing gtf file..")
  gtf = data.table::fread(file = gtf)
  colnames(gtf) = c("chr", "source", "feature", "start", "end", "ph", "strand", "ph2", "info")
  gtf[,chr := as.character(chr)]
  gtf[,start :=as.numeric(as.character(start))]
  gtf[,end := as.numeric(as.character(end))]
  data.table::setkey(x = gtf, chr, start, end)
  
  query = data.table::data.table(chr, start, end)
  data.table::setkey(x = query, chr, start, end)
  
  gene_models = data.table::foverlaps(x = query, y = gtf, type = "any", nomatch = NULL)
  
  if(nrow(gene_models) == 0){
    warning("No features found within the requested loci!")
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


#Test
# .run_test = function(){
#   bigWigs = list.files(path = "/Volumes/datadrive/bws/trackR/", pattern = "bw", full.names = TRUE)
#   qregion = "chr3:187,439,165-187,463,513"
#   region = data.frame(chr = c("chr3"), start = c(187452998), end = c(187454879))
#   gm = "~/Downloads/hg19.gtf.gz"
#   #with highlight
#   trackplot(bigWigs = bigWigs, loci = qregion, gene_model = gm, draw_gene_track = TRUE, show_axis = FALSE, mark_regions = region, gene = "BCL6", gene_track_width = 2)
# }