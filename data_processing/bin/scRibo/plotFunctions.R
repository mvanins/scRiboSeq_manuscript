#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-04-25
# R functions to return common dataframes and scRibo QC plots

library(dplyr)
library(tidyr)
library(forcats)
 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

reformatFeature <- function(plotlist){
  pltlist <- lapply(plotlist, 
                    function(x) x+
                      coord_equal()+
                      scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")+
                      theme(axis.line = element_blank(),
                            axis.title = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank())
                      )
  return(pltlist)
}

select_top_groups <- function(tbl, n, order_by, with_ties = TRUE){
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist

  keep = tbl %>%
    select(c(grps, {{ order_by }})) %>%
    distinct({{order_by}}) %>% 
    ungroup() %>%
    slice_max(order_by = {{order_by}},
              n=n, with_ties = with_ties) %>%
    select(grps)

  tbl %>% right_join(keep, by=grps) 
}

readCounts <- function(folder, sname = ""){
  require(Matrix)
  fn_mtx <- list.files(folder, pattern = "*matrix.mtx*")
  if(length(fn_mtx) > 1){
    stop(paste0("Too many count matrix *matrix.mtx* files found in folder: ", folder))
  }
  if(length(fn_mtx) < 1){
    stop(paste0("Cannot find count matrix *matrix.mtx* in folder: ", folder))
  }

  fn_features <- list.files(folder, pattern = "*features.tsv*")
  if(length(fn_features) > 1){
    stop(paste0("Too many feature files *features.tsv* found in folder: ", folder))
  }
  if(length(fn_features) < 1){
    stop(paste0("Cannot find feature tsv *features.tsv* in folder: ", folder))
  }

  fn_barcodes <- list.files(folder, pattern = "*barcodes.tsv*")
  if(length(fn_barcodes) > 1){
    stop(paste0("Too many barcode files *barcodes.tsv* found in folder: ", folder))
  }
  if(length(fn_barcodes) < 1){
    stop(paste0("Cannot find barcode tsv *barcodes.tsv* in folder: ", folder))
  }


  counts <- readMM(file.path(folder,fn_mtx))
  rownames(counts) <- read.table(file.path(folder, fn_features), 
                                 header = FALSE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE)$V1
  if(nchar(sname) > 0){
    colnames(counts) <- paste(sname, 
                              read.table(file.path(folder, fn_barcodes),
                                         header = FALSE,
                                         sep = "\t",
                                         stringsAsFactors = FALSE)$V1,
                              sep = "_")
  }else{
    colnames(counts) <- read.table(file.path(folder, fn_barcodes),
                                   header = FALSE,
                                   sep = "\t",
                                   stringsAsFactors = FALSE)$V1
  }

  return(counts)
}

readMergeCounts <- function(samples){
  require(Matrix.utils)
  counts <- mapply(readCounts,
                   folder = samples$counts,
                   sname = samples$names)
  allcounts <- Reduce(function(x,y) {
                        allgenes <- unique(c(rownames(x), rownames(y)))
                        notInX <- allgenes[!(allgenes %in% rownames(x))]
                        notInY <- allgenes[!(allgenes %in% rownames(y))]

                        if(length(notInX > 0)){
                          x <- rbind(x,
                                     Matrix(0, 
                                            nrow = length(notInX),
                                            ncol = ncol(x),
                                            dimnames = list(notInX, colnames(x)),
                                            sparse = TRUE))
                        }
                        if(length(notInY > 0)){
                          y <- rbind(y,
                                     Matrix(0, 
                                            nrow = length(notInY),
                                            ncol = ncol(y),
                                            dimnames = list(notInY, colnames(y)),
                                            sparse = TRUE))
                        }

                        mm <- merge.Matrix(x, y, 
                                           by.x = rownames(x),
                                           by.y = rownames(y),
                                           all.x = TRUE,
                                           all.y = TRUE)
                        return(mm)},
                        counts)
  #rownames(allcounts) <- gsub("_", "-", rownames(allcounts))

  return(allcounts)
}

# readMergeCounts <- function(samples){
# 		require(Matrix.utils)
# 		counts <- mapply(readCounts,
# 						 folder = samples$counts,
# 						 sname = samples$names)
# 		allcounts <- Reduce(function(x,y) {
# 									allgenes <- unique(c(rownames(x), rownames(y)))
# 									notInX <- allgenes[!(allgenes %in% rownames(x))]
# 									notInY <- allgenes[!(allgenes %in% rownames(y))]
# 
# 									if(length(notInX > 0)){
# 
# 										
# 									merge.Matrix(x, y, 
# 												 by.x = rownames(x),
# 												 by.y = rownames(y)) },
# 							counts)
# 		rownames(allcounts) <- gsub("_", "-", rownames(allcounts))
# 
# 		return(allcounts)
# }

readSoloRaw <- function(sname,fname, type = "raw"){
  require(Matrix)
  folder <- file.path(fname, "Gene", type)
  counts <- readMM(file.path(folder,"matrix.mtx.gz"))
  rownames(counts) <- read.table(file.path(folder, "features.tsv.gz"), 
                                 header = FALSE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE)$V1
  colnames(counts) <- paste(sname,
                            read.table(file.path(folder,"barcodes.tsv.gz"),
                                       header = FALSE,
                                       sep = "\t",
                                       stringsAsFactors = FALSE)$V1,
                            sep="_")

  return(counts)
}

readSolo <- function(samples, type = "raw"){
  require(Matrix.utils)
  counts <- mapply(readSoloRaw, sname=samples$names, fname=samples$solo, type=type)
  allcounts <- Reduce(function(x,y) { merge.Matrix(x, y, by.x = rownames(x), by.y = rownames(y)) }, 
                      counts)
  rownames(allcounts) <- gsub("_", "-", rownames(allcounts))

  return(allcounts)
}

columnFeatures <- function(columns,rows){
  # columns: numbers e.g., 1:24
  # rows: letters e.g., LETTERS[1:16]
  return(unlist(lapply(rows,function(x) {paste0(x,columns)})))
}


readReadsDT <- function(sname, fname){
  require(data.table)

  reads <- fread(fname, showProgress = FALSE)

  if("CB" %in% colnames(reads)){
    reads <- na.omit(reads, cols="CB")
    reads <- reads[CB != ""][, CB := paste(sname, CB, sep="_")]
  }
  reads <- reads[, sample := sname]

  return(reads)
}

importReads <- function(samples){
  require(data.table)

  reads <- rbindlist( mapply(readReadsDT, sname = samples$names, fname = samples$reads, SIMPLIFY = FALSE), fill=TRUE)
  return(reads)
}

importReadsFeather <- function(samples){
  require(feather)
  reads = data.frame()
  for(i in 1:nrow(samples)){
    tempreads <- read_feather(samples[i,"reads"])
    tempreads <- tempreads %>%
      filter(!is.na(CB)) %>%
      mutate(CB = paste(samples[i,"names"], CB, sep="_"),
             sample = samples[i,"names"])
      reads <- bind_rows(reads,tempreads)
      rm(tempreads)
  }

  return(reads)
}


readFACS <- function(samples, barcodes = NULL){
		if(!missing(barcodes)){
				rownames(barcodes) <- barcodes$well
		}
		facs <- data.frame()
		for(i in 1:nrow(samples)){
				tempfacs <- read.csv(samples$facs[i], header = TRUE, stringsAsFactors = FALSE)
				tempfacs$sample <- samples$names[i]

				if(!missing(barcodes)){
						tempfacs$CB <- paste(samples$names[i], barcodes[tempfacs$Well,"CB"], sep = "_")
				}
				facs <- bind_rows(facs, tempfacs)
		}
		return(facs)
}

cosine_similarity <- function(DF){
  Matrix <- as.matrix(DF)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  return(sim)
}

## wiggles
wiggleRegion <- function(reads, group = "single", site = 'cut5', reference = "cds_start", mindist = -40, maxdist = 60, name=NULL) {

  region <- reads %>%
    mutate(csd = .data[[site]] - .data[[reference]]) %>%
    filter(csd >= mindist,
           csd <= maxdist) %>%
    { if (startsWith(group,'b')) group_by(., csd) else if (startsWith(group, 's')) group_by(., csd, CB) } %>%
    summarize(n = n() ) %>%
    ungroup() %>%
    { if (startsWith(group,'b')) complete(., csd = full_seq(mindist:maxdist, 1)) else if (startsWith(group,'s')) complete(., csd = full_seq(mindist:maxdist, 1), CB) } %>%
    mutate(n = replace_na(n, 0)) %>%
    {if (is.character(name)) mutate(., name = name) else . }  

  return(region)
}

wiggleDF <- function(reads, group = 'bulk', site = 'cut5', upstream = 50 , downstream = 100){
	# returns data frame for start/stop wiggle 
	start <- wiggleRegion(reads, group = group, site = site, reference = "cds_start",  mindist = -upstream, maxdist = downstream, name = "start")
	stop  <- wiggleRegion(reads, group = group, site = site, reference = "cds_end", mindist = -upstream, maxdist = downstream, name = "stop")
	
	## combine
	wiggle <- bind_rows(start,stop) %>%
			mutate(site = site) %>%
			group_by(name, csd) %>%
			mutate(position_total = sum(n,na.rm=TRUE)) %>%
			ungroup()
	
	if (startsWith(group, 's')){
			wiggle <- wiggle %>%
					group_by(CB) %>%
					mutate(cell_total = sum(n)) %>%
					ungroup() %>%
					mutate(CB = fct_reorder(CB, cell_total),
						   cell_rank = dense_rank(desc(cell_total)))
	}
	return(wiggle)
}


## length metaheatmap
lengthMetaHeatmap <- function(reads, sites = c('cut5', 'cut3'), upstream = 60, downstream = 60, lengthRange = NULL){
	position = 0
	lengthmeta <- data.frame()
	if(is.null(lengthRange)){
			lengthRange <- c(min(reads$length),max(reads$length))
	}
	for(site in sites){
			start <- reads %>%
					mutate(csd = .data[[site]] - cds_start) %>%
					filter(csd > -upstream,
						   csd < downstream) %>%
					group_by(csd, length) %>%
					summarize(n=n()) %>%
					ungroup() %>%
					complete(csd = full_seq(c(-upstream,downstream),1),
						 length = full_seq(lengthRange,1)) %>%
					mutate(csd = replace_na(csd,0),
	   				   length = replace_na(length,0)) %>%
					mutate(set = paste0(site, " start codon"), 
						   pos = position)

			position = position + 1

			stop <- reads %>%
					mutate(csd = .data[[site]] - (cds_end)) %>%
					filter(csd > -downstream,
						   csd < upstream) %>%
					group_by(csd, length) %>%
					summarize(n=n()) %>%
					ungroup() %>%
					complete(csd = full_seq(c(-downstream,upstream),1),
							 length = full_seq(lengthRange,1)) %>%
					mutate(csd = replace_na(csd,0),
						   length = replace_na(length,0)) %>%
					mutate(set = paste0(site, " stop codon"), 
						   pos = position)
					
			position = position + 1
			
			lengthmeta <- bind_rows(lengthmeta, start, stop)
	}

	lengthmeta <- lengthmeta %>%
			mutate(set = fct_reorder(set,pos))

	plt_length_metaheatmap <- ggplot(lengthmeta,aes(x=csd,y=length,fill=n)) +
			geom_tile() +
			facet_wrap(~set, ncol=2) +
			coord_equal() + 
			scale_fill_gradientn(colours=rev(brewer.pal(7,"RdYlBu")),na.value=rev(brewer.pal(7,"RdYlBu"))[1])+
			scale_y_continuous(expand=c(0,0))+
		   	scale_x_continuous(expand=c(0,0))+
			theme( panel.border = element_rect(colour = "black", fill=NA))

	return(plt_length_metaheatmap)
}


# scaled DF
scaleDF <- function(reads, group = 'bulk', site = 'cut5', utr5 = 100, cds = 200, utr3 = 100){
  scaled <- reads %>%
    mutate( pos = case_when( .data[[site]] < l_utr5 ~ round( (.data[[site]]/l_utr5)*(utr5-1) )+1,
                            ((.data[[site]] >= l_utr5) & (.data[[site]] < (l_utr5+l_cds) )) ~ round( ((.data[[site]]-l_utr5)/l_cds)*(cds-1) ) + 1 + utr5,
                            (.data[[site]] >= (l_utr5 + l_cds)) ~ round( ((.data[[site]] - l_utr5 - l_cds)/l_utr3)*(utr3-1) ) + utr5 + cds + 1 ) ) 

  if (group == "single"){
    scaled <- scaled %>%
      group_by(CB,pos)
  }else if (group == "bulk"){
    scaled <- scaled %>%
      group_by(pos)
  }
  
  scaled <- scaled %>%
    summarize(n=n()) %>%
    ungroup()
  
  if (group == "single"){
    scaled <- scaled %>%
      complete(pos = full_seq(c(1,utr5 + cds + utr3),1), CB) %>%
      group_by(CB) %>%
      mutate(tot = sum(n, na.rm=TRUE)) %>%
      ungroup()
  
  } else if (group == "bulk") {
    scaled <- scaled %>%
      complete(pos = full_seq(c(1,utr5+cds+utr3),1)) %>%
      mutate(tot = sum(n, na.rm=TRUE)) 
  }
  
  scaled <- scaled %>%
    mutate(n=replace_na(n,0))
  
  return(scaled)
}

# Percent of reads in frame
frameDF <- function(reads, site = 'cut5'){
  frames <- reads %>%
    mutate(frame = (.data[[site]] - cds_start) %% 3) %>%
    group_by(frame) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    mutate(tot=sum(n),
           set = site)

    return(frames)
}


