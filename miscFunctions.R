#A few handy R-functions
#Author: Simon Forsberg (mostly)

#Get p-value from lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  if(!is.null(f)){
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
  }
  else
    p <- NA
  return(p)
}

#Called internally from the plot.manhattan functions
#Modified from GenABEL
sortmap <- function (chrom, map, delta = 1) 
{
  chnum <- as.numeric(as.factor(chrom))
  ix <- order(chnum, map)
  map <- map[ix]
  off <- c(0, map[1:(length(map) - 1)])
  off <- map - off
  off[which(off <= 0)] <- delta
  cummap <- cumsum(off)
  
  #To be used by sortmap.ranges, called from plot.manhattan_highlRegion
  chnum.off <- c(0, chnum[1:(length(chnum) - 1)])
  chnum.off <- chnum - chnum.off
  newChr <- which(chnum.off == 1)
  pos2cum <- cummap[newChr] #The starting position of every chromosome in the cumulative coordinates
  pos2cum[1] <- 0 #Switch from pos to cumulative pos: pos + pos2cum
  names(pos2cum) <- unique(chrom) 
  
  out <- list()
  out$ix <- ix
  out$cummap <- cummap
  out$chnum <- chnum
  out$pos2cum <- pos2cum 
  out
}

#Internal function called by plot.manhattan_highlRegions
#Modified from GenABEL
sortmap.ranges <- function (ranges, snp.map) 
{
  #ranges - A GRanges object with the ranges to be drawn on the plot
  #map - The sorted cumulative SNP positions returned by sortmap
  chr <- as.character(decode(ranges@seqnames))
  
  #From
  from <- start(ranges)
  from <- from + snp.map$pos2cum[chr]

  #To
  to <- end(ranges)
  to <- to + snp.map$pos2cum[chr]

  #Wrap up
  out <- list()
  out$cummap.from <- from
  out$cummap.to <- to
  out$chr <- chr
  out
}

plot.manhattan <- function(pos, chr, y, log10 = T, col = wes_palette(name = 'Darjeeling1', n = 2), 
                           cex.yaxis = 1.5, cex.xaxis = 2, ylim = NULL, zoom.chr = NULL, zoom.pos = NULL, 
                           padj.xaxis = NA, cex = .5, ...){
  #Partially recycled from GenABELs plot function
  #pos = vector with genomic positions
  #chr = vector with chromosome info matching pos
  #Input has to be sorted
  require(wesanderson)
  
  #Make continuous positions for plotting
  newmap <- sortmap(chr, pos)
  mymap <- newmap$cummap
  
  if(log10)
    y <- -log10(y)
  
  if(is.null(zoom.chr)){
    if(is.null(ylim))
      ylim <- c(0,max(y, na.rm = T))
    xlim = c(min(mymap), max(mymap))
    chrom.num <- as.numeric(as.factor(as.character(chr)))
    chind <- chrom.num%%length(col)
    idxCH <- which(chind == 0)
    
    plot(mymap[idxCH], y[idxCH], xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, col = col[length(col)], pch = 19, cex = cex, xlab = '', ylab = '', ...)
    axis(side = 2, cex.axis = cex.yaxis)
    for (colidx in c(1:(length(col) - 1))) {
      idxCH <- which(chind == colidx)
      points(mymap[idxCH], y[idxCH], col = col[colidx], pch = 19, cex = .5, ylim = ylim, xlim = xlim, ...)
    }
    
    #Axis
    chrom.uniq <- unique(chr)
    chpos <- c()
    for (j in 1:length(chrom.uniq)) 
      chpos[j] <- mean(mymap[chr == chrom.uniq[j]])
    axis(side = 1, at = chpos, labels = chrom.uniq, cex.axis = cex.xaxis, padj = padj.xaxis)
  }
  else{
    if(!(zoom.chr %in% chr))
      stop(paste('No chromosome named', zoom.chr, 'found in the chr vector', chr))

    y <- y[chr == zoom.chr]
    pos <- pos[chr == zoom.chr]
    chr <- chr[chr == zoom.chr]
    if(is.null(zoom.pos))
      zoom.pos <- c(min(pos) - 1, max(pos) + 1)
    
    plot(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], y[pos > zoom.pos[1] & pos < zoom.pos[2]], 
         xaxt = 'n', yaxt = 'n', col = col[1], pch = 19, cex = cex, xlab = '', ylab = '', ...)
    axis(side = 2, cex.axis = cex.yaxis)
    axis(side = 1, cex.axis = cex.xaxis, padj = padj.xaxis)
  }
}

plot.manhattan2 <- function(pos, chr, y, log10 = T, cex.yaxis = 1.5, cex.xaxis = 2, 
                            ylim = NULL, zoom.chr = NULL, zoom.pos = NULL, col = NULL,
                            xlab = '', ylab = '', xaxt = 's', yaxt = 's', cex = .5, ...){
  #Partially recycled from GenABELs plot function
  #Different design than plot.manhattan
  #pos = vector with genomic positions
  #chr = vector with chromosome info matching pos
  #y = vector OR a list with vectors. If a list, the results will be overlayed with different colors 
  #Input has to be sorted
  require(wesanderson)
  
  #Make continuous positions for plotting
  newmap <- sortmap(chr, pos)
  mymap <- newmap$cummap
  
  if(log10){
    if(class(y) == 'list')
      y <- lapply(y, function(x){-log10(x)})
    else
      y <- -log10(y)
  }
  
  #Multiple results to be overlaid?
  if(class(y) == 'list'){
    if(is.null(col))
      col = wes_palette(name = 'Darjeeling1', n = length(y))
    
    y.list <- y[2:length(y)]
    y <- y[[1]]
  }
  else if(is.null(col))
    col = wes_palette(name = 'Darjeeling1', n = 1)
  
  
  if(is.null(zoom.chr)){
    if(is.null(ylim))
      ylim <- c(0,max(y, na.rm = T))
    xlim = c(min(mymap), max(mymap))
    
    #Setup plot with highlighted chrs
    plot(mymap, y, xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, xlab = '', ylab = '', type = 'n')
    chrom.uniq <- unique(chr)
    chpos <- c()
    for (j in 1:length(chrom.uniq)){
      chpos[j] <- mean(mymap[chr == chrom.uniq[j]])
      
      #Draw rectangle for "even" chrs
      if(j %% 2 == 0){
        rect(xleft = min(mymap[chr == chrom.uniq[j]]), ybottom = ylim[1] - 10, xright = max(mymap[chr == chrom.uniq[j]]), ytop = ylim[2] + 10, 
             col = rgb(220,220,220, maxColorValue = 255, alpha = 120), border = NA)
      }
    }
    axis(side = 1, at = chpos, labels = chrom.uniq, cex.axis = cex.xaxis)
    
    #Draw plot
    par(new = T)
    plot(mymap, y, xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, col = col[1], pch = 19, cex = cex, xlab = xlab, ylab = ylab, ...)
    axis(side = 2, cex.axis = cex.yaxis)
    
    #Overlay more things?
    if(exists('y.list')){
      for(j in 1:length(y.list)){
        points(mymap, y.list[[j]], xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, col = col[j+1], pch = 19, cex = cex, xlab = '', ylab = '', ...)
      }
    }
  }
  else{
    if(!(zoom.chr %in% chr))
      stop(paste('No chromosome named', zoom.chr, 'found in the chr vector', chr))
    
    y <- y[chr == zoom.chr]
    pos <- pos[chr == zoom.chr]
    if(is.null(zoom.pos))
      zoom.pos <- c(min(pos) - 1, max(pos) + 1)
    
    #Draw plot
    plot(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], 
         y[pos > zoom.pos[1] & pos < zoom.pos[2]], 
         xaxt = 'n', yaxt = 'n', col = col[1], pch = 19, cex = cex, xlab = '', ylab = '', ylim = ylim, ...)
    
    #Overlay more things?
    if(exists('y.list')){
      for(j in 1:length(y.list)){
        y.overlay <- y.list[[j]][chr == zoom.chr]
        points(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], 
               y.overlay[pos > zoom.pos[1] & pos < zoom.pos[2]], 
               xaxt = 'n', yaxt = 'n', ylim = ylim, col = col[j+1], pch = 19, cex = cex, xlab = '', ylab = '', ...)
      }
    }
    
    if(yaxt != 'n')
      axis(side = 2, cex.axis = cex.yaxis)
    if(xaxt != 'n')
      axis(side = 1, cex.axis = cex.xaxis)
  }
}

plot.manhattan3 <- function(pos, chr, y, log10 = T, cex.yaxis = 1.5, 
                            cex.xaxis = 2, ylim = NULL, zoom.chr = NULL, 
                            zoom.pos = NULL, pch = 19, cex = .5,  
                            highlSNP = NULL, col.highl = 'red', cex.highl = 2, pch.highl = 17, 
                            col = 'gray', padj.xaxis = NA, ...){
  #Partially recycled from GenABELs plot function
  #The idea with this function is to customly color each SNP. For instance, by LD to a focal SNP
  #Input has to be sorted
  
  #pos = vector with genomic positions
  #chr = vector with chromosome info matching pos
  #y = vector with values to put on the y
  require(wesanderson)
  # stopifnot(!is.null(col))
  
  #Make continuous positions for plotting
  newmap <- sortmap(chr, pos)
  mymap <- newmap$cummap
  
  if(log10)
    y <- -log10(y)
  
  # pch <- rep(19, length(pos))
  # cex <- rep(.5, length(pos))
  # col <- rep(col, length(pos))
  # if(!is.null(highlSNP)){ #highlight SNP
  #   pch[highlSNP] <- pch.highl
  #   cex[highlSNP] <- cex.highl
  #   col[highlSNP] <- col.highl
  # }

  if(is.null(zoom.chr)){
    if(is.null(ylim))
      ylim <- c(0,max(y, na.rm = T))
    xlim = c(min(mymap), max(mymap))
    
    #Setup plot with highlighted chrs
    plot(mymap, y, xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, xlab = '', ylab = '', type = 'n')
    chrom.uniq <- unique(chr)
    chpos <- c()
    for (j in 1:length(chrom.uniq)){
      chpos[j] <- mean(mymap[chr == chrom.uniq[j]])
      
      #Draw rectangle for "even" chrs
      if(j %% 2 == 0){
        rect(xleft = min(mymap[chr == chrom.uniq[j]]), ybottom = ylim[1] - 10, xright = max(mymap[chr == chrom.uniq[j]]), ytop = ylim[2] + 10, 
             col = rgb(220,220,220, maxColorValue = 255, alpha = 120), border = NA)
      }
    }
    axis(side = 1, at = chpos, labels = chrom.uniq, cex.axis = cex.xaxis, padj = padj.xaxis)
    
    #Draw plot
    par(new = T)
    plot(mymap, y, xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, xlab = '', ylab = '', pch = pch, cex = cex, col = col, ...) 
    axis(side = 2, cex.axis = cex.yaxis)
    
    #The highlighted SNPs
    points(mymap[highlSNP], y[highlSNP], ylim = ylim, xlim = xlim, 
           col = col.highl, pch = pch.highl, cex = cex.highl)
  }
  else{
    if(!(zoom.chr %in% chr))
      stop(paste('No chromosome named', zoom.chr, 'found in the chr vector', chr))
    
    if(is.numeric(highlSNP))
      highlSNP <- 1:length(y) %in% highlSNP #Convert from indices to logical
    
    # pch <- pch[chr == zoom.chr]
    # cex <- cex[chr == zoom.chr]
    # col <- col[chr == zoom.chr]
    # chr <- chr[chr == zoom.chr]
    y <- y[chr == zoom.chr]
    pos <- pos[chr == zoom.chr]
    highlSNP <- highlSNP[chr == zoom.chr]
    
    if(is.null(zoom.pos))
      zoom.pos <- c(min(pos) - 1, max(pos) + 1)
    
    # plot(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], y[pos > zoom.pos[1] & pos < zoom.pos[2]], xaxt = 'n', yaxt = 'n', 
    #      col = col[pos > zoom.pos[1] & pos < zoom.pos[2]], 
    #      pch = pch[pos > zoom.pos[1] & pos < zoom.pos[2]], 
    #      cex = cex[pos > zoom.pos[1] & pos < zoom.pos[2]], xlab = '', ylab = '', ...)
    plot(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], y[pos > zoom.pos[1] & pos < zoom.pos[2]], 
         xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
         pch = pch, cex = cex, ...)
    axis(side = 2, cex.axis = cex.yaxis)
    axis(side = 1, cex.axis = cex.xaxis)
    
    #The highlighted SNPs
    points(pos[pos > zoom.pos[1] & pos < zoom.pos[2] & highlSNP], 
           y[pos > zoom.pos[1] & pos < zoom.pos[2] & highlSNP], 
           col = col.highl, pch = pch.highl, cex = cex.highl)
    
  }
}


#WARNING: Currently not working for multiple chromosomes. sortmap.ranges needs to be fixed to plot the ranges in the right places
#I believe I've fixed this
plot.manhattan_highlRegions <- function(pos, chr, y, log10 = T, cex.yaxis = 1.5, cex.xaxis = 2, 
                            ylim = NULL, zoom.chr = NULL, zoom.pos = NULL, col = NULL,
                            xlab = '', ylab = '', ranges, ranges.col = 'black', shadeRanges = F, ...){
  #Partially recycled from GenABELs plot function
  #Different design than plot.manhattan
  #pos = vector with genomic positions
  #chr = vector with chromosome info matching pos
  #y = vector OR a list with vectors. If a list, the results will be overlayed with different colors 
  #Input has to be sorted
  require(wesanderson)
  stopifnot(class(ranges) == 'GRanges')
  
  #Make continuous positions for plotting
  #The SNPs
  newmap <- sortmap(chr, pos)
  mymap <- newmap$cummap
  #The ranges
  newmap.ranges <- sortmap.ranges(ranges, newmap) #Fixed I believe

  if(log10){
    if(class(y) == 'list')
      y <- lapply(y, function(x){-log10(x)})
    else
      y <- -log10(y)
  }
  
  #Multiple results to be overlaid?
  if(class(y) == 'list'){
    if(is.null(col))
      col = wes_palette(name = 'Darjeeling1', n = length(y))
    
    y.list <- y[2:length(y)]
    y <- y[[1]]
  }
  else if(is.null(col))
    col = wes_palette(name = 'Darjeeling1', n = 1)
  
  
  if(is.null(zoom.chr)){
    if(is.null(ylim))
      ylim <- c(-2, max(y, na.rm = T))
    xlim = c(min(mymap), max(mymap))
    
    #Setup plot with highlighted chrs
    plot(mymap, y, xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, xlab = '', ylab = '', type = 'n')
    chrom.uniq <- unique(chr)
    chpos <- c()
    for (j in 1:length(chrom.uniq)){
      chpos[j] <- mean(mymap[chr == chrom.uniq[j]])
      
      #Draw rectangle for "even" chrs
      if(j %% 2 == 0){
        rect(xleft = min(mymap[chr == chrom.uniq[j]]), ybottom = ylim[1] - 10, xright = max(mymap[chr == chrom.uniq[j]]), ytop = ylim[2] + 10, 
             col = rgb(220,220,220, maxColorValue = 255, alpha = 120), border = NA)
      }
    }
    axis(side = 1, at = chpos, labels = chrom.uniq, cex.axis = cex.xaxis)
    
    #Draw plot
    par(new = T)
    plot(mymap, y, xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, col = col[1], pch = 19, cex = .5, xlab = xlab, ylab = ylab, ...)
    axis(side = 2, cex.axis = cex.yaxis)
    
    #Overlay more things?
    if(exists('y.list')){
      for(j in 1:length(y.list)){
        points(mymap, y.list[[j]], xaxt = 'n', yaxt = 'n', ylim = ylim, xlim = xlim, col = col[j+1], pch = 19, cex = .5, xlab = '', ylab = '', ...)
      }
    }
    
    #Draw the ranges
    rect(xleft = newmap.ranges$cummap.from, xright = newmap.ranges$cummap.to, ybottom = -2, ytop = -1, col = ranges.col)
    if(shadeRanges)
      rect(xleft = newmap.ranges$cummap.from, xright = newmap.ranges$cummap.to, ybottom = ylim[1] - 10, ytop = ylim[2] + 10, 
           col = rgb(220,220,220, maxColorValue = 255, alpha = 120), border = rgb(220,220,220, maxColorValue = 255, alpha = 120))
  }
  else{
    if(!(zoom.chr %in% chr))
      stop(paste('No chromosome named', zoom.chr, 'found in the chr vector'))
    
    y <- y[chr == zoom.chr]
    pos <- pos[chr == zoom.chr]
    if(is.null(ylim))
      ylim <- c(-2, max(y, na.rm = T))
    if(is.null(zoom.pos))
      zoom.pos <- c(min(pos) - 1, max(pos) + 1)
    
    #Draw plot
    plot(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], 
         y[pos > zoom.pos[1] & pos < zoom.pos[2]], 
         xaxt = 'n', yaxt = 'n', col = col[1], pch = 19, cex = .5, xlab = '', ylab = '', ylim = ylim, ...)
    
    #Overlay more things?
    if(exists('y.list')){
      for(j in 1:length(y.list)){
        y.overlay <- y.list[[j]][chr == zoom.chr]
        points(pos[pos > zoom.pos[1] & pos < zoom.pos[2]], 
               y.overlay[pos > zoom.pos[1] & pos < zoom.pos[2]], 
               xaxt = 'n', yaxt = 'n', ylim = ylim, col = col[j+1], pch = 19, cex = .5, xlab = '', ylab = '', ...)
      }
    }
    axis(side = 2, cex.axis = cex.yaxis)
    axis(side = 1, cex.axis = cex.xaxis)
    
    #Draw the ranges
    from <- start(ranges)
    to <- end(ranges)
    ranges.chr <- as.character(decode(ranges@seqnames))
    if(!(zoom.chr %in% ranges.chr))
      warning(paste('No chromosome named', zoom.chr, 'found in the ranges object'))
    else{
      rect(xleft = from[ranges.chr == zoom.chr], xright = to[ranges.chr == zoom.chr], ybottom = -2, ytop = -1, col = ranges.col)
      if(shadeRanges)
        rect(xleft = from[ranges.chr == zoom.chr], xright = to[ranges.chr == zoom.chr], ybottom = ylim[1] - 10, ytop = ylim[2] + 10, 
             col = rgb(220,220,220, maxColorValue = 255, alpha = 120), border = rgb(220,220,220, maxColorValue = 255, alpha = 120))
    }
  
    
    
  }
}




midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(',.*','',gsub('\\(|\\[|\\)|\\]','', x)))
  upper <- as.numeric(gsub('.*,','',gsub('\\(|\\[|\\)|\\]','', x)))
  return(round(lower+(upper-lower)/2, dp))
}


plot.freqVStime <- function(freqs, freqs.names = NULL, legendPos = 'topleft', cex.legend = 1){
  #A function just to plot allele freq VS timepoint in the hs selection experiment
  if(is.null(freqs.names)){
    freqs.names <- colnames(freqs)
  }
  freq.generation <- numeric(18)
  freq.treatment <- numeric(18)
  freq.generation[grep(pattern = 'G1_', freqs.names)] <- 1
  freq.generation[grep(pattern = 'G11_', freqs.names)] <- 2
  freq.generation[grep(pattern = 'G25_', freqs.names)] <- 3
  freq.treatment[grep(pattern = 'N1|N2|N3', freqs.names)] <- 'control'
  freq.treatment[grep(pattern = 'N4|N5|N6', freqs.names)] <- 'hs'
  
  pal <- wes_palette(name = 'BottleRocket2', 2)
  col <- character(18)
  col[freq.treatment == 'hs'] <- pal[1] 
  col[freq.treatment == 'control'] <- pal[2] 
  
  plot(freq.generation, unlist(freqs), pch = 19, col = col, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', cex = 2)
  axis(side = 1, at = 1:3, labels = c('G1', 'G11', 'G25'), cex.axis = 2)
  axis(side = 2, cex.axis = 1.5)
  mtext(text = 'Allele frequency', side = 2, cex = 2, line = 2.5)
  # abline(a = case1.scan$Intercept[nr], b = case1.scan$generation[nr], col = pal[2], lty = 2, lwd = 2)
  # abline(a = case1.scan$Intercept[nr] + case1.scan$treatment[nr], b = case1.scan$generation[nr] + case1.scan$`generation:treatment`[nr], col = pal[1], lty = 2, lwd = 2)
  legend(legendPos, c('HS', "control"), col = pal, pch = 19, cex = cex.legend)
}

plot.freqVStime_v2 <- function(freqs, freqs.names = NULL, legendPos = 'topleft', cex.points = 3,
                               cex.legend = 1, cex.lab = 2.5, cex.axis = 1.5, pal = wes_palette(name = 'BottleRocket2', 2),
                               lines = F, lwd = 1, legend = T, ...){
  #A function just to plot allele freq VS timepoint in the hs selection experiment
  require(wesanderson)
  if(is.null(freqs.names)){
    freqs.names <- colnames(freqs)
  }
  freq.generation <- numeric()
  freq.treatment <- numeric()
  
  freq.generation[grep(pattern = 'G1_', freqs.names)] <- 1
  freq.generation[grep(pattern = 'G11_', freqs.names)] <- 2
  freq.generation[grep(pattern = 'G25_', freqs.names)] <- 3
  freq.generation[grep(pattern = 'G100_', freqs.names)] <- 4
  
  freq.treatment[grep(pattern = 'N1|N2|N3', freqs.names)] <- 'control'
  freq.treatment[grep(pattern = 'N4|N5|N6', freqs.names)] <- 'hs'
  
  col <- character()
  col[freq.treatment == 'hs'] <- pal[1] 
  col[freq.treatment == 'control'] <- pal[2] 
  
  par(mar = c(5,5,4,2) + .1)
  plot(freq.generation, unlist(freqs), pch = 19, col = col, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', cex = cex.points, ...)
  axis(side = 1, at = 1:4, labels = c('1', '11', '25', '100'), cex.axis = cex.axis)
  mtext(text = 'Generation', side = 1, cex = cex.lab, line = 3)
  axis(side = 2, cex.axis = cex.axis)
  mtext(text = 'Allele frequency', side = 2, cex = cex.lab, line = 3)
  # abline(a = case1.scan$Intercept[nr], b = case1.scan$generation[nr], col = pal[2], lty = 2, lwd = 2)
  # abline(a = case1.scan$Intercept[nr] + case1.scan$treatment[nr], b = case1.scan$generation[nr] + case1.scan$`generation:treatment`[nr], col = pal[1], lty = 2, lwd = 2)
  if(legend)
    legend(legendPos, c('High Sugar', "Control"), col = pal, pch = 19, cex = cex.legend)
  
  if(lines){
    pops <- paste0('N', 1:6)
    for (i in 1:6) {
      idx <- grep(pops[i], freqs.names)
      pop.freqs <- freqs[, idx, with = F]
      pop.gen <- freq.generation[idx]
      pop.col <- unique(col[idx])
      stopifnot(length(pop.col) == 1)
      
      lines(x = pop.gen[order(pop.gen)], y = pop.freqs[, order(pop.gen), with = F], col = pop.col, lwd = lwd)
    }
  }
}

plot.contBoxes <- function(x, y, col = NULL, spacer1 = .01, spacer2 = 3, xlab = 'Simulated h2', ylab = 'Estimated h2', legend.text = NULL, legend.title = NULL, legend.pos = 'topleft', legend.cex = 1, ...){
  if(is.null(col)){
    library(wesanderson)
    cols <- wes_palette('Darjeeling2', ncol(y))
  }
  
  #Set up grouping variable
  tmp <- seq(from = 0, length.out = ncol(y), by = spacer1)
  tmp <- tmp - mean(tmp)
  grouping <- rep(x, ncol(y)) + rep(tmp, each = length(x))
  #Set up between "group" spacing
  box <- boxplot(unlist(y) ~ grouping, xaxt = 'n', yaxt = 'n', col = cols, plot = F)
  tmp <- rep(seq(0, length.out = length(unique(x)), by = 2), each = ncol(y))
  xPos <- rep(1:ncol(y), length(unique(x))) + tmp*spacer2
  #pos for x-axis
  start <- mean(xPos[1:ncol(y)])
  end <- mean(xPos[(length(xPos) - ncol(y) + 1):length(xPos)])
  xlab.pos <- seq(start, end, length.out = length(unique(x)))
  
  #plot
  boxplot(unlist(y) ~ grouping, xaxt = 'n', at = xPos, yaxt = 'n', col = cols, ...)
  grid(10,10)
  par(new = T)
  boxplot(unlist(y) ~ grouping, xaxt = 'n', at = xPos, yaxt = 'n', col = cols, ...)
  axis(side = 1, at = xlab.pos, labels = unique(x), las = 2, cex.axis = 2)
  if(par('usr')[4] < max(x))
    axis(side = 2, cex.axis = 2)
  else
    axis(side = 2, at = unique(x), labels = unique(x), cex.axis = 2)
  mtext(text = xlab, side = 1, line = 4, cex = 2, font = 2)
  mtext(text = ylab, side = 2, line = 2.5, cex = 2, font = 2)
  # lines(x = c(0,25.5), y = c(0, .9), lty = 2, lwd = 3)
  if(!is.null(legend.text))
    legend(legend.pos, legend.text, col = cols, title = legend.title, pch = 15, cex = legend.cex)
}


# sim_pop <- function(N = 200, M = 1000, Fst = 0.1, maf_max = 0.5, maf_min = 0.05, seed = 1){ #Stolen from https://variani.github.io/bigcov/vignettes/popstrat.html
#   set.seed(seed)
#   maf_values <- runif(M, maf_min, maf_max)
#   
#   freq1 <- sapply(1:M, function(i) rbeta(1, 
#                                          maf_values[i] * (1 - Fst) / Fst, 
#                                          (1 - maf_values[i]) * (1 - Fst) / Fst))
#   freq2 <- sapply(1:M, function(i) rbeta(1, 
#                                          maf_values[i] * (1 - Fst) / Fst, 
#                                          (1 - maf_values[i]) * (1 - Fst) / Fst))
#   
#   gdat1 <- sapply(1:M, function(i) sample(c(0, 1, 2), N, replace = TRUE,
#                                           prob = c(((1 - freq1[i])^2), (2 * freq1[i] * (1 - freq1[i])), (freq1[i]^2))))
#   gdat2 <- sapply(1:M, function(i) sample(c(0, 1, 2), N, replace = TRUE,
#                                           prob = c(((1 - freq2[i])^2), (2 * freq2[i] * (1 - freq2[i])), (freq2[i]^2))))
#   
#   gdat <- rbind(gdat1, gdat2)
#   return(gdat)
# }

sim_pop <- function(N = 200, M = 1000, Fst = 0.1, maf_max = 0.5, maf_min = 0.05, seed = 1, nrPops = 2){ #Stolen from https://variani.github.io/bigcov/vignettes/popstrat.html
  #My generalized version
  set.seed(seed)
  maf_values <- runif(M, maf_min, maf_max)
  
  gdat <- matrix(ncol = M, nrow = N*nrPops)
  k <- 1
  for(j in 1:nrPops){
    freq.pop <- sapply(1:M, function(i) rbeta(1, 
                                              maf_values[i] * (1 - Fst) / Fst, 
                                              (1 - maf_values[i]) * (1 - Fst) / Fst))
    
    gdat.pop <- sapply(1:M, function(i) sample(c(0, 1, 2), N, replace = TRUE,
                                               prob = c(((1 - freq.pop[i])^2), (2 * freq.pop[i] * (1 - freq.pop[i])), (freq.pop[i]^2))))
    gdat[k:(k+N-1), ] <- gdat.pop
    k <- k + N
  }
  return(gdat)
}


boxplot.snp.twoWay <- function(marker1.geno, marker2.geno, y, legend = F, names = c('SNP1', 'SNP2'), ...){
  box <- boxplot(y ~ as.matrix(marker1.geno)*as.matrix(marker2.geno), 
                 col = 'lightblue', las = 1, frame=F, cex.axis = .8, xaxt = "n", ...)
  grid()
  par(new = T)
  boxplot(y ~ as.matrix(marker1.geno)*as.matrix(marker2.geno), col = 'lightblue', las = 1, frame=F, cex.axis = .8, xaxt = "n", ...)
  axis(1, at=1:length(box$n), labels=paste("n =", box$n), line=2, lty=0, cex.axis = .8)
  
  labels.marker1 <- gsub(pattern = "(.*)\\..*", replacement = "\\1", x = box$names)
  labels.marker2 <- gsub(pattern = ".*\\.(.*)", replacement = "\\1", x = box$names)
  labels.marker2.nrPerClass <- table(labels.marker2)[1]
  tmp <- 1 + (labels.marker2.nrPerClass - 1)/2
  tmp2 <- sapply(X = 2:length(unique(labels.marker2)) - 1, FUN = function(x){tmp + x*labels.marker2.nrPerClass })
  labels.marker2.at <- c(tmp, tmp2)
  axis(1, at = 1:length(labels.marker1), labels = labels.marker1, cex.axis = .8)
  
  labels.marker2.table <- table(labels.marker2)
  j <- labels.marker2.table[1]
  for(i in 2:length(labels.marker2.table)){
    abline(v = j + .5)
    j <- j + labels.marker2.table[1]
  }
  
  if(legend){
    marker1.uniqGeno <- paste(unique(marker1.geno), collapse=" ")
    marker2.uniqGeno <- paste(unique(marker2.geno), collapse=" ")
    legend("topleft", c(paste(names[2], ":", marker2.uniqGeno, collapse=""), 
                        paste(names[1], ":", marker1.uniqGeno, collapse="")), 
           cex=.8, text.col = c("red", "black"))
    axis(3, at = labels.marker2.at, labels=unique(labels.marker2), cex.axis = .8, col.axis = "red")
  }
  else{
    axis(3, at = labels.marker2.at, labels=unique(labels.marker2), cex.axis = .8)
  }
}


fst <- function(p1, p2){ #Simplistic Fst calculation from two allele freqs. From Noah
  pbar=(p1+p2)/2
  num=((p1**2 + p2**2)/2 - pbar**2)
  den=(pbar*(1-pbar))
  fst=num/den
  return(fst)
}


LDheatmap.hacked <- function (gdat, genetic.distances = NULL, distances = "physical", 
                              LDmeasure = "r", title = "Pairwise LD", add.map = TRUE, add.key = TRUE, 
                              geneMapLocation = 0.15, geneMapLabelX = NULL, geneMapLabelY = NULL, 
                              SNP.name = NULL, color = NULL, newpage = TRUE, name = "ldheatmap", 
                              vp.name = NULL, pop = FALSE, flip = NULL, text = FALSE) 
{
  requireNamespace("grid")
  if (is.null(color)) {
    if (inherits(gdat, "LDheatmap")) 
      color <- gdat$color
    else color <- grey.colors(20)
  }
  if (is.null(flip)) {
    if (inherits(gdat, "LDheatmap") && !is.null(gdat$flipVP)) 
      flip <- TRUE
    else flip <- FALSE
  }
  if (is.null(genetic.distances)) {
    if (inherits(gdat, "data.frame")) 
      genetic.distances = 1000 * (1:ncol(gdat))
    else if (inherits(gdat, "matrix")) 
      genetic.distances = 1000 * (1:length(gdat[1, ]))
    else genetic.distances = gdat$genetic.distances
  }
  if (inherits(gdat, "SnpMatrix")) {
    if (!is.vector(genetic.distances)) {
      stop("Distance should be in the form of a vector")
    }
    o <- order(genetic.distances)
    genetic.distances <- genetic.distances[o]
    gdat <- gdat[, o]
    if (LDmeasure == "r") 
      LDmatrix <- snpStats::ld(gdat, depth = ncol(gdat) - 
                                 1, stats = "R.squared")
    else if (LDmeasure == "D") 
      LDmatrix <- snpStats::ld(gdat, depth = ncol(gdat) - 
                                 1, stats = "D.prime")
    else stop("Invalid LD measurement, choose r or D'.")
    LDmatrix <- as.matrix(LDmatrix)
    LDmatrix[lower.tri(LDmatrix, diag = TRUE)] <- NA
  }
  else if (inherits(gdat, "data.frame")) {
    for (i in 1:ncol(gdat)) {
      if (!genetics::is.genotype(gdat[, i])) 
        stop("column ", i, " is not a genotype object\n")
    }
    gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 
                             2))
    genetic.distances <- genetic.distances[gvars]
    gdat <- gdat[gvars]
    if (!is.vector(genetic.distances)) {
      stop("Distance should be in the form of a vector")
    }
    o <- order(genetic.distances)
    genetic.distances <- genetic.distances[o]
    gdat <- gdat[, o]
    myLD <- genetics::LD(gdat)
    if (LDmeasure == "r") 
      LDmatrix <- myLD[[LDmeasure]]^2
    else if (LDmeasure == "D'") 
      LDmatrix <- abs(myLD[[LDmeasure]])
    else stop("Invalid LD measurement, choose r or D'.")
  }
  else if (inherits(gdat, "LDheatmap")) {
    LDmatrix <- gdat$LDmatrix
    distances <- gdat$distances
  }
  else if (inherits(gdat, "matrix")) {
    if (nrow(gdat) != ncol(gdat)) 
      stop("The matrix of linkage disequilibrium measurements must be a square matrix")
    LDmatrix <- gdat
    LDmatrix[lower.tri(LDmatrix, diag = TRUE)] <- NA
  }
  else if (!missing(gdat)) 
    stop(paste("No method for an object of class", class(gdat)))
  else stop("Need to supply LD matrix or genotypes")
  heatmapVP <- viewport(width = unit(0.8, "snpc"), height = unit(0.8, 
                                                                 "snpc"), name = vp.name)
  flipVP <- viewport(width = unit(0.8, "snpc"), height = unit(0.8, 
                                                              "snpc"), y = 0.6, angle = -45, name = "flipVP")
  if (color[1] == "blueToRed") 
    color = rainbow(20, start = 4/6, end = 0, s = 0.7)[20:1]
  if (newpage) 
    grid.newpage()
  mybreak <- 0:length(color)/length(color)
  imgLDmatrix <- LDmatrix
  byrow <- ifelse(flip, FALSE, TRUE)
  colcut <- as.character(cut(1 - imgLDmatrix, mybreak, labels = as.character(color), 
                             include.lowest = TRUE))
  if (is.numeric(color)) 
    colcut <- as.integer(colcut)
  ImageRect <- LDheatmap:::makeImageRect(dim(LDmatrix)[1], dim(LDmatrix)[2], 
                             colcut, name = "heatmap", byrow)
  ImageText <- NULL
  if (text) 
    ImageText <- makeImageText(dim(LDmatrix)[1], dim(LDmatrix)[2], 
                               round(imgLDmatrix, digits = 2), name = "heatmaptext")
  title <- textGrob(title, 0.5, 1.05, gp = gpar(cex = 1), name = "title")
  if (flip) {
    ImageRect <- editGrob(ImageRect, vp = flipVP)
    if (text) {
      ImageText <- makeImageText(dim(LDmatrix)[1], dim(LDmatrix)[2], 
                                 round(imgLDmatrix, digits = 2), name = "heatmaptext", 
                                 flip = TRUE)
      textVal <- ImageText
      ImageText <- editGrob(ImageText, vp = flipVP, rot = 0, 
                            just = c("right", "top"))
    }
  }
  heatMap <- gTree(children = gList(ImageRect, ImageText, title), 
                   name = "heatMap")
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps - 1)
  # ind <- match(SNP.name, row.names(LDmatrix), nomatch = 0)
  ind <- 1:nrow(LDmatrix)
  geneMapVP <- NULL
  if (flip) 
    geneMapVP <- flipVP
  geneMap <- LDheatmap:::LDheatmapMapNew.add(nsnps, genetic.distances = genetic.distances, 
                                             geneMapLocation = geneMapLocation, add.map, geneMapLabelX = geneMapLabelX, 
                                             geneMapLabelY = geneMapLabelY, distances = distances, 
                                             vp = geneMapVP, SNP.name = SNP.name, ind = ind, flip = flip)
  if (add.key) 
    Key <- LDheatmap:::LDheatmapLegend.add(color, LDmeasure, heatmapVP)
  else Key <- NULL
  LDheatmapGrob <- gTree(children = gList(heatMap, geneMap, 
                                          Key), vp = heatmapVP, name = name, cl = "ldheatmap")
  grid.draw(LDheatmapGrob)
  if (pop) {
    downViewport(heatmapVP$name)
    popViewport()
  }
  ldheatmap <- list(LDmatrix = LDmatrix, LDheatmapGrob = LDheatmapGrob, 
                    heatmapVP = heatmapVP, flipVP = geneMapVP, genetic.distances = genetic.distances, 
                    distances = distances, color = color)
  class(ldheatmap) <- "LDheatmap"
  invisible(ldheatmap)
}

matrixPairsWhere <- function(x, expr, rNames = NULL, cNames = NULL, symmetric = F, ignoreDiag = T){
  #The function retrieve elements in x that fulfills expr and returns 'long' format like:
  #rowname  colname  value
  #   .       .      .
  #   .       .      .
  #Useful for instance when examining cor and distance matrices
  
  # x = A matrix
  # expr = An expression given as a string. Use x to refer to the matrix. Example: 'x < .5' or 'x > 0 & x < 2'
  if(is.null(rNames))
    rNames <- rownames(x)
  if(is.null(cNames))
    cNames <- colnames(x)
  
  if(symmetric){ #Intended for cor/distance matrices etc that are symmetric
    stopifnot(nrow(x) == ncol(x))
    x[lower.tri(x)] <- NA
    if(ignoreDiag)
      diag(x) <- NA #Skip diagonal. Self identity in the case of cor matrix
  }
  
  pairs <- which(eval(parse(text = expr))) #Evaluates the expression
  
  #Convert every idx to a pair of col/row idx
  row.idx <- pairs %% nrow(x)
  row.idx[row.idx == 0] <- nrow(x) #Elements on the last row
  col.idx <- pairs/nrow(x)
  col.idx[col.idx %% 1 != 0] <- floor(col.idx[col.idx %% 1 != 0] + 1) #Ignore elements on the last row
  
  #Wrap up
  result <- data.frame(rowname = rNames[row.idx], colname = cNames[col.idx], val = x[pairs], stringsAsFactors = F)
  # if(rmDupl){ #Useful to remove duplicates when x is symmetric
  #   names1 <- paste(result$rowname, result$colname, sep = '_')
  #   names2 <- paste(result$colname, result$rowname, sep = '_')
  #   dupl <- duplicated(names1)
  #   duplFlip <- duplicated(data.frame(names1, names2))
  #   
  #   if(rmSelfPair)
  #     return()
  #   else
  #     return(result[!(dupl | duplFlip)])
  # }
  # else
    return(result)
}

matrix2pairs <- function (m, samples) {
  #Extract pairwise values from a matrix and return in long format:
  #rowname  colname  value
  #   .       .      .
  #   .       .      .
  #Stolen from extract.val in the phylip package
  
  i <- cbind(match(samples[, 1], rownames(m)), match(samples[, 2], colnames(m)))
  m[i]
}

ll <- function(sort = F, units = 'auto'){
  if(sort)
    sort(sapply(ls(envir=.GlobalEnv), function(x){format(object.size(get(x)), units = units)}), decreasing = T)
  else
    sapply(ls(envir=.GlobalEnv), function(x){format(object.size(get(x)), units = units)})
}

plot.wgcna_softThreshold <- function(sft, cex.text = 0.9){
  #Code recycled from the WGCNA tutorials
  par(mfrow = c(1,2))
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
       main = "Scale independence")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=sft$fitIndices$Power, cex=cex.text, col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=sft$fitIndices$Power, cex=cex.text,col="red")
  
}

add.alpha <- function(col, alpha=1){
  #Add alpha (transparency) to a color vector
  #col = vector with colors (in a form accepted by col2rgb)
  #alpha = alpha value between 0 & 1. Can be one number or a vector matching col
  
  if(missing(col))
    stop("Please provide a vector of colors")
  col.rgb <- rbind(sapply(col, col2rgb)/255, alpha)
  apply(col.rgb, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=x[4]))  
}

plotGenoContTable <- function(geno, panelText = "frac", ...){
  if(ncol(geno) != 2)
    stop("geno should be a genotype matrix with two columns/markers")
  if(!require(lattice))
    stop("Lattice package not found")
  if(!panelText %in% c("frac", "counts"))
    stop("panelText needs to be \"frac\" or \"counts\"")
  if(class(geno) == "matrix")
    geno <- as.data.frame(geno)
  
  geno.table <- table(geno)
  #Expected table
  geno.table_margin1 <- margin.table(geno.table, margin = 1)
  geno.table_margin2 <- margin.table(geno.table, margin = 2)
  geno.table_exp <- geno.table_margin1 %*% t(geno.table_margin2) / sum(geno.table)
  
  if(panelText == "counts"){
    myPanel <- function(x, y, z, ...) { #stolen and modified from stack overflow: http://stackoverflow.com/questions/22827677/entering-cell-values-from-a-matrix-into-a-levelplot-made-in-lattice-in-r
      panel.levelplot(x,y,z,...)
      #Write the observed counts
      panel.text(as.numeric(x) - .1, as.numeric(y) + .1,  "obs:")
      panel.text(as.numeric(x) + .1, as.numeric(y) + .1,  t(geno.table)[cbind(x,y)])
      
      #Write the expected counts
      panel.text(as.numeric(x) - .1, as.numeric(y) - .1,  "exp:")
      panel.text(as.numeric(x) + .1, as.numeric(y) - .1,  round(t(geno.table_exp))[cbind(x,y)])
    }
  }
  else{
    myPanel <- function(x, y, z, ...) { #stolen and modified from stack overflow: http://stackoverflow.com/questions/22827677/entering-cell-values-from-a-matrix-into-a-levelplot-made-in-lattice-in-r
      panel.levelplot(x,y,z,...)
      panel.text(x, y,  format(t(geno.table)/t(geno.table_exp), digits = 2)[cbind(x,y)]) ## use handy matrix indexing
    }
  }
  print(levelplot(t(geno.table), panel = myPanel, ...))
}

#Use to color stuff by a numeric variable
num2color <- function(x, cols = NULL){
  if(is.null(cols)){
    library(RColorBrewer)
    cols <- brewer.pal(4, "Blues")
  }
  
  # Define colour palette
  pal <- colorRampPalette(cols)
  # Rank variable for colour assignment
  order <- findInterval(x, sort(x))
  # Create colors
  pal(length(x))[order]
}

#Mats function to plot a gene based on a GRanges object
plot_gtf <- function(gtf_gr, outfile = NULL, display_id  = F, x_lim = NULL, y_lim = NULL){
  if(!is.null(outfile)){
    pdf(file = outfile, width = 8, height = 4)
  }
  level_count <- max(rowSums(table(gtf_gr$gene_id, gtf_gr$transcript_id) > 0)) + 2
  if(is.null(x_lim)) x_lim <- c(min(start(gtf_gr)), max(end(gtf_gr)))
  if(is.null(y_lim)) y_lim <- c(0, level_count)
  plot(x = mean(start(gtf_gr)), y = level_count /2, xlim = x_lim, ylim = y_lim, type = "n", main = paste("Chr", unique(gtf_gr@seqnames)), xlab = "Position", ylab = "")
  genes <-  unique(gtf_gr$gene_id)
  for (gene in genes){
    current_lvl <- level_count
    gene_entry <- gtf_gr$type == "gene" & gtf_gr$gene_id == gene
    
    text(x = mean(c(start(gtf_gr)[gene_entry], end(gtf_gr)[gene_entry])), y = current_lvl - 0.6, labels = gene, cex = 0.8)
    current_lvl <- current_lvl - 1
    
    segments(x0 = start(gtf_gr)[gene_entry], x1 = end(gtf_gr)[gene_entry], y0 = current_lvl, col = "firebrick", lwd = 5)
    transcript_entries <- which(gtf_gr$type == "transcript" & gtf_gr$gene_id == gene)
    current_lvl <- current_lvl - 1
    
    for(transcript in transcript_entries){
      transcript_id <- gtf_gr$transcript_id[transcript]
      if(display_id){
        text(x = start(gtf_gr)[transcript], y = current_lvl + 0.4, labels = transcript_id, cex = 0.8, pos = 4)
      }
      segments(x0 = start(gtf_gr)[transcript], x1 = end(gtf_gr)[transcript], y0 = current_lvl, col = "black", lwd = 2)
      exon_entries <- gtf_gr$type == "CDS" & gtf_gr$transcript_id == transcript_id
      rect(xleft = start(gtf_gr)[exon_entries], xright = end(gtf_gr)[exon_entries ], ybottom = current_lvl - 0.2, ytop = current_lvl + 0.2, col = "black")
      if(any(gtf_gr$type == "five_prime_utr" & gtf_gr$transcript_id == transcript_id)){
        utr5_entries <- gtf_gr$type == "five_prime_utr" & gtf_gr$transcript_id == transcript_id
        rect(xleft = start(gtf_gr)[utr5_entries], xright = end(gtf_gr)[utr5_entries], ybottom = current_lvl - 0.2, ytop = current_lvl + 0.2, col = "firebrick1")
      }
      if(any(gtf_gr$type == "three_prime_utr" & gtf_gr$transcript_id == transcript_id)){
        utr3_entries <- gtf_gr$type == "three_prime_utr" & gtf_gr$transcript_id == transcript_id
        rect(xleft = start(gtf_gr)[utr3_entries], xright = end(gtf_gr)[utr3_entries], ybottom = current_lvl - 0.2, ytop = current_lvl + 0.2, col = "firebrick4")
      }
      
      current_lvl <- current_lvl - 1
    }
  }
  if(!is.null(outfile)){
    dev.off()
  }
}

#Mats function
plot_cnv_coverage <- function(cnv_list, pdf_file = "~/Projects/Herring/doc/cnv_v2.0.2/cnv_coverage.pdf", size_df  = Ch_v2.0.2_sizes){
  del_GR_list <- cnv_list$del
  dup_GR_list <- cnv_list$dup
  del_cov_v2.0.2 <- coverage(unlist(del_GR_list))
  del_cov_v2.0.2_df <- as.data.frame(ranges(del_cov_v2.0.2))
  del_cov_v2.0.2_df[,"cov"] <- unlist(runValue(del_cov_v2.0.2))
  del_cov_v2.0.2_df[,"global_start"] <- del_cov_v2.0.2_df[,"start"] + size_df[match(del_cov_v2.0.2_df[, "group_name"], size_df[,"name"]), "offset"]
  del_cov_v2.0.2_df[,"global_end"] <-  del_cov_v2.0.2_df[,"end"] + size_df[match(del_cov_v2.0.2_df[, "group_name"], size_df[,"name"]), "offset"]
  
  dup_cov_v2.0.2 <- coverage(unlist(dup_GR_list))
  dup_cov_v2.0.2_df <- as.data.frame(ranges(dup_cov_v2.0.2))
  dup_cov_v2.0.2_df[,"cov"] <- unlist(runValue(dup_cov_v2.0.2))
  dup_cov_v2.0.2_df[,"global_start"] <-  dup_cov_v2.0.2_df[,"start"] + size_df[match(dup_cov_v2.0.2_df[, "group_name"], size_df[,"name"]), "offset"]
  dup_cov_v2.0.2_df[,"global_end"] <-  dup_cov_v2.0.2_df[,"end"] + size_df[match(dup_cov_v2.0.2_df[, "group_name"], size_df[,"name"]), "offset"]
  
  pdf(pdf_file , width = 15, height = 6)
  par(mar = c(5, 7, 4, 2) + 0.1)
  plot(x=0, y= 0, xlim = c(0, max(del_cov_v2.0.2_df[,"global_end"])), ylim = c(0,80), type = "n", xlab = "", ylab = "Number of individuals", axes = F, main = "Whole genome")
  axis(2)
  axis(1, labels = F, at = size_df[1:27,"offset"], pos = -1)
  rect(xleft = size_df[1:26,"offset"], xright = size_df[2:27,"offset"], ybottom = 0, ytop = 80, col = c("grey80", "white"), border = NA)
  
  segment_filter <- grepl("chr", dup_cov_v2.0.2_df[, "group_name"]) # dup_cov_v2.0.2_df[, "cov"] > 0 #Optional
  segments(x0 = dup_cov_v2.0.2_df[segment_filter, "global_start"], x1 = dup_cov_v2.0.2_df[segment_filter, "global_end"], y0 = dup_cov_v2.0.2_df[segment_filter, "cov"], lwd = 1.5, col = "steelblue")
  segment_filter <- grepl("chr", del_cov_v2.0.2_df[, "group_name"]) # del_cov_v2.0.2_df[, "cov"] > 0 #Optional
  segments(x0 = del_cov_v2.0.2_df[segment_filter, "global_start"], x1 = del_cov_v2.0.2_df[segment_filter, "global_end"], y0 = del_cov_v2.0.2_df[segment_filter, "cov"], lwd = 1.5, col = "firebrick")
  
  
  mtext(paste("Chr", 1:26), 1, 0, at = rowMeans(cbind(size_df[1:26,"offset"],size_df[2:27, "offset"])), las = 2, cex = 1.6)
  
  for(chr in  size_df[1:26,"name"]){
    plot(x=0, y= 0, xlim = c(0,size_df[size_df[,"name"] == chr, "size"]), ylim = c(0,80), type = "n", ylab = "Number of individuals", xlab = "Position", main = chr)
    segment_filter <- dup_cov_v2.0.2_df[, "group_name"] == chr
    segments(x0 = dup_cov_v2.0.2_df[segment_filter, "start"], x1 = dup_cov_v2.0.2_df[segment_filter, "end"], y0 = dup_cov_v2.0.2_df[segment_filter, "cov"], lwd = 3, col = "steelblue")
    
    segment_filter <- del_cov_v2.0.2_df[, "group_name"] == chr # del_cov_v2.0.2_df[, "cov"] > 0 #Optional
    segments(x0 = del_cov_v2.0.2_df[segment_filter, "start"], x1 = del_cov_v2.0.2_df[segment_filter, "end"], y0 = del_cov_v2.0.2_df[segment_filter, "cov"], lwd = 3, col = "firebrick")
  }
  
  dev.off()
  return(invisible(list(del_df = del_cov_v2.0.2_df, dup_df = dup_cov_v2.0.2_df)))
}


# Function to plot color bar. Stolen and slightly modified from https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), cex.axis = 2, ...) {
  scale = (length(lut)-1)/(max-min)
  
  par(mar = c(1,5,4,1))
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', ...)
  axis(2, ticks, las=1, cex.axis = cex.axis, line = -1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


#Function to simulate a pair of polygenic traits. With pleiotropy etc
genBivarCov.simulatePair_outbred <- function(n_loci = 200, n_ind = 250, n_perInd = 4, 
                                             h2.add.t1 = 0.4, h2.add.t2 = 0.4, # The heritabilities of the two traits. If rhog > 0, h2.add.t2 will approach h2.add.t1
                                             h2.cov = 0.4,                     # The heritability of the intra-individual correlation/causality t1 -> t2
                                             rhog = 0,                         # Genetic correlation (pleiotropy). Here, it is defined as the fraction of overlapping causal loci. Hence rhog is proportional, but not equal to, genetic correlation
                                             mu1 = 0, mu2 = 0,                 # The intercepts/population means for the traits
                                             beta_12 = 0,                      # The population mean of the effect t1 -> t2
                                             geno = NULL,                      # If null, the genotype is simulated
                                             freqDistr = 'U',                  # Distribution to draw allele frequencies from 
                                             returnGRM = T, 
                                             returnGeno = F,
                                             nrPops = NULL, Fst = NULL,        # Controls the degree of genetic stratification in the simulated population
                                             postproc = FALSE,
                                             e.cov = 0,
                                             trait1.binary = F, 
                                             trait1.fixed = F,                 # trait1 the same across all individuals. For instance equal time points when trait2 is measured. This should affect var(trait2)
                                             qtls.t1 = NA, qtls.t2 = NA        # Number of large-ish QTLs to add per trait (above the polygenic effect). NOTE: h2.add.t1 and h2.add.t2 will no longer be correct
                                             )                 
{
  #Same scenario as in the function genBivarCov.simulatePairKin() but simulating genotypes in outbred populations, rather than nuclear families
  #The GRM is SNP based and hence more dense than in the *Kin functions
  #NOTE: rhog is the fraction of overlapping causal loci, which is proportional but not equal to the genetic correlation. When rhog > 0, h2.add.t2 and h2.cov should also be proportional but not equal to the true parameter values
  # Simulation:
  #   trait1 ~ g.add.t1 + e.t1
  #   trait2 ~ g.add.t2 + g.cov * trait1 + e.t2
  # Where:
  #   g.cov = polygenic effect on cor(trait1, trait2)
  #   g.add = additive polygenic effect on trait1/trait2. g.add per individual is the sum of all the SNP effects, where the fraction of overlapping causal SNPs is determined by rhog
  
  stopifnot(h2.add.t1 >= 0 & h2.add.t1 <= 1)
  totVar.t2 <- h2.cov + h2.add.t2 + e.cov
  stopifnot(totVar.t2 >= 0 & totVar.t2 <= 1)
  
  #### Simulate Genotypes ####
  if(is.null(geno)){
    if(!is.null(nrPops) & !is.null(Fst)){ #Simulate stratified population
      stopifnot(!n_ind %% nrPops)
      geno <- sim_pop(N = n_ind/nrPops, M = n_loci, Fst = Fst, nrPops = nrPops)
    }
    else{ #Simulate unrelated individuals
      #Allele frequency distribution
      if(freqDistr == 'U' | freqDistr == 'u'){
        #Draw alleles from from a U-shaped freq distribution
        tmp <- seq(from = .01, to = .99, by = .01)
        prob <- 1/(tmp*(1-tmp))
        prob <- prob/max(prob)
        p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
      }
      else if(freqDistr == 'uni' | freqDistr == 'uniform')
        p <- runif(n_loci, 0.1, 0.9)
      else if(freqDistr > 0 & freqDistr < 1)
        p <- rep(freqDistr, n_loci)
      else
        stop(paste('Did not recognize freqDistr:', freqDistr))
      
      #Assign genotypes according to HW
      geno <- matrix(nrow = n_ind, ncol = n_loci)
      for(i in 1:n_ind){
        # geno[i,] <- sapply(X = p, FUN = function(x){sample(x = -1:1, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))}) #Same thing as below but slower
        haplo1 <- as.numeric(runif(n_loci) < p)
        haplo2 <- as.numeric(runif(n_loci) < p)
        geno[i,] <- haplo1 + haplo2
      }
    }
    
    #Remove non-polymorphic loci
    fixed <- apply(geno, 2, var) == 0
    if(any(fixed)){
      geno <- geno[, !fixed]
      n_loci <- sum(!fixed)
    }
  }
  
  # Center & scale
  col_means <- colMeans(geno, na.rm = TRUE)
  col_sd <- apply(geno, 2, sd, na.rm = TRUE)
  # col_freq <- col_means / 2  # col_means = 2 * col_freq
  # col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  Z <- sweep(geno, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")
  
  if(returnGRM){
    Zg <- Z / sqrt(n_loci)
    G <- tcrossprod(Zg) # the same as tcrossprod(Z) / M
    rownames(G) <- colnames(G) <- 1:n_ind 
  }
  else
    G <- NA
  
  
  #### Simulate phenotypes ####
  dat <- data.frame(id = rep(c(1:n_ind), each = n_perInd),
                    rep = rep(1:n_perInd, n_ind), obs = seq(1, n_ind*n_perInd),
                    trait1 = 0.0, trait2 = 0.0,
                    g.add.t1 = 0.0, g.add.t2 = 0.0,
                    g.cov = 0.0)
  
  #Additive genetic effect per locus. Amount of overlap in causative loci determined by rhog. This is what induces the pleiotropy / genetic correlation
  b.trait1 <- rnorm(n_loci, 0, sqrt(h2.add.t1/n_loci))
  b.trait2 <- rnorm(n_loci, 0, sqrt(h2.add.t2/n_loci))
  #Add larger QTLs. NOTE: h2.add.t1 and h2.add.t2 will no longer be correct
  if(!is.na(qtls.t1)){
    qtls.t1 <- sample(1:n_loci, qtls.t1)
    b.trait1[qtls.t1] <- rnorm(qtls.t1, 1, .5)
  }
  if(!is.na(qtls.t2)){
    qtls.t2 <- sample(1:n_loci, qtls.t2)
    b.trait2[qtls.t2] <- rnorm(qtls.t2, 1, .5)
  }
  
  if(rhog > 0)
    b.trait2[1:round(rhog*n_loci)] <- b.trait1[1:round(rhog*n_loci)]
  
  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.add.t1 <- Z %*% b.trait1
  g.add.t2 <- Z %*% b.trait2
  dat$g.add.t1 <- g.add.t1[dat$id]
  dat$g.add.t2 <- g.add.t2[dat$id]
  
  #Genetic effect per locus on individual correlation
  b.cov <- rnorm(n_loci, 0, sqrt(h2.cov/n_loci)) 
  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.cov <- Z %*% b.cov
  dat$g.cov <- g.cov[dat$id]
  
  #Non genetic effect on individual slope
  dat$res.cov <- rep(rnorm(n_ind*n_perInd, sd = sqrt(e.cov)))
  
  #Individual phenotypes are the sum of the genetic effects + environmental noise
  if(trait1.fixed){
    t <- 1:n_perInd - mean(1:n_perInd)
    dat$trait1 <- rep(t, times = n_ind)
  }
  else{
    y1 <- mu1 + dat$g.add.t1 + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.add.t1))
    if(trait1.binary)
      dat$trait1 <- as.numeric(cut(y1, breaks = 2)) - 1 #Discretize trait1
    else
      dat$trait1 <- y1
  }
  
  dat$trait2 <- mu2 + dat$g.add.t2 + (beta_12 + dat$g.cov + dat$res.cov) * dat$trait1 + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.add.t2 - h2.cov - e.cov))
  
  if(postproc) {
    dat <- within(dat, {
      id <- as.character(id)
      rid <- id
      #trait1 <- as.numeric(scale(trait1))
      #trait2 <- as.numeric(scale(trait2))
    })  
  }
  
  if(!returnGeno)
    geno <- NULL
  
  return(list(pheno = dat, G = G, geno = geno, qtls.t1 = qtls.t1, qtls.t2 = qtls.t2))
}


#Modification of the clump.markers function from cgmisc. I've modified it to:
# - Take data.frame input rather than GenABEL objects
# - Take arbitrary metrics to base the clumping on, instead of just r2
#Authors: Marcin Kierczak, Simon Forsberg
clump.markers <- function (assoc, snp.cor, chr = 1, bp.dist = 250000, p1 = 1e-04, 
                           p2 = 0.01, r2 = 0.5, image = F, verbose = F) 
{
  #assoc = a data.frame (or table) with p-vals per SNP (pos, chr, p)
  #snp.cor = snp-by-snp matrix with metric to base clumping on. r^2 in the plink original 
  if(!all(colnames(assoc)[1:3] == c('pos', 'chr', 'p')))
    stop('The first three columns of assoc need to be: pos, chr, p')
  if(ncol(snp.cor) != nrow(snp.cor))
    stop('snp.cor is not square. It should be a snp-x-snp matrix with metric to clump on')
  if(nrow(snp.cor) != nrow(assoc))
    stop('Dimensions of assoc and snp.cor do not match')
  
  snpNames <- paste(assoc$pos, assoc$chr, sep = '_')
  duplSNPs <- duplicated(snpNames)
  if(any(duplSNPs)){
    warning('Data contains duplicated SNPs. Removing duplicates')
    snpNames <- snpNames[!duplSNPs]
    assoc <- assoc[!duplSNPs, ]
    snp.cor <- snp.cor[!duplSNPs, !duplSNPs]
  }
  rownames(assoc) <- snpNames
  rownames(snp.cor) <- colnames(snp.cor) <- snpNames
  
  assoc.chr <- assoc[assoc$chr == chr,]
  assoc.sorted <- assoc[order(assoc.chr$p), ]
  signif.p1 <- rownames(assoc.sorted)[assoc.sorted$p <= p1]
  signif.p2 <- rownames(assoc.chr)[assoc.chr$p <= p2]
  assoc.signif <- assoc.chr[rownames(assoc.chr) %in% signif.p2, ]
  snp.cor <- snp.cor[signif.p2, signif.p2]
  
  d <- as.matrix(dist(cbind(assoc.signif$pos, rep(0, times = length(signif.p2)))))
  clumpmatrix <- matrix(rep(0, times = length(signif.p2)^2), 
                        nrow = length(signif.p2), ncol = length(signif.p2))
  clumpmatrix[which(d <= bp.dist)] <- clumpmatrix[which(d <=  bp.dist)] + 1
  clumpmatrix[which(snp.cor >= r2)] <- clumpmatrix[which(snp.cor >=  r2)] + 3
  
  marker.names <- rownames(assoc.signif)
  rownames(clumpmatrix) <- colnames(clumpmatrix) <- marker.names
  used <- rep(0, times = length(signif.p2))
  clumps <- list()
  for (i in 1:length(signif.p1)) {
    marker <- signif.p1[i]
    marker.index <- which(signif.p2 == marker)
    if (used[marker.index] == 0) {
      used[marker.index] <- 1
      clump <- which(clumpmatrix[marker, ] == 4)
      unused <- which(used[clump] == 0)
      clump <- clump[unused]
      if (length(clump) > 0) {
        p <- paste("Marker ", marker, " clumps with markers: ", 
                   paste(marker.names[clump], collapse = ", "), 
                   sep = "")
        snpnames <- c(marker, marker.names[clump])
        clumps[[marker]] <- assoc.signif[marker.names %in% snpnames, ]
        if (verbose) {
          print(p)
        }
      }
      used[clump] <- 1
    }
  }
  if (image == T) {
    par(mfrow = c(1, 3))
    image(snp.cor, col = rev(heat.colors(100)), main = "r2 matrix")
    image(d, col = rev(heat.colors(100)), main = "distance matrix")
    image(clumpmatrix, col = c("cornsilk1", "blue", "tomato", "red"), main = "clumping matrix")
  }
  if (length(clumps) == 0) {
    warning("No clumps found in dataset!")
  }
  clumps
}


