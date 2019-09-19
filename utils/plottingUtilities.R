require(ggplot2)
require(scales)
require(gridExtra)
require(grid)

generateImage <- function(template, data) {
  # plot expression levels onto a template that defines the embryo ellipse
  # args:
  #  template: a numeric matrix with the same size as the desired image. An
  #   entry of 1 indicates embryo pixels and 0 indicates non-embryo
  #  data: a column vector of values to plot, with length equal to the number of
  #   1s in template
  embryo.idcs <- which(template == 1, arr.ind=TRUE)
  outer.idcs <- which(template == 0, arr.ind=TRUE)
  n.row <- nrow(template)
  n.col <- ncol(template)

  full.img <- matrix(0, nrow=n.row, ncol=n.col)
  for (j in 1:nrow(embryo.idcs)) {
    full.img[embryo.idcs[j,1], embryo.idcs[j,2]] <- data[j]
  }

  if (length(outer.idcs)) {
    for (j in 1:nrow(outer.idcs)) {
      full.img[outer.idcs[j,1], outer.idcs[j,2]] <- NA
    }
  }

  return(full.img)
}


plotImg <- function(x, template, legend, label) {
  # Main plotting function for single embryo image. Plots expression level data
  # in x onto embryo template.
  #
  # args:
  #  x: vector of gene expression levels by pixel. Must be equal in length to
  #   sum(template).
  #  template: binary matrix indicting embryo pixels
  #  legend: TRUE/FALSE should a legend showing expression levels be plotted?

  img <- t(generateImage(template, x)[nrow(template):1,])
  img.grid <- expand.grid(1:nrow(img), 1:ncol(img))
  img.grid <- as.data.frame(cbind(img.grid, c(img)))
  names(img.grid ) <- c('x', 'y', 'expr')
  p <- ggplot(img.grid, aes(x, y)) + geom_raster(aes(fill=expr))
  p <- p + scale_fill_gradient2()

  p <- p + ggtitle(label) + theme(axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 8))
  if (!legend) p <- p + theme(legend.position="none")
  return(p)
}

plotGroupImg <- function(x, template, vectorized=TRUE, legend=FALSE) {
  # Wrapper of plotImg to plot multiple images. These can be given as expression
  # level vecotrs (default) or already mapped on to the template.
  if(!is.matrix(x)) x <- as.matrix(x)
  if (!vectorized) x <- apply(x, MAR=3, function(z) z[template == 1])
  x.names <- colnames(x)
  x <- split(x, rep(1:ncol(x), each = nrow(x)))
  plot.list <- lapply(x, plotImg, template=template, legend=legend)

  for (i in 1:length(plot.list)) {
    plot.list[[i]] <- plot.list[[i]] + ggtitle(x.names[i])
  }
  do.call('grid.arrange', c(plot.list))
}
getNodePosition <- function(g) return(lapply(g@AgNode, function(n) c(n@center@x,
                                                                     n@center@y)))
getNodeLabel <- function(g) return(lapply(g@AgNode, function(n) n@name))

overlayImage <- function(image, center, label, color.pal, buffer=10) {
  im.width <- ncol(image)
  im.height <- nrow(image)
  img.col <- img2Color(image, color.pal)
  rasterImage(img.col, center[1] - im.width / 2, center[2] - im.height / 2,
              center[1] + im.width / 2, center[2] + im.height / 2)
  text(center[1], center[2] + im.height / 2 + buffer, label)
}

overlayImages <- function(graph, images, color.pal) {
  centers <- getNodePosition(graph)
  labels <- getNodeLabel(graph)
  mapply(function(cent, lab, img) overlayImage(img, cent, lab, color.pal,
                                               buffer=8),
         centers, labels, images)
}

getNodeImages <- function(graph, x.list, color.pal) {
  nodes <- getNodeLabel(graph)
  imgs <- lapply(nodes, function(n)
                 generateImage(x.list$template, as.matrix(x.list$x[,n])))
  return(imgs)
}

img2Color <- function(img, color.pal) {
  img.dim <- dim(img)
  n.cols <- length(color.pal)
  img.cut <- cut(img, n.cols)
  img.col <- matrix(color.pal[img.cut], nrow=img.dim[1])
  img.col[is.na(img.col)] <- '#A9A9A9'
  return(img.col)
}



plotPP <- function(x, dict, pp, template) {
  pp.pixels <- sapply(pp, function(p) which(dict[,p] > 0.1))
  pp.pixels <- unique(unlist(pp.pixels))
  x[-pp.pixels,] <- 0
  plotGroupImg(x, template)
}

shiftImg <- function(x, template, dir, by, vectorized=TRUE) {
  if (!vectorized) x <- x[template == 1]
  img <- t(generateImage(template, x)[nrow(template):1,])
  img[is.na(img)] <- 0
  if (dir == 'horiz') {
    img <- shiftHoriz(img, by)
  } else if (dir == 'vert') {
    img <- shiftVert(img, by)
  }
  img <- t(img)[nrow(template):1,]
  x <- img[template == 1]
  return(x)
}

shiftVert <- function(x, by) {
  magnitude <- abs(by)
  nrow.x <- nrow(x)
  ncol.x <- ncol(x)
  if (by > 0) {
    x.sub <- x[,(magnitude + 1):ncol.x]
    x <- cbind(x.sub, matrix(0, nrow=nrow.x, ncol=magnitude))
  } else if (by < 0) {
    x.sub <- x[,1:(ncol.x - magnitude)]
    x <- cbind(matrix(0, nrow=nrow.x, ncol=magnitude), x.sub)
  }
  return(x)
}

shiftHoriz <- function(x, by) {
  magnitude <- abs(by)
  nrow.x <- nrow(x)
  ncol.x <- ncol(x)
  if (by < 0) {
    x.sub <- x[(magnitude + 1):nrow.x,]
    x <- rbind(x.sub, matrix(0, nrow=magnitude, ncol=ncol.x))
  } else if (by > 0) {
    x.sub <- x[1:(nrow.x - magnitude),]
    x <- rbind(matrix(0, nrow=magnitude, ncol=ncol.x), x.sub)
  }
  return(x)
}

mapLate2Early <- function(late.atom, pp, early.dict, late.template, early.template, cors) {

  p.late <- plotImg(late.atom ^ 2, late.template, FALSE)
  p.late <- p.late + ggtitle(paste0('Stage 9-10 PP', pp))
  early.scored <- early.dict %*% cors

  if (sum(early.scored) > 0) {
    p.early <- plotImg(early.scored ^ 2, early.template, FALSE)
    p.early <- p.early + ggtitle('Associated stage 4-6 PP regions')
  } else {
    out.data <- data.frame(x=0, y=0, label='No strongly associated stage 4-6 PP regions')
    p.early <- ggplot(out.data, aes(x=x, y=y, label=label))
    p.early <- p.early + theme(axis.line=element_blank(),
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   panel.background=element_blank(),
                   panel.border=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   plot.background=element_blank())
    p.early <- p.early + geom_text()
  }
  grid.arrange(p.late, p.early)

}

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, widths=NULL, heights=NULL,
                      title=NULL, titlefont = "", titleface = 1, titlesize = 16) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (!is.null(title)) { # Add a narrow row at the top for the title
    layout <- rbind(rep(0,ncol(layout)),layout)
    if (is.null(heights)) {
      plotrows <- nrow(layout)-1
      rowheights <- c(0.1, rep(1,plotrows)/plotrows)
    } else {
      rowheights <- c(0.1, heights/sum(heights))
    }
  } else {
    if (is.null(heights)) {
      rowheights <- rep(1,nrow(layout))
    } else {
      rowheights <- heights
    }
  }

  if (is.null(widths)) {
    colwidths <- rep(1, cols)
  } else {
    colwidths <- widths
  }

  if (numPlots==1) {

    return(plots[[1]] + labs(title=title))

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),
                                               widths=colwidths,
                                               heights=rowheights)))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }

    if (!is.null(title)) {
      grid.text(title, vp = viewport(layout.pos.row = 1
                                     , layout.pos.col = 1:ncol(layout)),
                gp = gpar(fontfamily = titlefont, fontface = titleface,
                          fontsize = titlesize))
    }

  }
return(invisible(NULL))
}

plot_all_PPs <- function(dict, template, reference = NULL, ncols = 5){

  if(is.null(reference)){
    perm <- seq(1, dim(dict)[2], by= 1)
  }else{
    # we sort the dictionary by comparing it with a reference dictionary
    remove_col <- c()
    perm <- c()
    for(i in 1:dim(reference)[2]){
        diff <- t(as.matrix(dict)) %*% reference[, i]
        perm_ <- order(diff, decreasing = TRUE)

        perm_ <- perm_[!(perm_ %in% remove_col)]
        perm <- c(perm, perm_[1])
        remove_col <- c(remove_col, perm_[1])
    }

  }

    p <- list()
      for(i in 1:dim(dict)[2]){
        p[[i]] <- plotImg(dict[, perm[i]], late$template, FALSE, paste0('PP', i))
      }

      plot <- multiplot(plotlist = p, cols = ncols)

    return(perm)
}
