






####################
## manhattan.plot ##
####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Manhattan Plot
#'
#' Generate a Manhattan plot showing the association score values or p-values (y-axis)
#' for each locus (x-axis) tested by an association test.
#'
#' @param p.vals A numeric vector containing p-values or association score values for each genetic locus.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @import adegenet

########################################################################

## EG:
# manhattan.plot(p.vals = corr.dat,
#                col = "deepseasun", # c("red", "royalblue"),
#                transp = 0.5,
#                sig.thresh = c(0.68,0.71,0.75,0.80,0.81),
#                jitter.amount = 0.00001,
#                min.p = NULL,
#                log10=FALSE,
#                thresh.col="wasp",
#                ylab = "Terminal Correlation Score")

manhattan.plot <- function(p.vals,
                           x = c(1:length(p.vals)),
                           col = "funky",
                           transp = 0.25,
                           sig.thresh = NULL,
                           thresh.col="red",
                           snps.assoc = NULL,
                           snps.assoc.col = "red",
                           jitter.amount = 0.00001,
                           min.p = NULL,
                           log10 = FALSE,
                           ylab = NULL,
                           main.title="Manhattan plot"){


  # require(adegenet) # transp, col.pal

  pval <- p.vals
  if(is.null(x)) x <- c(1:length(p.vals))

  # Handle thresholds
  if(!is.null(sig.thresh)){
    if(class(sig.thresh) == "list") sig.thresh <- as.vector(unlist(sig.thresh))
    sig.thresh.complete <- sig.thresh
    sig.thresh <- unique(sig.thresh)
  }

  ## replace any "0" p.vals with min.p
  if(any(pval[!is.na(pval)] == 0)){
    pv <- pval[!is.na(pval)]
    toReplace <- which(!is.na(pval))[which(pv == 0)]
    if(is.null(min.p)) min.p <- 1/length(pval)
    pval[toReplace] <- min.p
  }

  ## PLOTTING METHOD--MANHATTAN PLOT

  # with Bonferroni:
  if(log10 == TRUE){
    log.pval <- -log10(pval)
  }else{
    log.pval <- pval
  }
  set.seed(1)
  if(jitter.amount > 0) log.pval <- jitter(log.pval, amount=jitter.amount)

  ## get colour scheme
  ## from colour palette?
  col.pals <- c("bluepal", "redpal", "greenpal", "greypal",
                "flame", "azur",
                "seasun", "lightseasun", "deepseasun",
                "spectral", "wasp", "funky")

  if(col %in% col.pals){
    ## get unit length (base 10) for coloured sections:
    l <- length(p.vals)
    n1 <- keepFirstN(l, 1)
    # n1 <- as.numeric(n1) + 1
    n1 <- as.numeric(n1)
    n0 <- nchar(l)-1
    n0 <- rep(0, n0)
    n0 <- paste0(n0, collapse="")
    N <- as.numeric(paste0(n1, n0, collapse=""))
    N <- round(N/10)
    if(N < 1) N <- 1

    myCol <- eval(parse(text=paste(col, "(", 10, ")")))
    myCol <- as.vector(unlist(sapply(c(1:length(myCol)),
                                     function(e)
                                       rep(myCol[e], N))))
    myCol <- myCol[c(1:length(pval))]
  }else{
    ## from vector of colours?
    myCol <- as.vector(unlist(sapply(c(1:length(col)),
                                     function(e)
                                       rep(col[e], N))))
    # myCol <- sort(myCol)
    myCol <- rep(myCol, ceiling(length(pval)/(N*(length(col)))))
    myCol <- myCol[1:length(pval)]
  }

  ## add transparency?
  if(!is.null(transp)){
    if(transp < 0 | transp > 1){
      transp <- 0.5
      warning("transp must be between 0 and 1;
              using default value 0.5.")
    }
    transp <- 1-transp
    myCol <- transp(myCol, alpha = transp)
    }

  ## check no myCol = ~white!
  cols <- col2rgb(unique(myCol))
  toChange <- which(sapply(c(1:ncol(cols)), function(e) length(which(cols[,e] >= 240))) == 3)
  if(length(toChange) > 0){
    cols.new <- cols[,toChange]
    if(length(toChange) == 1){
      cols.new <- rgb(red=(cols.new[1]-30)/255, green=(cols.new[2]-30)/255, blue=(cols.new[3]-30)/255)
    }else{
      cols.new <- sapply(c(1:ncol(cols.new)), function(e)
        rgb(red=(cols.new[1,e]-30)/255, green=(cols.new[2,e]-30)/255, blue=(cols.new[3,e]-30)/255))
    }
    for(e in 1:length(cols.new)){
      myCol[which(myCol %in% unique(myCol)[toChange[e]])] <- transp(cols.new[e], alpha = transp)
    }
  }

  ## Get y-max (set above max or sig.thresh):
  if(log10 == TRUE){
    ymax <- max(log.pval[!is.na(pval)])
    if(ymax > max(-log10(sig.thresh))){
      ymax <- ymax + (ymax/10)
    }else{
      ymax <- max(-log10(sig.thresh))
      ymax <- ymax + (ymax/10)
    }
  }else{
    ymax <- max(pval[!is.na(pval)])
    if(ymax > max(sig.thresh)){
      ymax <- ymax + (ymax/10)
    }else{
      ymax <- max(sig.thresh)
      ymax <- ymax + (ymax/10)
    }
  }


  ##p- Manhattan- Bonferroni
  if(log10 == TRUE){
    if(is.null(ylab)){
      ylab <- "Uncorrected -log10(p-value)"
    }
    plot(x = x,
         y = log.pval,
         col = myCol,
         pch = 19,
         cex = 1,
         main=main.title,
         xlab="genetic loci",
         ylab=ylab,
         cex.main=1,
         cex.axis=1.1,
         cex.lab=1.2,
         ylim=c(0,ymax))

  }else{
    if(is.null(ylab)){
      ylab <- "Uncorrected p-value"
    }
    plot(x = x,
         y = log.pval,
         col = myCol,
         pch = 19,
         cex = 1,
         main=main.title,
         xlab="genetic loci",
         ylab=ylab,
         cex.main=1,
         cex.axis=1.1,
         cex.lab=1.2,
         ylim=c(0,ymax))
  }




  ## overlay/highlight snps.assoc?
  if(!is.null(snps.assoc)){
    ## modify colour
    if(is.null(snps.assoc.col)) snps.assoc.col <- "red"
    myCol[snps.assoc] <- snps.assoc.col
    ## overlay points:
    points(x=snps.assoc,
           y=log.pval[snps.assoc],
           col = myCol[snps.assoc],
           pch = 19,
           cex = 1)
  }

  ## get significance threshold
  if(thresh.col %in% col.pals){
    thresh.col <- eval(parse(text=paste(thresh.col, "(", length(sig.thresh), ")")))
  }

  for(i in 1:length(sig.thresh)){
    if(log10 == TRUE){
      thresh <- -log10(sig.thresh[i])
    }else{
      thresh <- sig.thresh[i]
    }
    ## move sig thresh below nearest points by 0.05% of ylim
    thresh <- thresh-(max(log.pval[!is.na(log.pval)])*0.0005)
    ## draw threshold line on plot:
    # abline(h=thresh, col = thresh.col[i], lwd = 2)
    # lines(x=c(-400, length(log.pval+400)), y=c(thresh, thresh), col=thresh.col[i], lwd=2)
    lines(x=c(-400, max(x)+400), y=c(thresh, thresh), col=thresh.col[i], lwd=2)
  } # end for loop plotting thresh lines

} # end manhattan.plot



##########################
## qqman Manhattan Plot ##
##########################
# require(qqman)
#
# # manhattan(gwasResults,
# #           main = "Manhattan Plot",
# #           ylim = c(0, 10),
# #           cex = 0.6,
# #           cex.axis = 0.9,
# #           col = c("blue4", "orange3"),
# #           suggestiveline = F,
# #           genomewideline = F,
# #           chrlabs = c(1:20, "P", "Q"))
#
# BP <- as.numeric(dimnames(snps.ori)[[2]])
# if(is.null(BP)) BP <- c(1:ncol(snps.ori))
# CHR <- sort(rep(1:20, ncol(snps.ori)/20))
# P <- p.vals
# df <- data.frame(BP, CHR, P)
#
# manhattan(df,
#           main = "Manhattan Plot",
#           ylim = c(0, 10),
#           cex = 0.6,
#           cex.axis = 0.9,
#           col = c("blue4", "orange3"),
#           suggestiveline = F,
#           genomewideline = F,
#           chrlabs = c(1:20, "P", "Q"))









###################
## plot_sig_snps ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Plot null distribution and significant sites.
#'
#' Plot a histogram of the null distribution,
#' indicating the significance threshold and
#' the names and association scores of significant sites.
#'
#' @param arg Description.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)

########################################################################
# res <- out$res
# vals <- res$vals
#
# corr.dat <- vals$subsequent$corr.dat
# corr.sim <- vals$subsequent$corr.sim
#
# sig.thresh <- as.vector(unlist(res$thresh))
#
# snps.assoc <- out$performance[[1]]
# sig.corrs <- corr.dat[snps.assoc]
# sig.snps <- snps.assoc

plot_sig_snps <- function(corr.dat,
                          corr.sim,
                          corr.sim.subset = NULL,
                          sig.corrs = NULL,
                          sig.snps = NULL,
                          sig.thresh = NULL,
                          test = NULL,
                          sig.snps.col = "blue",
                          hist.col = rgb(0,0,1,0.5), # blue ## OR ## rgb(0.1,0.1,0.1,0.5) # darkgrey
                          hist.subset.col = rgb(1,0,0,0.5), # red ## OR ## rgb(0.8,0.8,0.8,0.5) # lightgrey
                          thresh.col="seasun",
                          snps.assoc = NULL,
                          snps.assoc.col = "red",
                          bg = "lightgray",
                          grid=TRUE,
                          freq = FALSE,
                          plot.null.dist=TRUE,
                          plot.dist=FALSE,
                          main.title=TRUE,
                          ...){

  ## Use Abs Vals (?)
  # corr.dat <- abs(corr.dat)
  # corr.sim <- abs(corr.sim)

  ###############################################################################
  ## Add threshold(s) and SNP annotations to histogram of (null) distributions ##
  ###############################################################################

  # thresh <- sig.thresh
  if(!is.null(sig.thresh)){
    if(class(sig.thresh) == "list") sig.thresh <- as.vector(unlist(sig.thresh))
    sig.thresh.complete <- sig.thresh
    sig.thresh <- unique(sig.thresh)
  }

  if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=TRUE)

  h.null <- h.null.subset <- ymax <- xmax <- y.lab <- NULL

  if(freq == TRUE){
    y.lab <- "Frequency"
  }
  if(freq == FALSE){
    y.lab <- "Density"
  }

  if(test == "fisher"){
    x.lab <- paste("-log10", test, "p-values", sep=" ")
  }else{
    x.lab <- paste(test, "score", sep=" ")
  }

  ##############################
  ## plot.null.dist ############
  ##############################

  ###########################################################
  ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
  ###########################################################
  ## plot correlations btw simulated SNPs and phenotype:

  if(plot.null.dist==TRUE){

    h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)

    ## add second (subset) histogram):
    if(!is.null(corr.sim.subset)){
      brks <- length(h.null$breaks)
      h.null.subset <- hist(as.vector(unlist(corr.sim.subset)), breaks = brks, plot=FALSE)
      ## if no subset col provided, choose 2:
      if(is.null(hist.subset.col)){
        # hist.col <- rgb(0,0,1,0.5) # blue
        # hist.subset.col <- rgb(1,0,0,0.5) # red
        # warning("No colour provided for subset:
        #         Choosing primary and subset colours.")
        hist.subset.col <- hist.col ## (makes more sense if expecting a SUBSET)
      }
    }

    ## Get plot limits:
    ## Get x-max:

    xmax <- max(corr.dat[!is.na(corr.dat)])+.3*max(corr.dat[!is.na(corr.dat)])
    if(test == "terminal" & xmax > 1){
      # xmax <- round(xmax, 1)
      xmax <- 1
    }else{
      xmax <- ceiling(xmax)
    }
    ## Get y-max (for overlaying hists, sig.thresh, sig.loci):
    if(freq == TRUE){
      yvals <- h.null$counts
      ymax <- ceiling(max(yvals)+(0.2*max(yvals)))
      # if(!is.null(h.null.subset)) ymax <- max(ymax, ceiling(max(h.null.subset$counts)+.05))
      if(!is.null(h.null.subset)) ymax <- max(ymax, ceiling(max(h.null.subset$counts)+(0.2*max(h.null.subset$counts))))
    }else{
      yvals <- h.null$density
      ymax <- max(yvals)+(0.2*max(yvals))
      # if(!is.null(h.null.subset)) ymax <- max(ymax, max(h.null.subset$density)+0.005)
      if(!is.null(h.null.subset)) ymax <- max(ymax, max(h.null.subset$density)+(0.2*max(h.null.subset$density)))
    }

    ## if the true correlation value for SNP i is >
    ## max bin, then extend the x-axis of the plot
    ## to accommodate annotation:
    if(xmax > max(h.null$breaks)){

      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## EXTENDING THE X-AXIS
      plot(h.null,
           main = NULL,
           xlab=x.lab,
           ylab = y.lab,
           xlim=c(min(h.null$breaks), xmax),
           ylim=c(0, ymax),
           freq=freq,
           col=hist.col,
           cex.axis=1.1,
           cex.lab=1.2)
      # ...)


      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)
          # par(bg=bg)

          # axis(1) ## add axes back
          # axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main = NULL,
               xlab = NULL,
               ylab = NULL,
               xaxt = "n",
               yaxt = "n",
               cex.axis = 0,
               xlim=c(min(h.null$breaks), xmax),
               ylim=c(0, ymax),
               col=hist.col,
               freq=freq,
               add = TRUE)

        }
      } # end background


      ## Overlay subset histogram:
      if(!is.null(h.null.subset)){
        plot(h.null.subset,
             main = NULL,
             xlim=c(min(h.null$breaks), xmax),
             ylim=c(0, ymax),
             xlab = NULL,
             ylab = NULL,
             xaxt = "n",
             yaxt = "n",
             cex.axis = 0,
             col=hist.subset.col,
             freq=freq,
             add=TRUE)
      }

    }else{

      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## WITHOUT EXTENDING THE X-AXIS
      plot(h.null,
           main = NULL,
           xlab=x.lab,
           ylab = y.lab,
           xaxt = "n",
           yaxt = "n",
           ylim=c(0, ymax),
           freq=freq,
           col=hist.col,
           ...)

      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)
          # par(bg=bg)

          # axis(1) ## add axes back
          # axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main = NULL,
               xlab = NULL,
               ylab = NULL,
               xaxt = "n",
               yaxt = "n",
               cex.axis = 0,
               ylim=c(0, ymax),
               col=hist.col,
               freq=freq,
               add = TRUE)

        }
      } # end background

      ## Overlay subset histogram:
      if(!is.null(h.null.subset)){
        plot(h.null.subset,
             main = NULL,
             xlab = NULL,
             ylab = NULL,
             xaxt = "n",
             yaxt = "n",
             cex.axis = 0,
             xlim=c(min(h.null$breaks), max(h.null$breaks)),
             ylim=c(0, ymax),
             col=hist.subset.col,
             freq=freq,
             add=TRUE)
      }

    }

    ## Add box around periphery:
    box()

    #################
    ## Add density curves?
    # d <- density(dat1) # from=min(hist(dat1, plot=FALSE)$breaks) ## from=0
    # ymax <- ceiling(max(d$y))
    # # hist(dat, freq=F, xlim=c(0,1), ylim=c(0,ymax))
    # # lines(d, col="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))
    # polygon(d, col=transp("red", 0.25), border="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))
    #################



    ## get significant loci (incl. snps.assoc)
    sig.loc <- NULL
    if(!is.null(sig.snps)) sig.loc <- sig.snps
    if(!is.null(snps.assoc)) sig.loc <- c(sig.loc, snps.assoc)
    sig.loc <- unique(sig.loc)

    # myCol <- rep("blue", length(sig.loc))

    ## sig.snps.col
    if(!is.null(sig.loc)){
      ## get colour
      if(is.null(sig.snps.col)) sig.snps.col <- "blue"
      myCol <- rep(sig.snps.col, length(sig.loc))
    }

    ## snps.assoc.col
    if(!is.null(snps.assoc)){
      ## modify colour
      if(is.null(snps.assoc.col)) snps.assoc.col <- "red"
      myCol[which(sig.loc %in% snps.assoc)] <- snps.assoc.col
    }

    ################################
    # plot(h.null,
    #      main = NULL,
    #      xlab=paste(test, "score", sep=" "),
    #      ylab = y.lab,
    #      xlim=c(min(h.null$breaks), xmax),
    #      ylim=c(0, ymax),
    #      freq=freq,
    #      col=hist.col,
    #      cex.axis=1.1,
    #      cex.lab=1.2)
    ##############################

    ## overlay/highlight sig.loc
    if(length(sig.loc) > 0){
      X <- corr.dat[sig.loc]

      ymin <- ymax/10 # 10%
      ymax <- ymax - ymax/10 # 90%
      Y <- seq(ymin, ymax, length.out = length(X))
      # Y <- c(seq(ymin, ymax, length.out = length(X)/2), seq(ymin, ymax, length.out = length(X)/2))

      ## ADD arrows pointing from each label to
      ## position on X-axis:
      arrows(x0=X , y0=Y ,
             x1=X , y1=0 , col=myCol,
             length=0.1, lwd=1)
      ## add annotation text labelling SNPs >
      ## threshold at their location on the x-axis:
      text(x=X, y=Y,
           # labels=sig.loc,  # sig.snps.names2,
           labels=sig.snps,
           col=myCol, font=1, pos=4, cex=1) # font=2 # cex=1

    }else{
      text(x=(max(h.null$breaks)*3/4),
           y=(max(yvals)*7/8),
           labels="no significant SNPs found",
           col="red", font=3, pos=3, cex = 1)
    }


    #######
    ## get significance threshold(s)
    ## get colour scheme
    ## from colour palette?
    col.pals <- c("bluepal", "redpal", "greenpal", "greypal",
                  "flame", "azur",
                  "seasun", "lightseasun", "deepseasun",
                  "spectral", "wasp", "funky")
    if(length(sig.thresh) > 0){

      if(thresh.col %in% col.pals){
        my.thresh.col <- eval(parse(text=paste(thresh.col, "(", length(sig.thresh), ")")))
      }else{
        my.thresh.col <- rep(thresh.col, length(sig.thresh))
      }
      thresh.col <- my.thresh.col

      for(i in 1:length(sig.thresh)){

        thresh <- sig.thresh[i]

        ## move sig thresh below nearest points?
        # thresh <- thresh-0.05

        ## draw threshold line on plot:
        # abline(v=thresh, col = my.thresh.col[i], lwd = 2)
        lines(x=c(thresh, thresh), y=c(0, max(yvals)), col=thresh.col[i], lwd=2)
      } # end for loop plotting thresh lines
    }
    #######





    ## Add plot title:
    if(!is.null(main.title)){
      if(main.title != FALSE){
        if(is.character(main.title)){
          title(main.title)
        }else{
          if(!is.null(test)){
            title(paste("Null distribution \n(", test, "score)"
                        # \n (with significant SNPs indicated)"
                        , sep=" "))
          }else{
            title(paste("Null distribution \n"
                        # \n (with significant SNPs indicated)"
                        , sep=" "))
          }
        }
      }
    }

  } # end plot.null.dist




  ################################
  ## plot.dist ###################
  ################################

  #################################################################
  ## plot histogram of correlations btw real SNPs and phenotype: ##
  #################################################################

  if(plot.dist == TRUE){

    h.null <- hist(as.vector(unlist(corr.dat)), plot=FALSE)

    xmax <- ceiling(max(corr.dat)+.35*max(corr.dat))
    ## Get y-max (for overlaying hists, sig.thresh, sig.loci):
    if(freq == TRUE){
      yvals <- h.null$counts
      ymax <- ceiling(max(yvals)+(0.2*max(yvals)))
      if(!is.null(h.null.subset)) ymax <- max(ymax, ceiling(max(h.null.subset$counts)+(0.2*max(h.null.subset$counts))))
    }else{
      yvals <- h.null$density
      ymax <- max(yvals)+(0.2*max(yvals))
      if(!is.null(h.null.subset)) ymax <- max(ymax, max(h.null.subset$density)+(0.2*max(h.null.subset$density)))
    }


    ## if the true correlation value for SNP i is >
    ## max bin, then extend the x-axis of the plot
    ## to accommodate annotation:
    if(xmax > max(h.null$breaks)){
      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## EXTENDING THE X-AXIS
      plot(h.null,
           main = NULL,
           xlab=paste(test, "score", sep=" "),
           ylab = y.lab,
           xlim=c(min(h.null$breaks), xmax),
           freq=freq,
           col=hist.col)
      # ,...)

      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)
          # par(bg=bg)

          # axis(1) ## add axes back
          # axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main = NULL,
               xlab = NULL,
               ylab = NULL,
               xlim=c(min(h.null$breaks), xmax),
               col=hist.col,
               freq=freq,
               add = TRUE)

        }
      } # end background

    }else{
      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## WITHOUT EXTENDING THE X-AXIS
      plot(h.null,
           main = NULL,
           xlab=paste(test, "score", sep=" "),
           ylab = y.lab,
           freq=freq,
           col=hist.col)
      # ,...)

      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)
          # par(bg=bg)

          # axis(1) ## add axes back
          # axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main = NULL,
               xlab = NULL,
               ylab = NULL,
               col=hist.col,
               freq=freq,
               add = TRUE)

        }
      } # end background
    }

    ## Add box around periphery:
    box()



    ## get significant loci (incl. snps.assoc)
    sig.loc <- NULL
    if(!is.null(sig.snps)) sig.loc <- sig.snps
    if(!is.null(snps.assoc)) sig.loc <- c(sig.loc, snps.assoc)
    sig.loc <- unique(sig.loc)

    # myCol <- rep("blue", length(sig.loc))

    ## sig.snps.col
    if(!is.null(sig.loc)){
      ## get colour
      if(is.null(sig.snps.col)) sig.snps.col <- "blue"
      myCol <- rep(sig.snps.col, length(sig.loc))
    }

    ## snps.assoc.col
    if(!is.null(snps.assoc)){
      ## modify colour
      if(is.null(snps.assoc.col)) snps.assoc.col <- "red"
      myCol[which(sig.loc %in% snps.assoc)] <- snps.assoc.col
    }


    ## overlay/highlight sig.loc
    if(length(sig.loc) > 0){
      X <- corr.dat[sig.loc]
      ymin <- max(yvals)/10 # 10%
      ymax <- max(yvals) - max(yvals)/10 # 90%

      Y <- seq(ymin, ymax, length.out = length(X))

      ## ADD arrows pointing from each label to
      ## position on X-axis:
      arrows(x0=X , y0=Y ,
             x1=X , y1=0 , col=myCol,
             length=0.1, lwd=1)
      ## add annotation text labelling SNPs >
      ## threshold at their location on the x-axis:
      text(x=X, y=Y, labels=sig.loc,
           col=myCol, font=1, pos=4, cex = 1)

    }else{
      text(x=(max(h.null$breaks)*3/4),
           y=(max(yvals)*3/4),
           labels="no significant SNPs found",
           col="red", font=3, pos=3, cex = 1)
    }



    #################

    ## get significance threshold(s)
    ## get colour scheme
    ## from colour palette?
    col.pals <- c("bluepal", "redpal", "greenpal", "greypal",
                  "flame", "azur",
                  "seasun", "lightseasun", "deepseasun",
                  "spectral", "wasp", "funky")
    if(thresh.col %in% col.pals){
      my.thresh.col <- eval(parse(text=paste(thresh.col, "(", length(sig.thresh), ")")))
    }else{
      my.thresh.col <- rep(thresh.col, length(sig.thresh))
    }

    for(i in 1:length(sig.thresh)){

      thresh <- sig.thresh[i]

      ## move sig thresh below nearest points
      # thresh <- thresh-0.05

      ## draw threshold line on plot:
      # abline(v=thresh, col = my.thresh.col[i], lwd = 2)
      lines(x=c(thresh, thresh), y=c(0, max(yvals)), col=thresh.col[i], lwd=2)
    } # end for loop plotting thresh lines

    #######



    ## Add plot title:
    if(!is.null(main.title)){
      if(main.title != FALSE){
        if(is.character(main.title)){
          title(main.title)
        }else{
          if(!is.null(test)){
            title(paste("Empirical distribution \n(", test, "score)"
                        # \n (with significant SNPs indicated)"
                        , sep=" "))
          }else{
            title(paste("Empirical distribution \n"
                        # \n (with significant SNPs indicated)"
                        , sep=" "))
          }
        }
      }
    }

  } # end plot.dist == TRUE


  if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=FALSE)

} # end plot_sig_snps




#################################
##  ENABLE ALTERNATE FN NAME:  ##
#################################
# plot.sig.snps <- function(corr.dat, corr.sim, ...){
#   return(plot_sig_snps(corr.dat, corr.sim,  ...))
# } # end plot.sig.snps




