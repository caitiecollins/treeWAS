
####################
## manhattan.plot ##
####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param arg Description.
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
                           col = "wasp",
                           transp = 0.75,
                           sig.thresh = 0.05,
                           thresh.col="seasun",
                           snps.assoc = NULL,
                           snps.assoc.col = "red",
                           jitter.amount = 0.00001,
                           min.p = NULL,
                           log10=TRUE,
                           ylab=NULL){

  # require(adegenet) # transp, col.pal


  pval <- p.vals

  # Handle thresholds
  if(!is.null(sig.thresh)){
    if(class(sig.thresh) == "list") sig.thresh <- as.vector(unlist(sig.thresh))
    sig.thresh.complete <- sig.thresh
    sig.thresh <- unique(sig.thresh)
  }

  ## replace any "0" p.vals with min.p
  if(any(pval == 0)){
    toReplace <- which(pval == 0)
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
    myCol <- eval(parse(text=paste(col, "(", ceiling(length(pval)/1000), ")")))
    myCol <- as.vector(unlist(sapply(c(1:length(myCol)),
                                     function(e)
                                       rep(myCol[e], 1000))))
    myCol <- myCol[c(1:length(pval))]
  }else{
    ## from vector of colours?
    myCol <- as.vector(unlist(sapply(c(1:length(col)),
                                     function(e)
                                       rep(col[e], 1000))))
    myCol <- sort(myCol)
    myCol <- rep(myCol, ceiling(length(pval)/(1000*(length(col)))))
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



  ##p- Manhattan- Bonferroni
  if(log10 == TRUE){
    if(is.null(ylab)){
      ylab <- "Uncorrected -log10(p-value)"
    }
    plot(log.pval,
         col = myCol,
         pch = 19,
         cex = 1,
         main="Manhattan plot",
         xlab="SNP loci",
         ylab=ylab,
         cex.main=1)

  }else{
    if(is.null(ylab)){
      ylab <- "Uncorrected p-value"
    }
    plot(log.pval,
         col = myCol,
         pch = 19,
         cex = 1,
         main="Manhattan plot",
         xlab="SNP loci",
         ylab=ylab,
         cex.main=1)
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
    ## move sig thresh below nearest points
    thresh <- thresh-0.05
    ## draw threshold line on plot:
    abline(h=thresh, col = thresh.col[i], lwd = 2)
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
## plot.sig.snps ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param arg Description.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @import adegenet ape phangorn

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

plot.sig.snps <- function(corr.dat,
                          corr.sim,
                          corr.sim.subset = NULL,
                          sig.corrs,
                          sig.snps,
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
                          plot.null.dist=TRUE,
                          plot.dist=FALSE){

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
      h.null.subset <- hist(as.vector(unlist(corr.sim.subset)), plot=FALSE)
      ## if no subset col provided, choose 2:
      if(is.null(hist.subset.col)){
        # hist.col <- rgb(0,0,1,0.5) # blue
        # hist.subset.col <- rgb(1,0,0,0.5) # red
        # warning("No colour provided for subset:
        #         Choosing primary and subset colours.")
        hist.subset.col <- hist.col ## (makes more sense if expecting a SUBSET)
      }
    }

    xmax <- ceiling(max(corr.dat)+.05)

    ## if the true correlation value for SNP i is >
    ## max bin, then extend the x-axis of the plot
    ## to accommodate annotation:
    if(xmax > max(h.null$breaks)){

      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## EXTENDING THE X-AXIS
      plot(h.null,
           main=paste("Null distribution of", test, "scores"
           # \n (with significant SNPs indicated)"
                      , sep=" "),
           xlab=paste(test, "score", sep=" "),
           xlim=c(min(h.null$breaks), xmax),
           col=hist.col)


      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)
          # rect(lim[1],  lim[3], lim[2], lim[4], col=bg, border=NA)

          axis(1) ## add axes back
          axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main=paste("Null distribution of", test, "scores"
               # \n (with significant SNPs indicated)"
                          , sep=" "),
               xlab=paste(test, "score", sep=" "),
               xlim=c(min(h.null$breaks), xmax),
               col=hist.col,
               add = TRUE)

        }
      } # end background


      ## Overlay subset histogram:
      if(!is.null(h.null.subset)){
        plot(h.null.subset,
             xlim=c(min(h.null$breaks), xmax),
             col=hist.subset.col,
             add=TRUE)
      }

    }else{
      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## WITHOUT EXTENDING THE X-AXIS
      plot(h.null,
           main=paste("Null distribution of", test, "scores"
           # \n (with significant SNPs indicated)"
                      , sep=" "),
           xlab=paste(test, "score", sep=" "),
           col=hist.col)

      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)

          axis(1) ## add axes back
          axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main=paste("Null distribution of", test, "scores"
                # \n (with significant SNPs indicated)"
                          , sep=" "),
               xlab=paste(test, "score", sep=" "),
               col=hist.col,
               add = TRUE)

        }
      } # end background

      ## Overlay subset histogram:
      if(!is.null(h.null.subset)){
        plot(h.null.subset,
             xlim=c(min(h.null$breaks), max(h.null$breaks)),
             col=hist.subset.col,
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
      thresh <- thresh-0.05
      ## draw threshold line on plot:
      abline(v=thresh, col = my.thresh.col[i], lwd = 2)
    } # end for loop plotting thresh lines

    #######

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
      ymin <- max(h.null$counts)/10 # 10%
      ymax <- max(h.null$counts) - max(h.null$counts)/10 # 90%

      Y <- seq(ymin, ymax, length.out = length(X))

      ## ADD arrows pointing from each label to
      ## position on X-axis:
      arrows(x0=X , y0=Y ,
             x1=X , y1=0 , col=myCol,
             length=0.1, lwd=1)
      ## add annotation text labelling SNPs >
      ## threshold at their location on the x-axis:
      text(x=X, y=Y, labels=sig.loc,
           col=myCol, font=2, pos=4, cex=0.7)

    }else{
      text(x=(max(h.null$breaks)*3/4),
           y=(max(h.null$counts)*3/4),
           labels="no significant SNPs found",
           col=myCol, font=2, pos=2, cex = 0.7)
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

    xmax <- ceiling(max(corr.dat)+.05)

    ## if the true correlation value for SNP i is >
    ## max bin, then extend the x-axis of the plot
    ## to accommodate annotation:
    if(xmax > max(h.null$breaks)){
      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## EXTENDING THE X-AXIS
      plot(h.null,
           main=paste("Real distribution of", test, "scores"
                      # \n (with significant SNPs indicated)"
                      , sep=" "),
           xlab=paste(test, "score", sep=" "),
           xlim=c(min(h.null$breaks), xmax),
           col=hist.col)

      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)

          axis(1) ## add axes back
          axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main=paste("Real distribution of", test, "scores"
                      # \n (with significant SNPs indicated)"
                          , sep=" "),
               xlab=paste(test, "score", sep=" "),
               xlim=c(min(h.null$breaks), xmax),
               col=hist.col,
               add = TRUE)

        }
      } # end background

    }else{
      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      ## WITHOUT EXTENDING THE X-AXIS
      plot(h.null,
           main=paste("Real distribution of", test, "scores"
                      # \n (with significant SNPs indicated)"
                      , sep=" "),
           xlab=paste(test, "score", sep=" "),
           col=hist.col)

      ## Add grey background? ##
      if(!is.null(bg)){
        if(bg != "white"){
          lim <- par("usr")
          rect(lim[1],  lim[3], lim[2], lim[4], col=bg)

          axis(1) ## add axes back
          axis(2)

          ## add grid:
          if(grid == TRUE){
            grid(col="white", lwd=1, lty=1)
          }

          box()   ## and the plot frame

          ## Re-plot original plot:
          plot(h.null,
               main=paste("Real distribution of", test, "scores"
                      # \n (with significant SNPs indicated)"
                          , sep=" "),
               xlab=paste(test, "score", sep=" "),
               col=hist.col,
               add = TRUE)

        }
      } # end background
    }

    ## Add box around periphery:
    box()

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
      thresh <- thresh-0.05
      ## draw threshold line on plot:
      abline(v=thresh, col = my.thresh.col[i], lwd = 2)
    } # end for loop plotting thresh lines

    #######

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
      ymin <- max(h.null$counts)/10 # 10%
      ymax <- max(h.null$counts) - max(h.null$counts)/10 # 90%

      Y <- seq(ymin, ymax, length.out = length(X))

      ## ADD arrows pointing from each label to
      ## position on X-axis:
      arrows(x0=X , y0=Y ,
             x1=X , y1=0 , col=myCol,
             length=0.1, lwd=1)
      ## add annotation text labelling SNPs >
      ## threshold at their location on the x-axis:
      text(x=X, y=Y, labels=sig.loc,
           col=myCol, font=2, pos=4, cex = 0.7)

    }else{
      text(x=(max(h.null$breaks)*3/4),
           y=(max(h.null$counts)*3/4),
           labels="no significant SNPs found",
           col=myCol, font=2, pos=2, cex = 0.7)
    }

  } # end plot.dist == TRUE


  if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=FALSE)

} # end plot.sig.snps








#' ###################
#' ## plot.sig.snps ##
#' ###################
#'
#' ########################################################################
#'
#' ###################
#' ## DOCUMENTATION ##
#' ###################
#'
#' #' Short one-phrase description.
#' #'
#' #' Longer proper discription of function...
#' #'
#' #' @param arg Description.
#' #'
#' #'
#' #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' #' @export
#' #'
#' #' @import adegenet ape phangorn
#'
#' ########################################################################
#'
#'
#' plot.sig.snps <- function(corr.dat, corr.sim,
#'                           sig.corrs, sig.snps,
#'                           sig.thresh, test,
#'                           plot.null.dist=TRUE, plot.dist=FALSE){
#'
#'   thresh <- sig.thresh
#'
#'   if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=TRUE)
#'
#'   ## plot.null.dist #####################################
#'
#'   ###################################################
#'   ## Add SNP annotations to histogram of null dist ##
#'   ###################################################
#'
#'
#'   ##############
#'   ## TERMINAL ##
#'   ##############
#'   if(test=="terminal"){
#'     ## plot.null.dist ##
#'
#'     ###########################################################
#'     ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
#'     ###########################################################
#'     ## plot correlations btw simulated SNPs and phenotype:
#'     h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
#'
#'     ## get alternate (null dist) label heights:
#'     X <- thresh
#'     if(length(sig.corrs) > 0) X <- sig.corrs
#'     Y.null <- runif(n=length(sig.snps), min=0.00001,
#'                     max=max(h.null$counts))
#'
#'
#'     if(plot.null.dist==TRUE){
#'       ## if the true correlation value for SNP i is >
#'       ## max bin, then extend the x-axis of the plot
#'       ## to accommodate annotation:
#'       if(max(X) > max(h.null$breaks)){
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Terminal Correlation Scores
#'              \n (with significant SNPs indicated)",
#'              xlab="Terminal Correlation Score",
#'              xlim=c(min(h.null$breaks), max(X)+.05))
#'       }else{
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## WITHOUT EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Terminal Correlation Scores
#'              \n (with significant SNPs indicated)",
#'              xlab="Terminal Correlation Score"
#'         )
#'       }
#'
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h.null$counts)),
#'            labels="significance threshold", pos=2,
#'            col="grey", font=4)
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to
#'         ## position on X-axis:
#'         arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
#'                x1=X , y1=0 , col="blue",
#'                length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y.null, labels=sig.snps,
#'              col="red", font=2, pos=2)
#'       }else{
#'         text(x=thresh, y=(max(h.null$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=2)
#'       }
#'     } # end plot.null.dist
#'   } # end test Terminal
#'
#'
#'   ##################
#'   ## SIMULTANEOUS ##
#'   ##################
#'   if(test=="simultaneous"){
#'     ## plot.null.dist ##
#'
#'     ###########################################################
#'     ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
#'     ###########################################################
#'     ## plot correlations btw simulated SNPs and phenotype:
#'     h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
#'
#'     ## get alternate (null dist) label heights:
#'     X <- thresh
#'     if(length(sig.corrs) > 0) X <- sig.corrs
#'     Y.null <- runif(n=length(sig.snps), min=0.00001,
#'                     max=max(h.null$counts))
#'
#'
#'     if(plot.null.dist==TRUE){
#'       ## if the true correlation value for SNP i is >
#'       ## max bin, then extend the x-axis of the plot
#'       ## to accommodate annotation:
#'       if(max(X) > max(h.null$breaks)){
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Simultaneous Scores
#'              \n (with significant SNPs indicated)",
#'              xlab="Simultaneous Score",
#'              xlim=c(min(h.null$breaks), max(X)+.05))
#'       }else{
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## WITHOUT EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Simultaneous Scores
#'              \n (with significant SNPs indicated)",
#'              xlab="Simultaneous Score"
#'         )
#'       }
#'
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h.null$counts)),
#'            labels="significance threshold", pos=2,
#'            col="grey", font=4)
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to
#'         ## position on X-axis:
#'         arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
#'                x1=X , y1=0 , col="blue",
#'                length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y.null, labels=sig.snps,
#'              col="red", font=2, pos=2)
#'       }else{
#'         text(x=thresh, y=(max(h.null$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=2)
#'       }
#'     } # end plot.null.dist
#'   } # end test Simultaneous
#'
#'   ################
#'   ## SUBSEQUENT ##
#'   ################
#'   if(test=="subsequent"){
#'     ## plot.null.dist ##
#'
#'     ###########################################################
#'     ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
#'     ###########################################################
#'     ## plot correlations btw simulated SNPs and phenotype:
#'     h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
#'
#'     ## get alternate (null dist) label heights:
#'     X <- thresh
#'     if(length(sig.corrs) > 0) X <- sig.corrs
#'     Y.null <- runif(n=length(sig.snps), min=0.00001,
#'                     max=max(h.null$counts))
#'
#'
#'     if(plot.null.dist==TRUE){
#'       ## if the true correlation value for SNP i is >
#'       ## max bin, then extend the x-axis of the plot
#'       ## to accommodate annotation:
#'       if(max(X) > max(h.null$breaks)){
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Subsequent Scores
#'              \n (with significant SNPs indicated)",
#'              xlab="Subsequent Score",
#'              xlim=c(min(h.null$breaks), max(X)+.05))
#'       }else{
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## WITHOUT EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Subsequent Scores
#'              \n (with significant SNPs indicated)",
#'              xlab="Subsequent Score"
#'         )
#'       }
#'
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h.null$counts)),
#'            labels="significance threshold", pos=2,
#'            col="grey", font=4)
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to
#'         ## position on X-axis:
#'         arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
#'                x1=X , y1=0 , col="blue",
#'                length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y.null, labels=sig.snps,
#'              col="red", font=2, pos=2)
#'       }else{
#'         text(x=thresh, y=(max(h.null$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=2)
#'       }
#'     } # end plot.null.dist
#'   } # end test Subsequent
#'
#'
#'   #################
#'   ## CORRELATION ##
#'   #################
#'   if(test=="cor"){
#'     ## plot.null.dist ##
#'
#'     ###########################################################
#'     ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
#'     ###########################################################
#'     ## plot correlations btw simulated SNPs and phenotype:
#'     h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
#'
#'     ## get alternate (null dist) label heights:
#'     X <- thresh
#'     if(length(sig.corrs) > 0) X <- sig.corrs
#'     Y.null <- runif(n=length(sig.snps), min=0.00001,
#'                     max=max(h.null$counts))
#'
#'
#'     if(plot.null.dist==TRUE){
#'
#'       ## if the true correlation value for SNP i is >
#'       ## max bin, then extend the x-axis of the plot
#'       ## to accommodate annotation:
#'       if(max(X) > max(h.null$breaks)){
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of correlations
#'              \n (with significant SNPs indicated)",
#'              xlab="Correlation",
#'              xlim=c(min(h.null$breaks), max(X)+.05))
#'       }else{
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## WITHOUT EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of correlations
#'              \n (with significant SNPs indicated)",
#'              xlab="Correlation"
#'         )
#'       }
#'
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h.null$counts)),
#'            labels="significance threshold", pos=2,
#'            col="grey", font=4)
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to
#'         ## position on X-axis:
#'         arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
#'                x1=X , y1=0 , col="blue", length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y.null, labels=sig.snps,
#'              col="red", font=2)
#'       }else{
#'         text(x=thresh, y=(max(h.null$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=2)
#'       }
#'     } # end plot.null.dist
#'     } # end test cor
#'
#'
#'   #########################
#'   ## FISHER'S EXACT TEST ##
#'   #########################
#'   if(test=="fisher"){
#'     ## plot.null.dist ##
#'
#'     ###########################################################
#'     ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
#'     ###########################################################
#'     ## plot correlations btw simulated SNPs and phenotype:
#'     h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
#'
#'     ## get alternate (null dist) label heights:
#'     X <- thresh
#'     if(length(sig.corrs) > 0) X <- sig.corrs
#'     Y.null <- runif(n=length(sig.snps), min=0.00001,
#'                     max=max(h.null$counts))
#'
#'     if(plot.null.dist==TRUE){
#'       ## if the true correlation value for SNP i is <
#'       ## min bin, then extend the x-axis of the plot to
#'       ## accommodate annotation:
#'       if(min(X) < min(h.null$breaks)){
#'         ## plot histogram of correlations btw real
#'         ## SNPs and phenotype: ##
#'         ## EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Fisher's exact test p-values
#'              \n (with significant SNPs indicated)",
#'              xlab="p-value",
#'              xlim=c((min(X)-.05), max(h.null$breaks)))
#'       }else{
#'         ## plot histogram of correlations btw real SNPs and phenotype: ##
#'         ## WITHOUT EXTENDING THE X-AXIS
#'         plot(h.null,
#'              main="Null distribution of Fisher's exact test p-values
#'              \n (with significant SNPs indicated)",
#'              xlab="p-value"
#'         )
#'       }
#'
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h.null$counts)),
#'            labels="significance threshold", pos=4,
#'            col="grey", font=4)
#'       ## only ask to draw arrows if sig snps exist
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to
#'         ## position on X-axis:
#'         arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
#'                x1=X , y1=0 , col="blue", length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y.null, labels=sig.snps,
#'              col="red", font=2, pos=4)
#'       }else{
#'         text(x=thresh, y=(max(h.null$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=4)
#'       }
#'     } # end plot.null.dist
#'   } # end test fisher
#'
#'
#'   ## plot.dist ###################
#'
#'   #################################################################
#'   ## plot histogram of correlations btw real SNPs and phenotype: ##
#'   #################################################################
#'
#'   ###########
#'   ## SCORE ##
#'   ###########
#'   if(test=="score"){
#'     ## plot.dist ##
#'
#'     #################################################################
#'     ## plot histogram of correlations btw real SNPs and phenotype: ##
#'     #################################################################
#'     ## get histogram of correlations btw real SNPs and phenotype:
#'     h <- hist(corr.dat, plot=FALSE)
#'     ## get X and Y coords for labelling positions of all sig SNPs
#'     X <- sig.corrs
#'     ## get the number of sig SNPs in each bin of the histogram
#'     ## counting correlations==upper limit of each bin
#'     ## as falling within that bin...
#'     sig.counts <- sapply(c(1:(length(h$breaks)-1)),
#'                          function(e) length(X[which(X[which(X <=
#'                                                               h$breaks[(e+1)])] > h$breaks[e])]))
#'     ## keep only counts > 0
#'     sig.counts <- sig.counts[which(sig.counts > 0)]
#'     ## get label heights:
#'     Y <- list()
#'     ## get average height of labels
#'     Y.avg <- max(h$counts)/4
#'     ## for bins with > 1 sig SNP, adjust height
#'     if(length(sig.counts)!=0){
#'       for(i in 1:length(sig.counts)){
#'         if(!.is.integer0(sig.counts[i])){
#'           if(sig.counts[i]==1){
#'             Y[[i]] <- Y.avg
#'           }else{
#'             ## divide up the space on the y-axis into
#'             ## increments (adding 1 to the denomenator
#'             ## s.t y-max not exceeded)
#'             increment <- max(h$counts)/(sig.counts[i]+1)
#'             Y[[i]] <- increment*c(1:sig.counts[i])
#'           }
#'         }
#'       } # end for loop
#'     }
#'     Y <- as.vector(unlist(Y))
#'
#'     if(plot.dist==TRUE){
#'
#'       ## plot histogram of correlations btw real SNPs and phenotype: ##
#'       plot(h, main="Distribution of SNP-trait correlations
#'            \n (with significant SNPs indicated if present)")
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h$counts)),
#'            labels="significance threshold",
#'            col="grey", pos=2, font=4)
#'
#'       ## Only ask to draw arrows if at least 1 significant SNP:
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to position on X-axis:
#'         arrows(x0=X , y0=(Y-(max(h$counts)/50)) ,
#'                x1=X , y1=0 ,
#'                col="blue", length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y, labels=sig.snps,
#'              col="red", font=2, pos=2)
#'       }else{
#'         text(x=thresh, y=(max(h$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=2)
#'       }
#'     } # end plot.dist
#'   } # end test score
#'
#'   #################
#'   ## CORRELATION ##
#'   #################
#'   if(test=="cor"){
#'     ## plot.dist ##
#'
#'     #################################################################
#'     ## plot histogram of correlations btw real SNPs and phenotype: ##
#'     #################################################################
#'     ## get histogram of correlations btw real SNPs and phenotype:
#'     h <- hist(corr.dat, plot=FALSE)
#'     ## get X and Y coords for labelling positions of all sig SNPs
#'     X <- sig.corrs
#'     ## get the number of sig SNPs in each bin of the histogram
#'     ## counting correlations==upper limit of each bin as falling
#'     ## within that bin...
#'     sig.counts <- sapply(c(1:(length(h$breaks)-1)),
#'                          function(e) length(X[which(X[which(X <=
#'                                                               h$breaks[(e+1)])] > h$breaks[e])]))
#'     ## keep only counts > 0
#'     sig.counts <- sig.counts[which(sig.counts > 0)]
#'     ## get label heights:
#'     Y <- list()
#'     ## get average height of labels
#'     Y.avg <- max(h$counts)/4
#'     ## for bins with > 1 sig SNP, adjust height
#'     if(length(sig.counts)!=0){
#'       for(i in 1:length(sig.counts)){
#'         if(!.is.integer0(sig.counts[i])){
#'           if(sig.counts[i]==1){
#'             Y[[i]] <- Y.avg
#'           }else{
#'             ## divide up the space on the y-axis into
#'             ## increments (adding 1 to the denomenator
#'             ## s.t y-max not exceeded)
#'             increment <- max(h$counts)/(sig.counts[i]+1)
#'             Y[[i]] <- increment*c(1:sig.counts[i])
#'           }
#'         }
#'       } # end for loop
#'     }
#'     Y <- as.vector(unlist(Y))
#'
#'     if(plot.dist==TRUE){
#'
#'       ## plot histogram of correlations btw real
#'       ## SNPs and phenotype: ##
#'       plot(h, main="Distribution of SNP-trait
#'            correlations \n (with significant SNPs indicated if present)")
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h$counts)),
#'            labels="significance threshold",
#'            col="grey", pos=2, font=4)
#'
#'       ## Only ask to draw arrows if at least 1 significant SNP:
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each
#'         ## label to position on X-axis:
#'         arrows(x0=X , y0=(Y-(max(h$counts)/50)) ,
#'                x1=X , y1=0 , col="blue", length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs >
#'         ## threshold at their location on the x-axis:
#'         text(x=X, y=Y, labels=sig.snps, col="red", font=2)
#'       }else{
#'         text(x=thresh, y=(max(h$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=2)
#'       }
#'     } # end plot.dist
#'   } # end test cor
#'
#'
#'   #########################
#'   ## FISHER'S EXACT TEST ##
#'   #########################
#'   if(test=="fisher"){
#'     ## plot.dist ##
#'
#'     ############################################################
#'     ## plot histogram of correlations btw real SNPs and phen: ##
#'     ############################################################
#'
#'     ## get histogram of correlations btw real SNPs and phenotype:
#'     h <- hist(corr.dat, plot=FALSE)
#'     ## get X and Y coords for labelling positions of all sig SNPs
#'     X <- thresh
#'     if(length(sig.corrs) > 0) X <- sig.corrs
#'     ## get the number of sig SNPs in each bin of the histogram
#'     ## counting correlations==upper limit of each bin as
#'     ## falling within that bin...
#'     sig.counts <- sapply(c(1:(length(h$breaks)-1)),
#'                          function(e) length(X[which(X[which(X <=
#'                                                               h$breaks[(e+1)])] > h$breaks[e])]))
#'     ## keep only counts > 0
#'     sig.counts <- sig.counts[which(sig.counts > 0)]
#'     ## get label heights:
#'     Y <- list()
#'     ## get average height of labels
#'     Y.avg <- max(h$counts)/4
#'     ## for bins with > 1 sig SNP, adjust height
#'     if(length(sig.counts)!=0){
#'       for(i in 1:length(sig.counts)){
#'         if(!.is.integer0(sig.counts[i])){
#'           if(sig.counts[i]==1){
#'             Y[[i]] <- Y.avg
#'           }else{
#'             ## divide up the space on the y-axis into
#'             ## increments (adding 1 to the denomenator s.t y-max not exceeded)
#'             increment <- max(h$counts)/(sig.counts[i]+1)
#'             Y[[i]] <- increment*c(1:sig.counts[i])
#'           }
#'         }
#'       } # end for loop
#'     }
#'     Y <- as.vector(unlist(Y))
#'
#'     if(plot.dist==TRUE){
#'
#'       ## plot histogram of correlations btw real SNPs and phenotype: ##
#'       plot(h, main="Distribution of SNP-trait correlations
#'            \n (with significant SNPs indicated if present)")
#'       ## ADD threshold line in red on x-axis where thresh hits...
#'       abline(v=thresh, col="grey", lwd=2, lty=2)
#'       ## label threshold line(?)
#'       text(x=thresh, y=(max(h$counts)),
#'            labels="significance threshold",
#'            col="grey", pos=4, font=4)
#'
#'       ## Only ask to draw arrows if at least 1 significant SNP:
#'       if(length(sig.snps) > 0){
#'         ## ADD arrows pointing from each label to
#'         ## position on X-axis:
#'         arrows(x0=X , y0=(Y-(max(h$counts)/50)) ,
#'                x1=X , y1=0 , col="blue", length=0.1, lwd=1)
#'         ## add annotation text labelling SNPs > threshold at their location on the x-axis:
#'         text(x=X, y=Y, labels=sig.snps,
#'              col="red", font=2, pos=4)
#'       }else{
#'         text(x=thresh, y=(max(h$counts)/4)*3,
#'              labels="no significant SNPs found",
#'              col="red", font=2, pos=4)
#'       }
#'     } # end plot.dist
#'   } # end test fisher
#'
#'   if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=FALSE)
#'
#' } # end plot.sig.snps
