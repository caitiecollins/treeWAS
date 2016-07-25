
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


plot.sig.snps <- function(corr.dat, corr.sim,
                          sig.corrs, sig.snps,
                          sig.thresh, test,
                          plot.null.dist=TRUE, plot.dist=FALSE){

  thresh <- sig.thresh

  if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=TRUE)

  ## plot.null.dist #####################################

  ###################################################
  ## Add SNP annotations to histogram of null dist ##
  ###################################################


  ##############
  ## TERMINAL ##
  ##############
  if(test=="terminal"){
    ## plot.null.dist ##

    ###########################################################
    ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
    ###########################################################
    ## plot correlations btw simulated SNPs and phenotype:
    h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)

    ## get alternate (null dist) label heights:
    X <- thresh
    if(length(sig.corrs) > 0) X <- sig.corrs
    Y.null <- runif(n=length(sig.snps), min=0.00001,
                    max=max(h.null$counts))


    if(plot.null.dist==TRUE){
      ## if the true correlation value for SNP i is >
      ## max bin, then extend the x-axis of the plot
      ## to accommodate annotation:
      if(max(X) > max(h.null$breaks)){
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Terminal Correlation Scores
             \n (with significant SNPs indicated)",
             xlab="Terminal Correlation Score",
             xlim=c(min(h.null$breaks), max(X)+.05))
      }else{
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## WITHOUT EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Terminal Correlation Scores
             \n (with significant SNPs indicated)",
             xlab="Terminal Correlation Score"
        )
      }

      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h.null$counts)),
           labels="significance threshold", pos=2,
           col="grey", font=4)
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to
        ## position on X-axis:
        arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
               x1=X , y1=0 , col="blue",
               length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y.null, labels=sig.snps,
             col="red", font=2, pos=2)
      }else{
        text(x=thresh, y=(max(h.null$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=2)
      }
    } # end plot.null.dist
  } # end test Terminal


  ##################
  ## SIMULTANEOUS ##
  ##################
  if(test=="simultaneous"){
    ## plot.null.dist ##

    ###########################################################
    ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
    ###########################################################
    ## plot correlations btw simulated SNPs and phenotype:
    h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)

    ## get alternate (null dist) label heights:
    X <- thresh
    if(length(sig.corrs) > 0) X <- sig.corrs
    Y.null <- runif(n=length(sig.snps), min=0.00001,
                    max=max(h.null$counts))


    if(plot.null.dist==TRUE){
      ## if the true correlation value for SNP i is >
      ## max bin, then extend the x-axis of the plot
      ## to accommodate annotation:
      if(max(X) > max(h.null$breaks)){
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Simultaneous Scores
             \n (with significant SNPs indicated)",
             xlab="Simultaneous Score",
             xlim=c(min(h.null$breaks), max(X)+.05))
      }else{
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## WITHOUT EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Simultaneous Scores
             \n (with significant SNPs indicated)",
             xlab="Simultaneous Score"
        )
      }

      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h.null$counts)),
           labels="significance threshold", pos=2,
           col="grey", font=4)
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to
        ## position on X-axis:
        arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
               x1=X , y1=0 , col="blue",
               length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y.null, labels=sig.snps,
             col="red", font=2, pos=2)
      }else{
        text(x=thresh, y=(max(h.null$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=2)
      }
    } # end plot.null.dist
  } # end test Simultaneous

  ################
  ## SUBSEQUENT ##
  ################
  if(test=="subsequent"){
    ## plot.null.dist ##

    ###########################################################
    ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
    ###########################################################
    ## plot correlations btw simulated SNPs and phenotype:
    h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)

    ## get alternate (null dist) label heights:
    X <- thresh
    if(length(sig.corrs) > 0) X <- sig.corrs
    Y.null <- runif(n=length(sig.snps), min=0.00001,
                    max=max(h.null$counts))


    if(plot.null.dist==TRUE){
      ## if the true correlation value for SNP i is >
      ## max bin, then extend the x-axis of the plot
      ## to accommodate annotation:
      if(max(X) > max(h.null$breaks)){
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Subsequent Scores
             \n (with significant SNPs indicated)",
             xlab="Subsequent Score",
             xlim=c(min(h.null$breaks), max(X)+.05))
      }else{
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## WITHOUT EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Subsequent Scores
             \n (with significant SNPs indicated)",
             xlab="Subsequent Score"
        )
      }

      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h.null$counts)),
           labels="significance threshold", pos=2,
           col="grey", font=4)
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to
        ## position on X-axis:
        arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
               x1=X , y1=0 , col="blue",
               length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y.null, labels=sig.snps,
             col="red", font=2, pos=2)
      }else{
        text(x=thresh, y=(max(h.null$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=2)
      }
    } # end plot.null.dist
  } # end test Subsequent


  #################
  ## CORRELATION ##
  #################
  if(test=="cor"){
    ## plot.null.dist ##

    ###########################################################
    ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
    ###########################################################
    ## plot correlations btw simulated SNPs and phenotype:
    h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)

    ## get alternate (null dist) label heights:
    X <- thresh
    if(length(sig.corrs) > 0) X <- sig.corrs
    Y.null <- runif(n=length(sig.snps), min=0.00001,
                    max=max(h.null$counts))


    if(plot.null.dist==TRUE){

      ## if the true correlation value for SNP i is >
      ## max bin, then extend the x-axis of the plot
      ## to accommodate annotation:
      if(max(X) > max(h.null$breaks)){
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of correlations
             \n (with significant SNPs indicated)",
             xlab="Correlation",
             xlim=c(min(h.null$breaks), max(X)+.05))
      }else{
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## WITHOUT EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of correlations
             \n (with significant SNPs indicated)",
             xlab="Correlation"
        )
      }

      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h.null$counts)),
           labels="significance threshold", pos=2,
           col="grey", font=4)
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to
        ## position on X-axis:
        arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
               x1=X , y1=0 , col="blue", length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y.null, labels=sig.snps,
             col="red", font=2)
      }else{
        text(x=thresh, y=(max(h.null$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=2)
      }
    } # end plot.null.dist
    } # end test cor


  #########################
  ## FISHER'S EXACT TEST ##
  #########################
  if(test=="fisher"){
    ## plot.null.dist ##

    ###########################################################
    ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
    ###########################################################
    ## plot correlations btw simulated SNPs and phenotype:
    h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)

    ## get alternate (null dist) label heights:
    X <- thresh
    if(length(sig.corrs) > 0) X <- sig.corrs
    Y.null <- runif(n=length(sig.snps), min=0.00001,
                    max=max(h.null$counts))

    if(plot.null.dist==TRUE){
      ## if the true correlation value for SNP i is <
      ## min bin, then extend the x-axis of the plot to
      ## accommodate annotation:
      if(min(X) < min(h.null$breaks)){
        ## plot histogram of correlations btw real
        ## SNPs and phenotype: ##
        ## EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Fisher's exact test p-values
             \n (with significant SNPs indicated)",
             xlab="p-value",
             xlim=c((min(X)-.05), max(h.null$breaks)))
      }else{
        ## plot histogram of correlations btw real SNPs and phenotype: ##
        ## WITHOUT EXTENDING THE X-AXIS
        plot(h.null,
             main="Null distribution of Fisher's exact test p-values
             \n (with significant SNPs indicated)",
             xlab="p-value"
        )
      }

      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h.null$counts)),
           labels="significance threshold", pos=4,
           col="grey", font=4)
      ## only ask to draw arrows if sig snps exist
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to
        ## position on X-axis:
        arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) ,
               x1=X , y1=0 , col="blue", length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y.null, labels=sig.snps,
             col="red", font=2, pos=4)
      }else{
        text(x=thresh, y=(max(h.null$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=4)
      }
    } # end plot.null.dist
  } # end test fisher


  ## plot.dist ###################

  #################################################################
  ## plot histogram of correlations btw real SNPs and phenotype: ##
  #################################################################

  ###########
  ## SCORE ##
  ###########
  if(test=="score"){
    ## plot.dist ##

    #################################################################
    ## plot histogram of correlations btw real SNPs and phenotype: ##
    #################################################################
    ## get histogram of correlations btw real SNPs and phenotype:
    h <- hist(corr.dat, plot=FALSE)
    ## get X and Y coords for labelling positions of all sig SNPs
    X <- sig.corrs
    ## get the number of sig SNPs in each bin of the histogram
    ## counting correlations==upper limit of each bin
    ## as falling within that bin...
    sig.counts <- sapply(c(1:(length(h$breaks)-1)),
                         function(e) length(X[which(X[which(X <=
                                                              h$breaks[(e+1)])] > h$breaks[e])]))
    ## keep only counts > 0
    sig.counts <- sig.counts[which(sig.counts > 0)]
    ## get label heights:
    Y <- list()
    ## get average height of labels
    Y.avg <- max(h$counts)/4
    ## for bins with > 1 sig SNP, adjust height
    if(length(sig.counts)!=0){
      for(i in 1:length(sig.counts)){
        if(!.is.integer0(sig.counts[i])){
          if(sig.counts[i]==1){
            Y[[i]] <- Y.avg
          }else{
            ## divide up the space on the y-axis into
            ## increments (adding 1 to the denomenator
            ## s.t y-max not exceeded)
            increment <- max(h$counts)/(sig.counts[i]+1)
            Y[[i]] <- increment*c(1:sig.counts[i])
          }
        }
      } # end for loop
    }
    Y <- as.vector(unlist(Y))

    if(plot.dist==TRUE){

      ## plot histogram of correlations btw real SNPs and phenotype: ##
      plot(h, main="Distribution of SNP-trait correlations
           \n (with significant SNPs indicated if present)")
      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h$counts)),
           labels="significance threshold",
           col="grey", pos=2, font=4)

      ## Only ask to draw arrows if at least 1 significant SNP:
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to position on X-axis:
        arrows(x0=X , y0=(Y-(max(h$counts)/50)) ,
               x1=X , y1=0 ,
               col="blue", length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y, labels=sig.snps,
             col="red", font=2, pos=2)
      }else{
        text(x=thresh, y=(max(h$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=2)
      }
    } # end plot.dist
  } # end test score

  #################
  ## CORRELATION ##
  #################
  if(test=="cor"){
    ## plot.dist ##

    #################################################################
    ## plot histogram of correlations btw real SNPs and phenotype: ##
    #################################################################
    ## get histogram of correlations btw real SNPs and phenotype:
    h <- hist(corr.dat, plot=FALSE)
    ## get X and Y coords for labelling positions of all sig SNPs
    X <- sig.corrs
    ## get the number of sig SNPs in each bin of the histogram
    ## counting correlations==upper limit of each bin as falling
    ## within that bin...
    sig.counts <- sapply(c(1:(length(h$breaks)-1)),
                         function(e) length(X[which(X[which(X <=
                                                              h$breaks[(e+1)])] > h$breaks[e])]))
    ## keep only counts > 0
    sig.counts <- sig.counts[which(sig.counts > 0)]
    ## get label heights:
    Y <- list()
    ## get average height of labels
    Y.avg <- max(h$counts)/4
    ## for bins with > 1 sig SNP, adjust height
    if(length(sig.counts)!=0){
      for(i in 1:length(sig.counts)){
        if(!.is.integer0(sig.counts[i])){
          if(sig.counts[i]==1){
            Y[[i]] <- Y.avg
          }else{
            ## divide up the space on the y-axis into
            ## increments (adding 1 to the denomenator
            ## s.t y-max not exceeded)
            increment <- max(h$counts)/(sig.counts[i]+1)
            Y[[i]] <- increment*c(1:sig.counts[i])
          }
        }
      } # end for loop
    }
    Y <- as.vector(unlist(Y))

    if(plot.dist==TRUE){

      ## plot histogram of correlations btw real
      ## SNPs and phenotype: ##
      plot(h, main="Distribution of SNP-trait
           correlations \n (with significant SNPs indicated if present)")
      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h$counts)),
           labels="significance threshold",
           col="grey", pos=2, font=4)

      ## Only ask to draw arrows if at least 1 significant SNP:
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each
        ## label to position on X-axis:
        arrows(x0=X , y0=(Y-(max(h$counts)/50)) ,
               x1=X , y1=0 , col="blue", length=0.1, lwd=1)
        ## add annotation text labelling SNPs >
        ## threshold at their location on the x-axis:
        text(x=X, y=Y, labels=sig.snps, col="red", font=2)
      }else{
        text(x=thresh, y=(max(h$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=2)
      }
    } # end plot.dist
  } # end test cor


  #########################
  ## FISHER'S EXACT TEST ##
  #########################
  if(test=="fisher"){
    ## plot.dist ##

    ############################################################
    ## plot histogram of correlations btw real SNPs and phen: ##
    ############################################################

    ## get histogram of correlations btw real SNPs and phenotype:
    h <- hist(corr.dat, plot=FALSE)
    ## get X and Y coords for labelling positions of all sig SNPs
    X <- thresh
    if(length(sig.corrs) > 0) X <- sig.corrs
    ## get the number of sig SNPs in each bin of the histogram
    ## counting correlations==upper limit of each bin as
    ## falling within that bin...
    sig.counts <- sapply(c(1:(length(h$breaks)-1)),
                         function(e) length(X[which(X[which(X <=
                                                              h$breaks[(e+1)])] > h$breaks[e])]))
    ## keep only counts > 0
    sig.counts <- sig.counts[which(sig.counts > 0)]
    ## get label heights:
    Y <- list()
    ## get average height of labels
    Y.avg <- max(h$counts)/4
    ## for bins with > 1 sig SNP, adjust height
    if(length(sig.counts)!=0){
      for(i in 1:length(sig.counts)){
        if(!.is.integer0(sig.counts[i])){
          if(sig.counts[i]==1){
            Y[[i]] <- Y.avg
          }else{
            ## divide up the space on the y-axis into
            ## increments (adding 1 to the denomenator s.t y-max not exceeded)
            increment <- max(h$counts)/(sig.counts[i]+1)
            Y[[i]] <- increment*c(1:sig.counts[i])
          }
        }
      } # end for loop
    }
    Y <- as.vector(unlist(Y))

    if(plot.dist==TRUE){

      ## plot histogram of correlations btw real SNPs and phenotype: ##
      plot(h, main="Distribution of SNP-trait correlations
           \n (with significant SNPs indicated if present)")
      ## ADD threshold line in red on x-axis where thresh hits...
      abline(v=thresh, col="grey", lwd=2, lty=2)
      ## label threshold line(?)
      text(x=thresh, y=(max(h$counts)),
           labels="significance threshold",
           col="grey", pos=4, font=4)

      ## Only ask to draw arrows if at least 1 significant SNP:
      if(length(sig.snps) > 0){
        ## ADD arrows pointing from each label to
        ## position on X-axis:
        arrows(x0=X , y0=(Y-(max(h$counts)/50)) ,
               x1=X , y1=0 , col="blue", length=0.1, lwd=1)
        ## add annotation text labelling SNPs > threshold at their location on the x-axis:
        text(x=X, y=Y, labels=sig.snps,
             col="red", font=2, pos=4)
      }else{
        text(x=thresh, y=(max(h$counts)/4)*3,
             labels="no significant SNPs found",
             col="red", font=2, pos=4)
      }
    } # end plot.dist
  } # end test fisher

  if(plot.null.dist==TRUE & plot.dist==TRUE) par(ask=FALSE)

} # end plot.sig.snps
