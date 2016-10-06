

###################
## beeswarm.plot ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Beeswarm-and-Box-Plot.
#'
#' Wrapper combining the beeswarm and box plot functions from packages \code{beeswarm} and \code{ggplot2}.
#'
#' @param y A character string specifying the label of the (numeric) column
#' in data frame \code{df} to be plotted along the y-axis.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' ## load data
#' data(dist)
#' str(dist)
#'
#' ## basic use of fn
#' fn(arg1, arg2)
#'
#' #' ## more elaborate use of fn
#' fn(arg1, arg2)
#'
#' @import ggplot2
#' @importFrom beeswarm beeswarm
#'
#' @export

########################################################################

###################
## beeswarm.plot ##
###################

beeswarm.plot <- function(y="sensitivity", x="test", df, y.lab=NULL,
                          pt.size=4, x.text=FALSE){

  require(beeswarm)
  require(Hmisc)
  # library(plyr)
  # library(ggplot2)

  if(is.null(y.lab)) y.lab <- y

  ## Y ~ X ??
  fm <- as.formula(paste(y, x, sep=" ~ "))

  beeswarm.ori <- beeswarm::beeswarm(fm,
                                     data = df,
                                     #method="swarm", # swarm square hex center
                                     #priority="descending", ## ONLY for SWARM method...
                                     method="center", # swarm square hex center
                                     #priority="descending", ## ONLY for SWARM method...
                                     pwcol = eval(parse(text=x)),
                                     #col = myCol, ## to set w funky colours (INSTEAD of pwcol = test)
                                     ylim = c(-0.001,1), # otherwise ggplot can't plot ZERO values --> NAs
                                     las=2,
                                     cex=0.8,
                                     corral = "omit",
                                     do.plot = FALSE) # none gutter wrap omit
  # head(beeswarm)

  ######################################################
  ## Find and Replace OUTLIERS(' symbols in plot...): ##
  ######################################################
  outliers <- outlier.vals <- PCH <- list()

  if(!is.factor(df[,x])) df[,x] <- as.factor(df[,x])
  if(!all(as.character(beeswarm.ori$col) %in% levels(df[,x]))){
    foo <- beeswarm.ori$col
    foo <- levels(df[,x])[foo]
    beeswarm.ori$col <- factor(foo, levels=levels(df[,x]))
  }else{
    beeswarm.ori$col <- factor(beeswarm.ori$col)
  }

  noms <- as.character(levels(beeswarm.ori$col))

  ## FOR LOOP ##
  for(i in 1:length(noms)){
    #i <- 1
    # get vals for variable (and boxplot)
    val <- beeswarm.ori$y[which(beeswarm.ori$col==noms[i])]
    #boxplot(val, ylim=c(-0.001, 1))
    if(length(val) == 0){
      PCH[[i]] <- NULL
      outliers[[i]] <- NULL
    }else{
      PCH[[i]] <- rep(16, length(val)) # standard filled circle...

      ## get median
      M <- as.numeric(quantile(val, 0.5))
      # get lower 25 of box
      Q25 <- as.numeric(quantile(val, 0.25))
      # get upper 75 of box
      Q75 <- as.numeric(quantile(val, 0.75))
      # get box length
      box <- Q75-Q25

      if(box == 0) box <- 0.0000001

      # with a coef of 1.5 (the default for boxplots), identify outlying values
      outliers[[i]] <- c(which(val < Q25-(1.5*box)), which(val > Q75+(1.5*box)))
      # get values of outliers
      if(length(outliers[[i]]) > 0){
        outlier.vals[[i]] <- val[outliers[[i]]]
        PCH[[i]] <- replace(PCH[[i]], outliers[[i]], 17) # replace with triangle...
      }else{
        outlier.vals[[i]] <- NULL
      }
    }
  } # end for loop

  #outliers
  PCH <- as.vector(unlist(PCH))
  # PCH


  #########################################################################################################
  ######################
  ## plots, layers... ##
  ######################

  if(x.text == FALSE){

    ################
    ## NO X-TEXT: ##
    ################

    beeswarm.plot1 <- ggplot(beeswarm.ori, aes(x, y)) +
      xlab("") +
      guides(fill=FALSE) +
      scale_x_discrete(drop=FALSE) +
      scale_y_continuous(y.lab, limits=c(0,1))  # expression("char")

    beeswarm.plot2 <- beeswarm.plot1 +
      guides(fill=FALSE) +
      geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
                   outlier.shape = 17,
                   outlier.size=0) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y = element_text(size=13),
            # axis.title.y=element_text(size=18),
            axis.title.y=element_blank(),
            legend.position="none")

    beeswarm.plot3 <- beeswarm.plot2 +
      geom_point(data=beeswarm.ori, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
      guides(fill=FALSE) +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y = element_text(size=13),
            # axis.title.y=element_text(size=18),
            axis.title.y=element_blank(),
            legend.position="none")

    beeswarm.plot4 <- beeswarm.plot3 +
      guides(fill=FALSE) +
      geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
                   outlier.shape = 17,
                   outlier.size=0) +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y = element_text(size=13),
            # axis.title.y=element_text(size=18),
            axis.title.y=element_blank(),
            legend.position="none")
  }else{

    ##################
    ## WITH X-TEXT: ##
    ##################

    #############################
    ## NUMERIC x-axis labels?? ##
    #############################
    if(all.is.numeric(levels(df[, x]))){
      beeswarm.plot1 <- ggplot(beeswarm.ori, aes(x, y)) +
        xlab("") +
        guides(fill=FALSE) +
        scale_x_discrete(drop=FALSE) +
        scale_y_continuous(y.lab, limits=c(0,1)) # expression("char")

      beeswarm.plot2 <- beeswarm.plot1 +
        guides(fill=FALSE) +
        geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
                     outlier.shape = 17,
                     outlier.size=0) +
        theme(axis.text.x = element_text(size=13),
              # axis.text.x = element_text(angle=35, hjust=1, vjust=0.95, size=10),
              axis.text.y = element_text(size=13),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              legend.position="none")

      beeswarm.plot3 <- beeswarm.plot2 +
        geom_point(data=beeswarm.ori, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
        guides(fill=FALSE) +
        theme(axis.text.x = element_text(size=13),
              # axis.text.x = element_text(angle=35, hjust=1, vjust=0.95, size=10),
              axis.text.y = element_text(size=13),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              legend.position="none")

      beeswarm.plot4 <- beeswarm.plot3 +
        guides(fill=FALSE) +
        geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
                     outlier.shape = 17,
                     outlier.size=0) +
        theme(axis.text.x = element_text(size=13),
              # axis.text.x = element_text(angle=35, hjust=1, vjust=0.95, size=10),
              axis.text.y = element_text(size=13),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              legend.position="none")
    }else{

      ######################################
      ## CHARACTER STRING x-axis labels?? ##
      ######################################
      beeswarm.plot1 <- ggplot(beeswarm.ori, aes(x, y)) +
        xlab("") +
        guides(fill=FALSE) +
        scale_x_discrete(drop=FALSE) +
        scale_y_continuous(y.lab, limits=c(0,1)) # expression("char")

      beeswarm.plot2 <- beeswarm.plot1 +
        guides(fill=FALSE) +
        geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
                     outlier.shape = 17,
                     outlier.size=0) +
        theme(axis.text.x = element_text(angle=35, hjust=1, vjust=0.95, size=10),
              # axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              legend.position="none")

      beeswarm.plot3 <- beeswarm.plot2 +
        geom_point(data=beeswarm.ori, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
        guides(fill=FALSE) +
        theme(axis.text.x = element_text(angle=35, hjust=1, vjust=0.95, size=10),
              # axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              legend.position="none")

      beeswarm.plot4 <- beeswarm.plot3 +
        guides(fill=FALSE) +
        geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
                     outlier.shape = 17,
                     outlier.size=0) +
        theme(axis.text.x = element_text(angle=35, hjust=1, vjust=0.95, size=10),
              # axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              legend.position="none")
    }
  }


  ## PRINT PLOT ##
  plot(beeswarm.plot4)

} # end beeswarm.plot
