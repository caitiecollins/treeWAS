
#################
## gap.boxplot ##
#################

## from pkg plotrix on github


## args:
# x <- with(dat, boxplot(FPR ~ test)) # , ylim=c(0,1)

# x <- boxplot(dat$FPR ~ dat$test)
# x <- dat$FPR ~ dat$test
#
# formula <- FPR ~ test
# data <- dat
#
# col=transp(rainbow(32), 0.4)
#
# gap=list(top=c(NA, NA), bottom=c(0.0000, 0.0000))
#
# range=1.5
# width=NULL
# varwidth=FALSE
# notch=FALSE
# outline=TRUE
# names=NULL
# xlim=NA
# ylim=c(0,1)
# plot=TRUE
# border=par("fg")
#
# log=""
# axis.labels=NULL
# axes=TRUE
# pars=list(boxwex=0.8,staplewex=0.5,outwex=0.5)
# horizontal=FALSE
# add=FALSE
# at=NULL
# main=NULL


gap.boxplot<-function (x,...,gap=list(top=c(NA,NA),bottom=c(NA,NA)),
                       range=1.5,width=NULL,varwidth=FALSE,notch=FALSE,outline=TRUE,names,
                       xlim=NA,ylim=NA,plot=TRUE,border=par("fg"),col=NULL,log="",
                       axis.labels=NULL,axes=TRUE,pars=list(boxwex=0.8,staplewex=0.5,outwex=0.5),
                       horizontal=FALSE,add=FALSE,at=NULL,main=NULL) {

  if(!is.na(gap$top[1]))
    if(gap$top[1] > gap$top[2]) gap$top<-rev(gap$top)

    if(!is.na(gap$bottom[1]))
      if(gap$bottom[1] > gap$bottom[2]) gap$bottom<-rev(gap$bottom)

      if(is.na(ylim[1])) {
        bxpt<-boxplot(x, range=range,plot=FALSE)
        ylim<-range(c(bxpt$stats,bxpt$out))
      }else{
        bxpt<-boxplot(x, ylim=ylim,range=range,plot=FALSE)
        }
      bxgap<-bxpt


      if(!is.na(gap$top[1])) {
        bxgap$stats[bxgap$stats > gap$top[1] & bxgap$stats < gap$top[2]]<-NA
        if(any(is.na(bxgap$stats)))
          stop("gap cannot include the median, interquartiles or the staples")
        topdiff<-diff(gap$top)
        bxgap$stats[bxgap$stats > gap$top[2]]<-
          bxgap$stats[bxgap$stats > gap$top[2]]-topdiff
        intopgap<-bxgap$out > gap$top[1] & bxgap$out < gap$top[2]
        bxgap$out[intopgap]<-NA
        abovetop<-which(bxgap$out > gap$top[2])
        bxgap$out[abovetop]<-bxgap$out[abovetop]-topdiff
        rangetop<-gap$top[1]
        ylim[2]<-ylim[2]-topdiff
      }else{
        rangetop<-ylim[2]
      }


      if(!is.na(gap$bottom[1])) {
        bxgap$stats[bxgap$stats > gap$bottom[1] & bxgap$stats < gap$bottom[2]]<-NA
        if(any(is.na(bxgap$stats)))
          stop("gap cannot include the median, interquartiles or the staples")
        bottomdiff<-diff(gap$bottom)
        bxgap$stats[bxgap$stats < gap$bottom[1]]<-
          bxgap$stats[bxgap$stats < gap$bottom[1]]+bottomdiff
        bxgap$out[bxgap$out > gap$bottom[1] & bxgap$out < gap$bottom[2]] <- NA
        belowbottom<-which(bxgap$out < gap$bottom[1])
        bxgap$out[belowbottom]<-bxgap$out[belowbottom]+bottomdiff
        rangebottom<-gap$bottom[2]
        ylim[1]<-ylim[1]+bottomdiff
      }else{
        rangebottom<-ylim[1]
      }

      if(any(is.na(bxgap$out)))
        warning("At least one outlier falls into a gap")
      nboxes<-dim(bxgap$stats)[2]
      if(is.na(xlim[1])) {
        xlim<-c(0.5,nboxes+0.5)
        at<-1:nboxes
      }
      bxgap$group<-at
      plot(0,xlim=xlim,ylim=ylim,type="n",axes=FALSE,xlab="",ylab="",main=main)
      plotlim<-par("usr")
      box()

      if(axes) axis(1,labels=bxpt$names,at=at)
      midticks<-pretty(c(rangebottom,rangetop))
      if(axes) axis(2,at=midticks[midticks > rangebottom & midticks < rangetop])
      if(is.null(width)) width<-pars$boxwex
      rect(at-width/2,bxgap$stats[2,],at+width/2,
           bxgap$stats[4,],border=border,col=col)

      if(notch) {
        ymult<-getYmult()
        if(is.null(col)){
          boxcol<-"white"
        }else{
          boxcol<-col
        }
        rect(at-width/1.95,bxgap$conf[1,],at+width/1.95,
             bxgap$conf[2,],border=NA,col=boxcol)
        insets<-(bxgap$conf[2,]-bxgap$conf[1,])*pars$boxwex/ymult
        median.left<-((at-width/2)+insets)
        median.right<-((at+width/2)-insets)

        # display the notches
        segments(at-width/2,bxgap$conf[1,],median.left,bxgap$stats[3,],col=border)
        segments(at-width/2,bxgap$conf[2,],median.left,bxgap$stats[3,],col=border)
        segments(median.right,bxgap$stats[3,],at+width/2,bxgap$conf[1,],col=border)
        segments(median.right,bxgap$stats[3,],at+width/2,bxgap$conf[2,],col=border)
      }else{
        median.left<-at-width/2
        median.right<-at+width/2
      }


      # draw the median line
      segments(median.left,bxgap$stats[3,],median.right,bxgap$stats[3,],
               lwd=2,col=border)
      segments(at,bxgap$stats[1,],at,bxgap$stats[2,],lty=2,col=border)
      segments(at,bxgap$stats[4,],at,bxgap$stats[5,],lty=2,col=border)
      segments(at-pars$staplewex*width/2,bxgap$stats[1,],
               at+pars$staplewex*width/2,bxgap$stats[1,],col = border)
      segments(at-pars$staplewex*width/2,bxgap$stats[5,],
               at+pars$staplewex*width/2,bxgap$stats[5,],col=border)
      if(!is.na(gap$top[1])){
        topadjust<-diff(gap$top)
      }else{
        topadjust<-0
      }
      if(!is.na(gap$bottom[1])){
        bottomadjust<-diff(gap$bottom)
      }else{
        bottomadjust<-0
      }
      if(!is.null(axis.labels)) axis(2,labels=axis.labels,
                                     at=c(axis.labels[1]+bottomadjust,axis.labels[2]-topadjust))
      if(!is.na(gap$top[1])) axis.break(2,gap$top[1],style="gap")
      if(!is.na(gap$bottom[1]))
        axis.break(2,gap$bottom[2]-diff(plotlim[3:4])*0.02,style="gap")
      if(length(bxgap$group) == length(bxgap$out)) points(bxgap$group,bxgap$out)
      invisible(bxgap)
}
