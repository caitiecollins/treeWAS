
##################
## Pagel 94 fns ##
##################


# #' #######################################################################
# #'
# #' ##################
# #' # DOCUMENTATION ##
# #' ##################
# #'
# #' Short one-phrase description.
# #'
# #' Longer proper discription of function...
# #'
# #' @param tree A phylo object.
# #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
# #' @examples
# #'
# #' @import phytools geiger ape phangorn
# #' @export
# #'
# #' #######################################################################
#
# get.n.subs <- function(n){
#   Q <- matrix(c(-1/n, 1/n, 1/n, -1/n), 2, 2)
#   rownames(Q) <- colnames(Q) <- letters[1:2]
#   tt1 <- sim.history(tree, Q)
#   ## Done simulation(s).
#
#   maps <- tt1$maps
#   n.subs <- length(as.vector(unlist(maps)))-length(maps)
#   return(n.subs)
# } # end get.n.subs




###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' @import phytools geiger ape phangorn
#' @export

########################################################################

##########
## mexp ##
##########
# library(rmutil) #??
mexp <- function(x, type="spectral decomposition", t=1, n=20, k=3){

  if(!is.matrix(x))stop("x must be a matrix")
  if(length(dim(x))!=2)stop("x must be a two dimensional matrix")
  if(dim(x)[1]!=dim(x)[2])stop("x must be a square matrix")

  type <- match.arg(type,c("spectral decomposition","series approximation"))
  d <- ncol(x)

  if(type=="spectral decomposition"){
    z <- eigen(t*x,sym=F)
    p <- z$vectors%*%diag(exp(z$values))%*%solve(z$vectors)

  }else{
    xx <- x*t/2^k
    p <- diag(d)
    q <- p
    for(r in 1:n){
      q <- xx%*%q/r
      p <- p+q
    }
    for(i in 1:k)
      p <- p%*%p
  }
  p
} # end mexp

# test1 <- t(matrix(c(
#   4, 2, 0,
#   1, 4, 1,
#   1, 1, 4), 3, 3))
#
# library(PSM)
# matexp(test1) # PSM fn
# mexp(test1) # my fn (modified from online sources...)
# identical(matexp(test1), mexp(test1)) # FALSE ...
# identical(round(matexp(test1), 4), round(mexp(test1), 4)) # TRUE!

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' @import phytools geiger ape phangorn
#' @export

########################################################################

#########
## Q2P ##
#########
Q2P <- function(Q, type=NULL){

  #requires fn mexp

  ## GET PROBABILITY MATRIX:
  P <- mexp(0.05 * Q)
  row.names(P) <- row.names(Q)
  colnames(P) <- colnames(Q)

  ## SPECIFIC COMPONENTS FOR fitPagel x|y VARIABLES:
  ## (ie. if Q is a rate matrix for
  ## a single variable that combines 2 variables as x|y,
  ## for example, as used/returned by fitPagel)
  if(!is.null(type)){
    ## INDEPENDENT RATES
    if(type=="independent"){
      ## CHECK: if type=="idependent"
      ## Make sure pairs of cells add up to same probabilities.
      #       if(
      #         any(c(
      #         !identical(round(sum(P[1,1:2]), 10), round(sum(P[2,1:2]), 10)),
      #         !identical(round(sum(P[1,3:4]), 10), round(sum(P[2,3:4]), 10)),
      #         !identical(round(sum(P[3,1:2]), 10), round(sum(P[4,1:2]), 10)),
      #         !identical(round(sum(P[3,3:4]), 10), round(sum(P[4,3:4]), 10))))
      #         ){
      #         warning("If type == independent,
      #                  pairs of cells in quadrants of
      #                  the probability matrix
      #                  should sum to the same thing.
      #                  They do not.")
      #       }

      P00 <- sum(P[1, 1:2])
      P01 <- sum(P[1, 3:4])
      P11 <- sum(P[3, 1:2])
      P10 <- sum(P[3, 3:4])

      probs <- c(P00, P01, P10, P11)
      names(probs) <- c("P00", "P01", "P10", "P11")
    }else{
      ## DEPENDENT RATES
      ## NOTE-- Only set up for 4x4 matrices right now.
      ## ... Could generalise later if useful to user.
      if(ncol(P) == 4){
        ## remove diagonals
        probs <- as.vector(unlist(t(P)))[c(2,3,5,8,9,12,14,15)]
        ## reorder in terms of var1 first
        probs <- probs[c(1,3,2,4,5,7,6,8)]

        v1 <- as.character(c(0, 0, 1, 1))
        v2 <- as.character(c(0, 1, 0, 1))

        noms <- list()
        noms[[1]] <- paste("P", v1[1], v1[2], "|", v2[1], v2[2], sep="")
        noms[[2]] <- paste("P", v1[2], v1[1], "|", v2[2], v2[1], sep="")
        noms[[3]] <- paste("P", v1[1], v1[3], "|", v2[1], v2[3], sep="")
        noms[[4]] <- paste("P", v1[2], v1[4], "|", v2[2], v2[4], sep="")
        noms[[5]] <- paste("P", v1[3], v1[1], "|", v2[3], v2[1], sep="")
        noms[[6]] <- paste("P", v1[4], v1[2], "|", v2[4], v2[2], sep="")
        noms[[7]] <- paste("P", v1[3], v1[4], "|", v2[3], v2[4], sep="")
        noms[[8]] <- paste("P", v1[4], v1[3], "|", v2[4], v2[3], sep="")
        noms <- as.vector(unlist(noms))

        names(probs) <- noms
      }
    }
    ## return both P matrix and individual probs vector:
    out <- list(P, probs)
    names(out) <- c("P", "probs")
  }else{

    ## if not a fitPagel xy variable, just return P, the whole prob matrix
    out <- list(P)
    names(out) <- c("P")
  }

  return(out)

} # end Q2P



########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#' @param Q A matrix containing transition rates (as returned )
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' @import phytools
#' @export

########################################################################

################
## get.n.subs ##
################

## Get n.subs for a given SNP
## from rate.mat Q (fitPagel iQ)
## & tree height + branch lengths

## ASSUMPTION (guess):
## Rates from fitPagel/ace in Q are in units of
## number of substitutions (in a given direction)
## per tree height.
## Therefore, the expected value of the n.subs in the tree
## is the sum(rates)*(sum(branch.lengths)/tree.height)

get.n.subs <- function(Q, tree){

  require(phytools)

  ## get height of tree:
  H <- nodeHeights(tree)
  tree.height <- max(H)

  ## get branch lengths in tree:
  branch.lengths <- tree$edge.length

  ## get the number of units of time contained in the tree:
  n.unit <- sum(branch.lengths)/tree.height
  # n.unit # 6.88 (for snps[,5])

  ## estimate the expected number of substitutions
  ## given rates and n.unit time:
  if(ncol(Q) == 4){
    rate01 <- Q[1,3]
    rate10 <- Q[3,1]
    n.subs01 <- rate01*n.unit
    n.subs10 <- rate10*n.unit
    n.subs <- c(n.subs01, n.subs10)
    names(n.subs) <- c("n.subs01", "n.subs10")
  }else{
    n.subs <- Q*n.unit
  }

  return(n.subs)

} # end get.n.subs



###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

##################################################################################
## FUNCTIONS FROM: ###############################################################
##################################################################################
## http://blog.phytools.org/2014/12/r-function-for-pagels-1994-correlation.html ##
## Posted here: http://www.phytools.org/fitPagel/v0.1/fitPagel.R #################
##################################################################################

## function fits Pagel '94 model of correlated evolution of two binary characters
## uses ape::ace or geiger::fitDiscrete internally
## written by Liam J. Revell 2014

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' @import phytools geiger ape phangorn
#' @export

########################################################################

fitPagel <- function(tree,x,y, equal=FALSE, ...){
  if(hasArg(method)){
    method <- list(...)$method
  }else{
    method <- "ace"
  }
  if(method=="fitDiscrete"){
    chk <- .check.pkg("geiger")
    if(!chk){
      cat("  method = \"fitDiscrete\" requires the package \"geiger\"\n")
      cat("  Defaulting to method = \"ace\"\n\n")
      method <- "ace"
      fitDiscrete <- function(...) NULL
    }
  }
  noms.x <- names(x)
  if(!is.factor(x)) x <- as.factor(x)
  levels.x <- levels(x)
  names(x) <- noms.x
  if(!is.factor(y)) y <- as.factor(y)
  levels.y <- levels(y)
  y <- y[names(x)]
  if(length(levels.x)!=2||length(levels.y)!=2)
    stop("Only binary characters for x & y currently permitted.")
  if(any(!levels.x %in% levels.y)){
    ## convert both to binary variables:
    warning("Levels of x != levels of y;
            converting both to binary variables.")
    x <- as.numeric(x)-1
    y <- as.numeric(y)-1
    if(!is.factor(x)) x <- as.factor(x)
    levels.x <- levels(x)
    if(!is.factor(y)) y <- as.factor(y)
    levels.y <- levels(y)
    y <- y[names(x)]
  }
  xy <- setNames(factor(paste(x,y,sep="|"),
                      levels=sapply(levels.x,paste,levels.y,sep="|")),
               names(x))

  ## fit independent model
  if(equal == FALSE){
    ## original (allows for diff rates fwd & bckwd)
    iQ <- matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
  }
  if(equal == TRUE){
    ## w EQUAL rates (for SNP only???) ## NOT INVERTIBLE
    #### iQ <- matrix(c(0,1,2,0,3,0,0,2,2,0,0,1,0,2,3,0),4,4,byrow=TRUE)

    ## w EQUAL rates (for SNP & PHEN both)
    iQ <- matrix(c(0,1,2,0,1,0,0,2,2,0,0,1,0,2,1,0),4,4,byrow=TRUE)
  }
  rownames(iQ) <- colnames(iQ) <- levels(xy)
  fit.iQ <- if(method=="fitDiscrete") fitDiscrete(tree,xy,model=iQ) else ace(xy,tree,type="discrete",model=iQ)

  ## fit dependendent model
  if(equal ==FALSE){
    dQ <- matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
  }
  if(equal == TRUE){
    ## w EQUAL rates (for SNP only???) ## NOT INVERTIBLE
    ### dQ <- matrix(c(0,1,2,0,3,0,0,4,2,0,0,5,0,4,6,0),4,4,byrow=TRUE)

    ## w EQUAL rates (for SNP & PHEN both)
    dQ <- matrix(c(0,1,2,0,1,0,0,3,2,0,0,4,0,3,4,0),4,4,byrow=TRUE)
  }
  rownames(dQ) <- colnames(dQ) <- levels(xy)
  if(method=="fitDiscrete"){
    fit.dQ  <- fitDiscrete(tree,xy,model=dQ)
  }else{
    fit.dQ  <- try(ace(xy,tree,type="discrete",model=dQ), silent=TRUE)
    if(class(fit.dQ) == "try-error") fit.dQ <- NULL
  }

  ## back translate independent model
  if(method=="fitDiscrete"){
    iQ <- geiger:::.Qmatrix.from.gfit(fit.iQ)
  }else{
    I <- fit.iQ$index.matrix
    I[I==0] <- NA
    iQ <- apply(I,2,function(i,x) x[i],x=fit.iQ$rates)
    iQ[is.na(iQ)] <- 0
    diag(iQ) <- -rowSums(iQ)
    rownames(iQ) <- colnames(iQ)
  }

  ## dependent model
  if(!is.null(fit.dQ)){
    if(method=="fitDiscrete"){
      dQ <- geiger:::.Qmatrix.from.gfit(fit.dQ)
    }else{
      I <- fit.dQ$index.matrix
      I[I==0] <- NA
      dQ <- apply(I,2,function(i,x) x[i],x=fit.dQ$rates)
      dQ[is.na(dQ)] <- 0
      diag(dQ) <- -rowSums(dQ)
      rownames(dQ) <- colnames(dQ)
    }

    ## assemble object to return
    obj <- list(independent.Q=iQ,
                dependent.Q=dQ,
                independent.lik.anc = fit.iQ$lik.anc,
                dependent.lik.anc = fit.dQ$lik.anc,
                independent.logL=logLik(fit.iQ),
                dependent.logL=logLik(fit.dQ),
                lik.ratio=2*(logLik(fit.dQ)-logLik(fit.iQ)),
                P=pchisq(2*(logLik(fit.dQ)-logLik(fit.iQ)),
                         df=length(levels(x))+length(levels(y)),
                         lower.tail=FALSE))
    class(obj) <- "fitPagel"
  }else{
    ## assemble object to return
    obj <- list(independent.Q=iQ,
                dependent.Q=NULL, # dQ,
                independent.lik.anc = fit.iQ$lik.anc,
                dependent.lik.anc = NULL, # fit.dQ$lik.anc,
                independent.logL=logLik(fit.iQ),
                dependent.logL=NULL, # logLik(fit.dQ),
                lik.ratio=NULL,
                #               2*(logLik(fit.dQ)-logLik(fit.iQ)),
                #               P=pchisq(2*(logLik(fit.dQ)-logLik(fit.iQ)),
                #                        df=length(levels(x))+length(levels(y)),
                #                        lower.tail=FALSE)
                P = 1
                )
  }

  obj

} # end fitPagel


###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' @import phytools geiger ape phangorn
#' @export

########################################################################

print.fitPagel <- function(x,...){
  cat("\n  Pagel's binary character correlation test:\n")
  cat("\nIndepedent model rate matrix:\n")
  print(x$independent.Q)
  cat("\nDependent model rate matrix:\n")
  print(x$dependent.Q)
  cat("\nIndependent ancestral state likelihoods:\n")
  print(x$independent.lik.anc)
  cat("\nDependent ancestral state likelihoods:\n")
  print(x$dependent.lik.anc)
  cat("\nModel fit:\n")
  obj <- matrix(c(x$independent.logL,x$dependent.logL),2,1)
  rownames(obj) <- c("independent","dependent")
  colnames(obj) <- "log-likelihood"
  print(obj)
  cat("\nHypothesis test result:\n")
  cat(paste("  likelihood-ratio: ",signif(x$lik.ratio,7),"\n"))
  cat(paste("  p-value: ",signif(x$P,7),"\n"))
} # end print.fitPagel

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' @import phytools geiger ape phangorn
#' @export

########################################################################

## fn from phytools github, March 17, 2016.
## https://github.com/liamrevell/phytools/blob/master/R/fitMk.R

fitMk <- function(tree,x,model="SYM",fixedQ=NULL,...){
  if(hasArg(output.liks)) output.liks<-list(...)$output.liks
  else output.liks<-FALSE
  N<-Ntip(tree)
  M<-tree$Nnode
  if(is.matrix(x)){
    x<-x[tree$tip.label,]
    m<-ncol(x)
    states<-colnames(x)
  } else {
    ##############
    ## CC EDITS ##
    ##############

    ## issue: if x variable is a factor with N levels,
    ## and the model (eg. as created within fitPagel) is based on all N levels
    ## BUT some of these levels are empty in x,
    ## fitMk complains that model does not have the right number of columns.

    ## soln: adding empty level to matrix x.

    if(is.factor(x) && any(table(x) == 0)){
      ## identify empty level(s):
      levels.x <- levels(x)
      missing <- which(table(x) == 0)
      ## make matrix of input x:
      x<-to.matrix(x,sort(unique(x)))
      ## make dummy matrix, bind x matrix and empty columns:
      mat <- matrix(NA, nrow=nrow(x), ncol=length(levels.x))
      mat[, c(1:length(levels.x))[-missing]] <- x
      mat[, missing] <- rep(0, nrow(x))
      colnames(mat) <- levels.x
      rownames(mat) <- rownames(x)
      x <- mat
      x<-x[tree$tip.label,]
      m<-ncol(x)
      states<-colnames(x)

    }else{
      ## make matrix of input x:
      x<-to.matrix(x,sort(unique(x)))
    } ############## # end CC edits

    x<-x[tree$tip.label,]
    m<-ncol(x)
    states<-colnames(x)
  }
  if(hasArg(pi)) pi<-list(...)$pi
  else pi<-"equal"
  if(pi[1]=="equal") pi<-setNames(rep(1/m,m),states)
  else if(pi[1]=="estimated"){
    pi<-if(!is.null(fixedQ)) statdist(fixedQ) else statdist(summary(fitMk(tree,x,model),quiet=TRUE)$Q)
    cat("Using pi estimated from the stationary distribution of Q assuming a flat prior.\npi =\n")
    print(round(pi,6))
    cat("\n")
  }
  else pi<-pi/sum(pi)
  if(is.null(fixedQ)){
    if(is.character(model)){
      rate<-matrix(NA,m,m)
      if(model=="ER"){
        k<-rate[]<-1
        diag(rate)<-NA
      } else if(model=="ARD"){
        k<-m*(m-1)
        rate[col(rate)!=row(rate)]<-1:k
      } else if(model=="SYM"){
        k<-m*(m-1)/2
        ii<-col(rate)<row(rate)
        rate[ii]<-1:k
        rate<-t(rate)
        rate[ii]<-1:k
      }
    } else {
      if(ncol(model)!=nrow(model))
        stop("model is not a square matrix")
      if(ncol(model)!=ncol(x))
        stop("model does not have the right number of columns")
      rate<-model
      k<-max(rate)
    }
    Q<-matrix(0,m,m)
  } else {
    rate<-matrix(NA,m,m)
    k<-m*(m-1)
    rate[col(rate)!=row(rate)]<-1:k
    Q<-fixedQ
  }
  index.matrix<-rate
  tmp<-cbind(1:m,1:m)
  rate[tmp]<-0
  rate[rate==0]<-k+1
  liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
  pw<-reorder(tree,"pruningwise")
  lik<-function(pp,output.liks=FALSE,pi){
    if(any(is.nan(pp))||any(is.infinite(pp))) return(1e50)
    comp<-vector(length=N+M,mode="numeric")
    Q[]<-c(pp,0)[rate]
    diag(Q)<--rowSums(Q)
    parents<-unique(pw$edge[,1])
    root<-min(parents)
    for(i in 1:length(parents)){
      anc<-parents[i]
      ii<-which(pw$edge[,1]==parents[i])
      desc<-pw$edge[ii,2]
      el<-pw$edge.length[ii]
      v<-vector(length=length(desc),mode="list")
      for(j in 1:length(v))
        v[[j]]<-matexpo(Q*el[j])%*%liks[desc[j],]
      vv<-if(anc==root) Reduce('*',v)[,1]*pi else Reduce('*',v)[,1]
      comp[anc]<-sum(vv)
      liks[anc,]<-vv/comp[anc]
    }
    if(output.liks)return(liks[1:M+N,,drop=FALSE])
    logL<--sum(log(comp[1:M+N]))
    return(if(is.na(logL)) Inf else logL)
  }
  if(is.null(fixedQ)){
    fit<-nlminb(rep(0.1,k),function(p) lik(p,pi=pi),lower=rep(0,k),upper=rep(1e50,k))
    obj<-list(logLik=-fit$objective,
              rates=fit$par,
              index.matrix=index.matrix,
              states=states,
              pi=pi)
    if(output.liks) obj$lik.anc<-lik(obj$rates,TRUE,pi=pi)
  } else {
    fit<-lik(Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],pi=pi)
    obj<-list(logLik=-fit,
              rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
              index.matrix=index.matrix,
              states=states,
              pi=pi)
    if(output.liks) obj$lik.anc<-lik(obj$rates,TRUE,pi=pi)
  }
  class(obj)<-"fitMk"
  return(obj)
} # end fitMk
