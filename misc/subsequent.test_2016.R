
####################################################
## subsequent.test ## ORIGINAL (w integral score) ##
####################################################

## SEE BELOW FOR NEW/CURRENT VERSION(s) !!!!!!!!!!!!!!
#####################
#
# subsequent.test <- function(snps.reconstruction,
#                             phen.reconstruction,
#                             tree){
#
#   snps.rec <- snps.reconstruction
#   phen.rec <- phen.reconstruction
#
#   ## get tree edges:
#   edges <- tree$edge
#
#   #########################
#   ## Get UNIQUE snps.rec ##
#   #########################
#   temp <- get.unique.matrix(snps.rec, MARGIN=2)
#   snps.rec.unique <- temp$unique.data
#   index <- temp$index
#
#   if(ncol(snps.rec.unique) == ncol(snps.rec)){
#     all.unique <- TRUE
#   }else{
#     all.unique <- FALSE
#   }
#
#   ## work w only unique snps:
#   snps.rec.ori <- snps.rec
#   snps.rec <- snps.rec.unique
#
#   ###############################
#   ## GET SCORE ACROSS BRANCHES ##
#   ###############################
#
#   score <- list()
#   score.p <- list()
#
#   ##############
#   ## FOR LOOP ##
#   ##############
#   for(i in 1:ncol(snps.rec)){
#
#     score3 <- list()
#
#     # ## TEMP -- try w new score3 for all nodes (ie no edges/lengths):
#     # for(e in 1:length(phen.rec)){
#     #
#     #   pi <- phen.rec[e]
#     #   si <- snps.rec[e, i]
#     #
#     #   score3[[e]] <- get.score3.1(Pi=pi, Si=si)
#     # }
#
#     ##############
#     ## FOR LOOP ##
#     ##############
#     for(e in 1:nrow(edges)){
#
#       pa <- phen.rec[edges[e,1]]
#       pd <- phen.rec[edges[e,2]]
#       sa <- snps.rec[edges[e,1], i]
#       sd <- snps.rec[edges[e,2], i]
#       bl <- tree$edge.length[e]
#
#       score3[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)
#
#     } # end (e) for loop
#
#     score3 <- as.vector(unlist(score3))
#
#     score[[i]] <- abs(sum(score3))
#
#     ## TEMP -- get positive-/negative-only scores?
#     score.p[[i]] <- max(abs(sum(score3[score3 > 0])), abs(sum(score3[score3 < 0])))
#
#
#   } # end (i) for loop
#
#   score <- as.vector(unlist(score))
#   names(score) <- colnames(snps.rec)
#
#   ## TEMP -- positive/neg scores:
#   score.p <- as.vector(unlist(score.p))
#   names(score.p) <- colnames(snps.rec)
#
#   # ## check out snps.assoc vs snps:
#   # snps.assoc.ori.ori <- snps.assoc
#   #
#   # snps.ass <- snps.index[snps.assoc]
#   # print("SCORE"); print(score[snps.ass])
#   # print("SCORE P"); print(score.p[snps.ass])
#
#   ################################################
#   ## get values for duplicate snps.rec columns: ##
#   ################################################
#
#   ## get reconstruction for all original sites
#   if(all.unique == TRUE){
#     score.complete <- score
#   }else{
#     score.complete <- score[index]
#     names(score.complete) <- colnames(snps.rec.ori)
#   }
#
#   score <- score.complete
#
#   return(score)
#
# } # end subsequent.test







#'
#' ################
#' ## get.score3 ##
#' ################
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
#' #' @param Pa A numeric value containing either the state,
#' #' or the probability of the state, of the phenotype at a given \emph{ancestral} node.
#' #' @param Pd A numeric value containing either the state,
#' #' or the probability of the state, of the phenotype at a given \emph{descendant} node.
#' #' @param Sa A numeric value containing either the state,
#' #' or the probability of the state, of SNPi at a given \emph{ancestral} node.
#' #' @param Sd A numeric value containing either the state,
#' #' or the probability of the state, of SNPi at a given \emph{descendant} node.
#' #' @param l A numeric value specifying the length of the branch in the phylogenetic tree
#' #' that joins the ancestral and descendant node.
#' #'
#' #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' #' @export
#' #' @examples
#' #'
#' #' ## basic use of fn
#' #' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#' #'
#'
#' ########################################################################
#'
# get.score3 <- function(Pa, Pd, Sa, Sd, l){
#
#   score3 <- NULL
#
#   ## CHECKS:
#   if(length(Pa) > 1) stop("Pa must be of length one
#                                     (i.e., (the probability of)
#                           the phenotypic state of ONE ancestral node.")
#   if(length(Pd) > 1) stop("Pd must be of length one
#                                     (i.e., (the probability of)
#                           the phenotypic state of ONE descendant node.")
#   if(length(Sa) > 1) stop("Sa must be of length one
#                                     (i.e., (the probability of)
#                           the SNPi state of ONE ancestral node.")
#   if(length(Sd) > 1) stop("Sd must be of length one
#                                     (i.e., (the probability of)
#                           the SNPi state of ONE descendant node.")
#
#   ## Original integral-based score:
#   score3 <- l*(((4/3)*Pa*Sa) +
#                  ((2/3)*Pa*Sd) +
#                  ((2/3)*Pd*Sa) +
#                  ((4/3)*Pd*Sd) -
#                  Pa -
#                  Pd -
#                  Sa -
#                  Sd +
#                  1)
#
#
#   ## Simple score1-like score (across edges):
# #   score3 <- (l/2)*(
# #     Pa*Sa + (1 - Pa)*(1 - Sa) -
# #       Pa*(1 - Sa) - (1 - Pa)*Sa +
# #       Pd*Sd + (1 - Pd)*(1 - Sd) -
# #       Pd*(1 - Sd) - (1 - Pd)*Sd)
#
#
#   return(score3)
#
# } # end get.score3



##################
## get.score3.1 ##
##################
## TEMPORARY FN-- TESTING SCORE3 = SCORE1 FOR ALL NODES:
# get.score3.1 <- function(Pi, Si){
#
#   score3 <- NULL
#
#   ## CHECKS:
#   if(length(Pi) > 1) stop("Pi must be of length one
#                                     (i.e., (the probability of)
#                           the phenotypic state of ONE ancestral node.")
#
#   if(length(Si) > 1) stop("Si must be of length one
#                                     (i.e., (the probability of)
#                           the SNPi state of ONE ancestral node.")
#
#
#   ## Simple score1-like score (for all/individual nodes):
#   score3 <- Pi*Si + (1 - Pi)*(1 - Si) -
#       Pi*(1 - Si) - (1 - Pi)*Si
#
#
#   return(score3)
#
# } # end get.score3.1



####################################################################################################################################

#############################
## NEW Pagel-like SCORE3 ? ##
#############################

## INCLUDE:

## Pr(state11 | state11) # maintained
## Pr(state01 | state01) # simultaneous
## Pr(state01 | state11) # subsequent

## And the REVERSE x2:
## ie.
## 01 & 10
## AND SNP|Phen & Phen|SNP

# snps <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_snps.Rdata"))
# snps.rec <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_res.Rdata"))$dat$snps.rec
# snps.assoc <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_performance.Rdata"))$snps.assoc
# phen <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_phen.Rdata"))
# phen.rec <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_res.Rdata"))$dat$phen.rec
# tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_tree.Rdata"))

# snps <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_30_snps.Rdata"))
# snps.rec <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_30_res.Rdata"))$dat$snps.rec
# snps.assoc <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_30_performance.Rdata"))$snps.assoc
# phen <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_30_phen.Rdata"))
# phen.rec <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_30_res.Rdata"))$dat$phen.rec
# tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_30_tree.Rdata"))


###################################################################################
## NEW (temporary/replacement) FUNCTION for subsequent test more like PAGEL test ##
###################################################################################

#####################
## subsequent.test ##
#####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param snps.reconstruction A matrix containing the terminal and reconstructed
#' ancestral states of SNPs for all nodes in the tree.
#' @param phen.reconstruction A vector containing the terminal and reconstructed
#' ancestral states of the phenotype for all nodes in the tree.
#' @param tree A phylo object containing the tree representing the ancestral relationships
#' between the individuals for which snps and phen are known.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export


########################################################################

subsequent.test <- function(snps.reconstruction,
                            phen.reconstruction,
                            tree){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction

  ## get tree edges:
  edges <- tree$edge

  #########################
  ## Get UNIQUE snps.rec ##
  #########################
  temp <- get.unique.matrix(snps.rec, MARGIN=2)
  snps.rec.unique <- temp$unique.data
  index <- temp$index

  if(ncol(snps.rec.unique) == ncol(snps.rec)){
    all.unique <- TRUE
  }else{
    all.unique <- FALSE
  }

  ## work w only unique snps:
  snps.rec.ori <- snps.rec
  snps.rec <- snps.rec.unique

###############################
## GET SCORE ACROSS BRANCHES ##
###############################

score <- SCORE3.1 <- SCORE3.2 <- SCORE3.3 <- SCORE3.4 <- SCORE3.5 <- SCORE3.raw <- SCORE3.p <- SCORE3.n <- SCORE3.RAW <- list()
SCORE3.integral.sum <- SCORE3.integral.max <- SCORE3.1.edges.sum <- SCORE3.1.edges.max <- list()
Q.ew <- Q.uw <- list()

## HANDLE PHEN.REC:
## get rid of any 0.5 values (make NA):
phen.rec.ori <- phen.rec
phen.rec <- replace(phen.rec, which(! phen.rec %in% c("A", "B", 0, 1)), NA)
## make numeric 0/1:
phen.rec <- as.numeric(as.factor(phen.rec)) - 1

##############
## FOR LOOP ##
##############
for(i in 1:ncol(snps.rec)){

  score3 <- list()
  score3.integral <- score3.1.edges <- list()

  #####################
  ## get ~ Q matrix: ##
  #####################
  ## HANDLE SNPS.REC:
  snp.rec <- snps.rec[,i]
  snp.rec <- replace(snp.rec, which(! snp.rec %in% c(0, 1)), NA)

  ## For each EDGE in EDGES, PASTE ANC SNP.PHEN w DEC SNP.PHEN as SaSd|PaPd:
  ##############
  ## FOR LOOP ##
  ##############
  edges <- tree$edge
  for(e in 1:nrow(edges)){

    Pa <- phen.rec[edges[e,1]]
    Pd <- phen.rec[edges[e,2]]
    Sa <- snp.rec[edges[e,1]]
    Sd <- snp.rec[edges[e,2]]
    l <- tree$edge.length[e]

    if(!any(is.na(Sa), is.na(Sd), is.na(Pa), is.na(Pd))){
      score3[[e]] <- paste(Sa, Pa, "|", Sd, Pd, sep="")
    }else{
      score3[[e]] <- NA
    }

    # score3[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)

    ## Original integral-based score:
    score3.integral[[e]] <- l*(((4/3)*Pa*Sa) +
                              ((2/3)*Pa*Sd) +
                              ((2/3)*Pd*Sa) +
                              ((4/3)*Pd*Sd) -
                              Pa -
                              Pd -
                              Sa -
                              Sd +
                              1)


    ## Simple score1-like score (across edges):
    score3.1.edges[[e]] <- (l/2)*(
                            Pa*Sa + (1 - Pa)*(1 - Sa) -
                              Pa*(1 - Sa) - (1 - Pa)*Sa +
                              Pd*Sd + (1 - Pd)*(1 - Sd) -
                              Pd*(1 - Sd) - (1 - Pd)*Sd)

  } # end (e) for loop

  ## Combine all edges into one FACTOR for this snps.rec column:
  score3 <- as.vector(unlist(score3))
  names(score3) <- c(1:length(score3))

  ## Combine alternative scores into one per snp: #############################

  ## integral-based score
  score3.integral <- as.vector(unlist(score3.integral))
  ## save BOTH the sum and the max-only score:
  ## sum:
  SCORE3.integral.sum[[i]] <- abs(sum(score3.integral))
  ## max -- positive-/negative-only scores?
  SCORE3.integral.max[[i]] <- max(abs(sum(score3.integral[score3.integral > 0])), abs(sum(score3.integral[score3.integral < 0])))

  ## score1-like cross-edges score:
  score3.1.edges <- as.vector(unlist(score3.1.edges))
  ## save BOTH the sum and the max-only score:
  ## sum:
  SCORE3.1.edges.sum[[i]] <- abs(sum(score3.1.edges))
  ## max -- positive-/negative-only scores?
  SCORE3.1.edges.max[[i]] <- max(abs(sum(score3.1.edges[score3.1.edges > 0])), abs(sum(score3.1.edges[score3.1.edges < 0])))

  ###########################################################################

  ## Use relative freq of each factor level as prob (or incorporate edge length to get ~ rate)...
  # table(score3)

  ## get separate positive and negative scores for concordant (00,11) and discordant (01,10) associations:
  score3.p.noms <- c("00|00", "11|11", "00|11", "11|00", "01|00", "10|00", "01|11", "10|11")
  score3.n.noms <- c("00|01", "00|10", "11|01", "11|10", "01|01", "10|10", "01|10", "10|01")
  score3.p <- replace(score3, which(!score3 %in% score3.p.noms), NA)
  score3.n <- replace(score3, which(!score3 %in% score3.n.noms), NA)

  tab.p <- table(score3.p)
  tab.n <- table(score3.n)

  SCORE3.raw[[i]] <- score3
  SCORE3.p[[i]] <- score3.p
  SCORE3.n[[i]] <- score3.n

  ################################################################################################

  ##################################
  ## Get UN-edge-weighted SCORES: ##
  ##################################

  # i <- 1369

  ##################
  ## SIMPLE score ## (all equal +/- 1): # p = 1.76x n
  ##################
  score3.p.u <- sum(tab.p) # 125
  score3.n.u <- sum(tab.n) # 71

  ## get unewighted simple score as max of pos and neg scores:
  SCORE3.1[[i]] <- score3.u <- max(score3.p.u, score3.n.u) # 125

  ############################
  ## COMPLEX weighted score ## (simultaneous > maintained > subsequent (?)): # p = 2.1x n
  ############################
  sim.p <- tab.p[which(names(tab.p) %in% c("00|11", "11|00"))]
  maint.p <- tab.p[which(names(tab.p) %in% c("00|00", "11|11"))]
  subsq.p <- tab.p[which(names(tab.p) %in% c("01|00", "01|11", "10|00", "10|11"))]
  ## positive complex score:
  score3.p.w.u <- sum(sim.p*3, maint.p*2, subsq.p*1) # 202

  sim.n <- tab.n[which(names(tab.n) %in% c("01|10", "10|01"))]
  maint.n <- tab.n[which(names(tab.n) %in% c("01|01", "10|10"))]
  subsq.n <- tab.n[which(names(tab.n) %in% c("00|01", "00|10", "11|01", "11|10"))]
  ## negative complex score:
  score3.n.w.u <- sum(sim.n*3, maint.n*2, subsq.n*1) # 96

  ## get weighted (un-edge-weighted) score as mx of pos and neg scores:
  SCORE3.2[[i]] <- score3.w.u <- max(score3.p.w.u, score3.n.w.u) # 202


  ###########################################
  ## Get SCORES weighted by EDGE LENGTH??: ##   ###   ###   ###   ###   ###   ###   ###   ###
  ###########################################

  ##################
  ## SIMPLE score ## (all equal +/- 1): # p = 2.06x n
  ##################
  score3.p.ew <- sum(tree$edge.length[which(score3 %in% score3.p.noms)]) # 6.7
  score3.n.ew <- sum(tree$edge.length[which(score3 %in% score3.n.noms)]) # 3.25

  ## Keep only higher of the two scores:
  SCORE3.3[[i]] <- score3.ew <- max(score3.p.ew, score3.n.ew) # 6.7

  ############################
  ## COMPLEX weighted score ## (simultaneous > maintained > subsequent (?)): # p = 2.5x n
  ############################
  sim.p.ew <- tree$edge.length[which(score3 %in% c("00|11", "11|00"))]
  maint.p.ew <- tree$edge.length[which(score3 %in% c("00|00", "11|11"))]
  subsq.p.ew <- tree$edge.length[which(score3 %in% c("01|00", "01|11", "10|00", "10|11"))]
  ## positive complex score:
  score3.p.w.ew <- sum(sim.p.ew*3, maint.p.ew*2, subsq.p.ew*1) # 11.03

  sim.n.ew <- tree$edge.length[which(score3 %in% c("01|10", "10|01"))]
  maint.n.ew <- tree$edge.length[which(score3 %in% c("01|01", "10|10"))]
  subsq.n.ew <- tree$edge.length[which(score3 %in% c("00|01", "00|10", "11|01", "11|10"))]
  ## negative complex score:
  score3.n.w.ew <- sum(sim.n.ew*3, maint.n.ew*2, subsq.n.ew*1) # 4.37

  ##############

  # Q.noms <- paste(c(0,0,1,1), c(0,1,0,1), sep="|")
  # ## (SNPphen.anc|SNPphen.dec)
  # m <- c("00|00", "00|01", "00|10", "00|11",
  #        "01|00", "01|01", "01|10", "01|11",
  #        "10|00", "10|01", "10|10", "10|11",
  #        "11|00", "11|01", "11|10", "11|11")
  # ## Get the sum of the edge lengths over which that type of association happens,
  # ## divided by the sum of edge lengths for which we have a score at this SNP
  # ## (ie. the edges for which neither the ancestor nor the descendant has
  # ## an unknown reconstructed SNP/phen (ie. value of 0.5))
  # m.corr <- c(2, 0.75, 0.75, 1, 3, 0.5, 0.25,  3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2)
  # mat.corr <- matrix(m.corr, ncol=4, nrow=4, byrow=TRUE)
  # ## sum to 1:
  # mat.corr <- mat.corr/sum(mat.corr)
  # ## (SNP|phen) -- anc in rows, dec in cols:
  # rownames(mat.corr) <- colnames(mat.corr) <- c("0|0", "0|1", "1|0", "1|1")
  # Q <- mat.corr

  ##############
  ## get rate matrix:
  # s3 <- SCORE3$corr.dat[[6]]
  #
  # for(s in 1:length(snps.assoc)){
  # score3 <- s3[[snps.assoc[s]]]

  Q.noms <- paste(c(0,0,1,1), c(0,1,0,1), sep="|")
  ## (SNPphen.anc|SNPphen.dec)
  m <- c("00|00", "00|01", "00|10", "00|11",
         "01|00", "01|01", "01|10", "01|11",
         "10|00", "10|01", "10|10", "10|11",
         "11|00", "11|01", "11|10", "11|11")
  ## Get the sum of the edge lengths over which that type of association happens,
  ## divided by the sum of edge lengths for which we have a score at this SNP
  ## (ie. the edges for which neither the ancestor nor the descendant has
  ## an unknown reconstructed SNP/phen (ie. value of 0.5))
  m.ew <- sapply(c(1:length(m)),
                 function(e)
                   sum(tree$edge.length[which(score3 == m[e])])) / sum(tree$edge.length[!is.na(score3)])
  mat.ew <- matrix(m.ew, ncol=4, nrow=4, byrow=TRUE)
  ## (SNP|phen) -- anc in rows, dec in cols:
  rownames(mat.ew) <- colnames(mat.ew) <- c("0|0", "0|1", "1|0", "1|1")
  # print(mat.ew)

  ## OR ##
  ## Get unweighed "probs"??
  m.uw <- sapply(c(1:length(m)),
                 function(e)
                   length(which(score3 == m[e]))) / length(score3[!is.na(score3)])
  mat.uw <- matrix(m.uw, ncol=4, nrow=4, byrow=TRUE)

  ## (SNP|phen) -- anc in rows, dec in cols:
  rownames(mat.uw) <- colnames(mat.uw) <- c("0|0", "0|1", "1|0", "1|1")
  # print(mat.uw)

  # sum(mat[,c(1,4)])
  # sum(mat[,c(2,3)])

  # print(sum(mat.ew[,c(1,4)]))
  # print(sum(mat.ew[,c(2,3)]))
  #
  # print(sum(mat.uw[,c(1,4)]))
  # print(sum(mat.uw[,c(2,3)]))
  # }

  Q.ew[[i]] <- mat.ew
  Q.uw[[i]] <- mat.uw

  # load("/home/caitiecollins/treeWAS/misc/SCORE3.set1_31.Rdata")
  # snps.assoc <- out$performance[[1]]$snps.assoc
  # tree <- out$tree[[1]]
  # snps.rec <- out$res$set1_31$dat$snps.rec
  # phen.rec <- out$res$set1_31$dat$phen.rec

  ##
  # var.rec <- snps.rec[,snps.assoc[1]]
  # var.rec <- round(var.rec)
  # var.rec <- replace(var.rec, which(var.rec == 0), "A")
  # var.rec <- replace(var.rec, which(var.rec == 1), "B")
  # if(any(!var.rec %in% c("A", "B"))) var.rec <- replace(var.rec, which(!var.rec %in% c("A", "B")), NA)
  # plot.phen(tree, phen.nodes=var.rec)
  # title("SNP 304", line=-1)


  # Q.ew <- SCORE3$corr.dat[[9]]
  # Q.uw <- SCORE3$corr.dat[[10]]
  #
  # Q.ew[[snps.assoc[1]]]
  # Q.uw[[snps.assoc[1]]]
  #
  # var <- snps.rec[,snps.assoc[1]][-which(!snps.rec[,snps.assoc[1]] %in% c(0,1))]
  # table(var)/length(var)
  # sum(Q.uw[[snps.assoc[1]]][, c(1,2)])
  # sum(Q.uw[[snps.assoc[1]]][, c(3,4)])

  ## sum.mat <- Reduce("+", list(mat, mat, mat)) # add matrices togeter by cell.
  ############

  ## KEEP ONLY THE HIGHER OF THE TWO SCORES!  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ## NOTE -- ORRRRRRRRRRRRRRR SHOULD WE KEEP THE ABS VAL OF THE DIFFERENCE BETWEEN THE TWO SCORES??????????????????????????    <----- (??)
  SCORE3.4[[i]] <- score3.w.ew <- max(score3.p.w.ew, score3.n.w.ew) # 11.03

  SCORE3.5[[i]] <- abs(score3.p.w.ew - score3.n.w.ew)

  ## for now, using weighted complex as subsequent score:
  ## trying w unweighted simple score (makes more sense, given sim is using probs without factoring in branch lengths..)
  score[[i]] <- SCORE3.1[[i]]
  } # end (i) for loop
  ################################################################################################

  score <- as.vector(unlist(score))
  names(score) <- colnames(snps.rec)

  ## TEMP!: ##
  ############
  ## name SCORE3 AND re-duplicate all elements AND save it to compare to new score....
  SCORE3.1 <- as.vector(unlist(SCORE3.1))
  names(SCORE3.1) <- colnames(snps.rec)

  SCORE3.2 <- as.vector(unlist(SCORE3.2))
  names(SCORE3.2) <- colnames(snps.rec)

  SCORE3.3 <- as.vector(unlist(SCORE3.3))
  names(SCORE3.3) <- colnames(snps.rec)

  SCORE3.4 <- as.vector(unlist(SCORE3.4))
  names(SCORE3.4) <- colnames(snps.rec)

  SCORE3.5 <- as.vector(unlist(SCORE3.5))
  names(SCORE3.5) <- colnames(snps.rec)

  ################################################
  ## get values for duplicate snps.rec columns: ##
  ################################################

  ## get reconstruction for all original sites
  if(all.unique == TRUE){
    score.complete <- score
  }else{
    score.complete <- score[index]
    names(score.complete) <- colnames(snps.rec.ori)

    SCORE3.1 <- SCORE3.1[index]
    names(SCORE3.1) <- colnames(snps.rec.ori)

    SCORE3.2 <- SCORE3.2[index]
    names(SCORE3.2) <- colnames(snps.rec.ori)

    SCORE3.3 <- SCORE3.3[index]
    names(SCORE3.3) <- colnames(snps.rec.ori)

    SCORE3.4 <- SCORE3.4[index]
    names(SCORE3.4) <- colnames(snps.rec.ori)

    SCORE3.5 <- SCORE3.5[index]
    names(SCORE3.5) <- colnames(snps.rec.ori)

    ## alternative scores:
    ## score 1-like cross-edges
    SCORE3.1.edges.sum <- SCORE3.1.edges.sum[index]
    names(SCORE3.1.edges.sum) <- colnames(snps.rec.ori)

    SCORE3.1.edges.max <- SCORE3.1.edges.max[index]
    names(SCORE3.1.edges.max) <- colnames(snps.rec.ori)

    ## integral-based score:
    SCORE3.integral.sum <- SCORE3.integral.sum[index]
    names(SCORE3.integral.sum) <- colnames(snps.rec.ori)

    SCORE3.integral.max <- SCORE3.integral.max[index]
    names(SCORE3.integral.max) <- colnames(snps.rec.ori)
  }


  SCORE3 <- list("SCORE3.1"=SCORE3.1, "SCORE3.2"=SCORE3.2, "SCORE3.3"=SCORE3.3, "SCORE3.4"=SCORE3.4, "SCORE3.5"=SCORE3.5,
                 "SCORE3.raw"=SCORE3.raw, "SCORE3.p"=SCORE3.p, "SCORE3.n"=SCORE3.n,
                 "SCORE3.integral.sum"=SCORE3.integral.sum, "SCORE3.integral.max"=SCORE3.integral.max,
                 "SCORE3.1.edges.sum"=SCORE3.1.edges.sum, "SCORE3.1.edges.max"=SCORE3.1.edges.max,
                 "Q.ew"=Q.ew, "Q.uw"=Q.uw)

  score <- score.complete

  SCORE <- list("score" = score,
                "SCORE3" = SCORE3)

  # return(score)
  return(SCORE)

} # end subsequent.test
##



# ## load saved SCORE3 list
# corr.dat <- SCORE3$corr.dat
# corr.sim <- SCORE3$corr.sim
#
# snps.assoc <- out$performance[[1]]$snps.assoc
#
# hist(corr.sim[[1]])
# cd <- corr.dat[[1]][snps.assoc]
# cd
# sapply(c(1:length(cd)), function(e) length(which(corr.sim[[1]] > cd[e])))
#
# hist(corr.sim[[2]])
# cd <- corr.dat[[2]][snps.assoc]
# cd
# sapply(c(1:length(cd)), function(e) length(which(corr.sim[[2]] > cd[e])))
#
# hist(corr.sim[[3]])
# cd <- corr.dat[[3]][snps.assoc]
# cd
# sapply(c(1:length(cd)), function(e) length(which(corr.sim[[3]] > cd[e])))
#
# hist(corr.sim[[4]])
# cd <- corr.dat[[4]][snps.assoc]
# cd
# sapply(c(1:length(cd)), function(e) length(which(corr.sim[[4]] > cd[e])))
#
# hist(corr.sim[[5]])
# cd <- corr.dat[[5]][snps.assoc]
# cd
# sapply(c(1:length(cd)), function(e) length(which(corr.sim[[5]] > cd[e])))





#
#









#





#
