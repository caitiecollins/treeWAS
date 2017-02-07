


##############
## ace.test ##
##############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param snps A matrix containing binary snps for all individuals.
#' @param phen A factor containing the phenotype for which to test for association.
#' @param tree A phylo object containing the tree representing the ancestral relationships
#' between the individuals for which snps and phen are known.
#' @param test A character string containing one or both (the default) of "simultaneous" and "subsequent",
#' specifying which test(s) of association to perform on the data.
#' @param snps.reconstruction Either a character string containing one of "parsimony" or "ace",
#' specifying whether the ancestral states of \code{snps} should be reconstructed by
#' parsimony (faster; the default for now) or ancestral character estimation (as performed in \emph{ape}, slower),
#' or an object containing a parsimonious reconstruction previously performed by the user.
#' @param phen.reconstruction Either a character string containing one of "parsimony" or "ace",
#' specifying whether the ancestral states of \code{phen} should be reconstructed by
#' parsimony (the default for now) or ancestral character estimation (as performed in \emph{ape}),
#' or an object containing a parsimonious reconstruction previously performed by the user.
#' @param method A character string specifying the type of ACE method to implement (only used if
#' \code{snps.rec} or \code{phen.rec} is set to "ace").
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## example: ##
#'
#' # load data
#' data("snps.ace")
#' data("phen.ace")
#' data("tree.ace")
#' snps <- snps.ace
#' phen <- phen.ace
#' tree <- tree.ace
#'
#' ## run ace.test on subset of data (otherwise SLOW!)
#' out <- ace.test(snps[,1:10], phen, tree, method="discrete")
#'
#' ## examine output
#' round(out, 4)
#'
#' @import phangorn

########################################################################

## Eg:
#
# data("snps.ace")
# data("phen.ace")
# data("tree.ace")
#
# snps <- snps.ori <- snps.ace
# phen <- phen.ori <- phen.ace
# tree <- tree.ori <- tree.ace
# test <- c("simultaneous", "subsequent")
# snps.reconstruction <- "parsimony"
# phen.reconstruction <- "parsimony"
# method <- "discrete"

##############################################################################################################################################################

###########################
## OLD ace.test FUNCTION ##
###########################
ace.test <- function(snps, phen, tree,
                     test = c("simultaneous", "subsequent"),
                     snps.reconstruction = "parsimony",
                     phen.reconstruction = "parsimony",
                     method = "discrete"){


  ##############
  ## OPTIONS: ##
  ##############

  ## TEST: ##
  ###########
  ## (1) Simultaneous Test for simultaneous substitutions,
  ## allowing for complementary pathways.
  ## (2) Subsequent Test for the effect of one variable's state on the maintenance or
  ## substitution of the other variable's state.

  ## Allow partial matching of argument names:
  test <- match.arg(arg = test,
                    choices =  c("simultaneous", "subsequent"),
                    several.ok = TRUE)
  snps.reconstruction <- match.arg(snps.reconstruction)
  phen.reconstruction <- match.arg(phen.reconstruction)

  ## POSSIBLE SUB-OPTIONS: ##
  ###########################
  ## (1) Parsimony on SNPs and/or Phen
  ## (2) ACE on SNPs and/or Phen
  ## (3) User data on SNPs and/or Phen
  ## (?????????????? ARE USERS LIKELY TO WANT TO PROVIDE ONLY PARSIMONIOUS RECONSTRUCTIONS OR POSSIBLY ACE AS WELL???? AND HOW TO SPECIFY IF SO ???????????????)
  ## NOTE: Regardless of combination of sub-options selected, these will be required to be the same choices for all association tests run.

  edges <- tree$edge

  ###############################
  ## get unique SNPs patterns: ##
  ###############################
  temp <- get.unique.matrix(snps, MARGIN=2)
  snps.unique <- temp$unique.data
  index <- temp$index
  if(ncol(snps.unique) == ncol(snps)){
    all.unique <- TRUE
  }else{
    all.unique <- FALSE
  }

  ## work w only unique snps:
  snps.ori <- snps
  snps <- snps.unique

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  ########################################
  ## RUN PARSIMONY on SNPS and/or PHEN: ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ########################################

  ## SNPS parsimony ##
  if(snps.reconstruction == "parsimony"){
    ## run get.ancestral.pars
    snps.pars <- get.ancestral.pars(var=snps, tree=tree)

    ## get elements of output
    snps.rec <- snps.pars$var.rec
    snps.subs.edges <- snps.pars$subs.edges
  }

  ## PHEN parsimony ##
  if(phen.reconstruction == "parsimony"){
    ## run get.ancestral.pars
    phen.pars <- get.ancestral.pars(var=phen, tree=tree)

    ## get elements of output
    phen.rec <- phen.pars$var.rec
    phen.subs.edges <- phen.pars$subs.edges
  }

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


  ##################################
  ## RUN ACE on SNPS and/or PHEN: ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ##################################

  ## SNPS ace ##

  if(snps.reconstruction == "ace"){
    ##################################
    ## get diffs for (unique) SNPS: ##
    ##################################
    snps.diffs <- snps.ACE <- list()
    # system.time == 105.042 elapsed
    # w unique ncol(snps) == 2500
    for(i in 1:ncol(snps)){
      ## get variable
      var <- snps[,i]
      var.terminal <- var
      snps.ACE[[i]] <- ace(var, tree, type=method)
      var.internal <- snps.ACE[[i]]$lik.anc[,2]
      var <- c(var.terminal, var.internal)

      ## get differences for this variable's ace likelihoods
      snps.diffs[[i]] <- get.ace.diffs(var, edges)
    }
  }



  ## PHEN ace ##

  if(phen.reconstruction == "ace"){
    #########################
    ## get diffs for PHEN: ##
    #########################
    ## do we need to check phen is numeric??
    phen.terminal <- phen
    phen.ACE <- ace(phen, tree, type=method)
    phen.internal <- phen.ACE$lik.anc[,2]
    var <- c(phen.terminal, phen.internal)

    phen.diffs <-  get.ace.diffs(var, edges)
  }


  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  #############################
  ## RUN ASSOCIATION TEST(s) ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  #############################

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


  #######################
  ## SIMULTANEOUS TEST ##
  #######################
  if("simultaneous" %in% test){

    ## OPTIONS: ##

    ## (1) #####################
    ## SNPS & Phen PARSIMMONY ##
    ############################
    ## Test for %branches w subs (in a given direction) shared btw SNPs & Phen.

    ## Similar to Corr. Score, but only for snp.subs,
    ## w no penalisation for phen.subs unaccompanied by snp.subs
    ## (ie. allowing for complementary pathways).
    ## So, no denominator term...

    if(snps.reconstruction == "parsimony" & phen.reconstruction == "parsimony"){
      ace.score <- list()
      phen.pos <- phen.subs.edges$pos
      phen.neg <- phen.subs.edges$neg
      phen.total <- phen.subs.edges$total
      for(i in 1:length(snps.subs.edges)){
        snp.pos <- snps.subs.edges[[i]]$pos
        snp.neg <- snps.subs.edges[[i]]$neg

        ace.score[[i]] <- abs(length(c(which(snp.pos %in% phen.pos),
                                       which(snp.neg %in% phen.neg)))
                              -
                                length(c(which(snp.pos %in% phen.neg),
                                         which(snp.neg %in% phen.pos)))) ##/ length(phen.total)

      } # end for loop
      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score2.1 <- ace.score
    }

    ## (2) ##################################
    ## SNPs/Phen ACE & Phen/SNPs PARSIMONY ##
    #########################################
    ## Test for relative diffs btw anc. & dec. $lik.anc (in ACE var) on branches w a substitution
    ## (in PARSIMONY var).

    ## (2.a) SNPs PARSIMONY & Phen ACE ##
    #####################################
    if(snps.reconstruction == "parsimony" & phen.reconstruction == "ace"){
      ace.score <- list()
      for(i in 1:length(snps.subs.edges)){
        snps.pos <- snps.subs.edges[[i]]$pos
        snps.neg <- snps.subs.edges[[i]]$neg
        snps.total <- snps.subs.edges[[i]]$total
        ## get the absolute value of
        ## ((the sum of all the positive (0-->1) snps subs' phen.diffs)
        ## MINUS
        ## (the sum of all the negative (1-->0) snps subs' phen.diffs))
        ace.score[[i]] <- abs(sum(phen.diffs[snps.pos])
                              -
                                sum(phen.diffs[snps.neg])) ##/ length(snps.total)
      } # end for loop
      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score2.2a <- ace.score
    }

    ## (2.b) SNPs ACE & Phen PARSIMONY ##
    #####################################
    if(snps.reconstruction == "ace" & phen.reconstruction == "parsimony"){
      ace.score <- list()
      phen.pos <- phen.subs.edges$pos
      phen.neg <- phen.subs.edges$neg
      phen.total <- phen.subs.edges$total
      for(i in 1:length(snps.diffs)){
        ## get the absolute value of
        ## ((the sum of all the positive (0-->1) phen subs' snps.diffs)
        ## MINUS
        ## (the sum of all the negative (1-->0) phen subs' snps.diffs))
        ace.score[[i]] <- abs(sum(snps.diffs[[i]][phen.pos])
                              -
                                sum(snps.diffs[[i]][phen.neg])) ##/ length(phen.total)
      } # end for loop
      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score2.2b <- ace.score
    }


    ## (3) ##############
    ## SNPs & Phen ACE ##
    #####################
    ## Multiplicative test: ((ACE.diffs_SNP) * (ACE.diffs_phen))
    if(snps.reconstruction == "ace" & phen.reconstruction == "ace"){
      ace.score <- list()
      for(i in 1:length(snps.diffs)){
        snp.diffs <- snps.diffs[[i]]
        sp.diffs <- snp.diffs * phen.diffs

        ace.score[[i]] <- abs(sum(sp.diffs)) ##
        ## NOTE: unlike the Simultaneous Test variants 1, 2a, 2b, this ace*ace version (3) does NOT
        ## have an obvious natural denominator, so it cannot be divided to be forced btw 0 and 1...
        ## In light of this, might it be better to be consistent and not divide the other 3 variants above..?
      }
      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score2.3 <- ace.score
    }

    ## temp: Check results for EG dataset:
    #   sort(ace.score1, decreasing=TRUE)[1:10]
    #   sort(ace.score2a, decreasing=TRUE)[1:10]
    #   # which(names(sort(ace.score2a, decreasing=TRUE)) == "3") # "3" == 25th
    #   sort(ace.score2b, decreasing=TRUE)[1:10]
    #   sort(ace.score3, decreasing=TRUE)[1:10]

    ############################################
    ## get values for duplicate snps columns: ##
    ############################################
    if(all.unique == TRUE){
      ace.score.complete <- ace.score
    }else{
      ace.score.complete <- rep(NA, ncol(snps.ori))
      for(i in 1:ncol(snps.unique)){
        ace.score.complete[which(index == i)] <- ace.score[i]
      }
      names(ace.score.complete) <- colnames(snps.ori)
    }

  } # end "simultaneous" %in% test
  ## END SIMULTANEOUS TEST ##

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  #########################
  ## SUBSEQUENT TEST ##
  #########################
  if("subsequent" %in% test){

    ## OPTIONS: ##

    ## (1) #####################
    ## SNPS & Phen PARSIMMONY ##
    ############################
    ## ...

    if(snps.reconstruction == "parsimony" & phen.reconstruction == "parsimony"){

      ##############
      ## FOR LOOP ##
      ##############
      ace.score <- list()

      for(i in 1:ncol(snps.rec)){

        score <- list()

        ##############
        ## FOR LOOP ##
        ##############
        for(e in 1:nrow(edges)){

          pa <- phen.rec[edges[e,1]]
          pd <- phen.rec[edges[e,2]]
          sa <- snps.rec[edges[e,1], i]
          sd <- snps.rec[edges[e,2], i]
          bl <- tree$edge.length[e]

          score[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)

        } # end (e) for loop

        ace.score[[i]] <- abs(sum(as.vector(unlist(score))))

      } # end (i) for loop

      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score3.1 <- ace.score

    } # end SNPs & phen by PARSIMONY


    ## (2) ##################################
    ## SNPs/Phen ACE & Phen/SNPs PARSIMONY ##
    #########################################

    ## (2.a) SNPs PARSIMONY & Phen ACE ##
    #####################################
    if(snps.reconstruction == "parsimony" & phen.reconstruction == "ace"){

      #######################
      ## Handle phen probs ##
      #######################
      ## CAREFUL--NOTE: names of phen (and snp row.names) are NOT in order 1:100
      ## (SHOULD THEY BE RE-ORDERED???) !!!!!!!!!!!!!!!!!!!!!!!!!!?????
      ## CAREFUL 2--NOTE: code below ASSUMES phen is numeric and binary!
      ## (Should add check for this above if not already there!!!!)

      ## get phen "lik.anc" (ie. states) for terminal nodes:
      phen.temp <- phen
      levels.phen <- unique(phen.temp)
      if(length(levels.phen) != 2) stop("At present, phenotype must be binary with two levels.")
      if(!is.numeric(phen.temp) | !all(levels.phen %in% c(0,1))){
        phen.temp <- as.numeric(as.factor(phen.temp))
        names(phen.temp) <- names(phen)
        phen.temp <- phen.temp - 1
      }
      phen.lik.term <- phen.temp

      ## get $lik.anc for internal nodes:
      phen.lik.int <- phen.ACE$lik.anc[,2]

      ## get phen.lik total:
      phen.lik <- c(phen.lik.term, phen.lik.int)

      ##############
      ## FOR LOOP ##
      ##############
      ace.score <- list()

      for(i in 1:ncol(snps.rec)){

        score <- list()

        ##############
        ## FOR LOOP ##
        ##############
        for(e in 1:nrow(edges)){

          pa <- phen.lik[edges[e,1]]
          pd <- phen.lik[edges[e,2]]
          sa <- snps.rec[edges[e,1], i]
          sd <- snps.rec[edges[e,2], i]
          bl <- tree$edge.length[e]

          score[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)

        } # end (e) for loop

        ace.score[[i]] <- abs(sum(as.vector(unlist(score))))

      } # end (i) for loop

      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score3.2a <- ace.score

    } # end SNPs by PARSIMONY & phen by ACE



    ## (2.b) SNPs ACE & Phen PARSIMONY ##
    #####################################
    if(snps.reconstruction == "ace" & phen.reconstruction == "parsimony"){

      ##############
      ## FOR LOOP ##
      ##############
      ace.score <- list()

      for(i in 1:length(snps.ACE)){

        #######################
        ## Handle SNPS probs ##
        #######################
        ## get $lik.anc for internal nodes:
        snp.lik.int <- snps.ACE[[i]]$lik.anc[,2]

        ## get "lik.anc" (ie. states) for terminal nodes:
        snp.lik.term <- snps[,i]

        ## get snp.lik total:
        snp.lik <- c(snp.lik.term, snp.lik.int)

        score <- list()

        ##############
        ## FOR LOOP ##
        ##############
        for(e in 1:nrow(edges)){

          pa <- phen.rec[edges[e,1]]
          pd <- phen.rec[edges[e,2]]
          sa <- snp.lik[edges[e,1]]
          sd <- snp.lik[edges[e,2]]
          bl <- tree$edge.length[e]

          score[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)

        } # end (e) for loop

        ace.score[[i]] <- abs(sum(as.vector(unlist(score))))

      } # end (i) for loop

      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score3.2b <- ace.score

    } # end SNPs by ACE & phen by PARSIMONY


    ## (3) ##############
    ## SNPs & Phen ACE ##
    #####################
    if(snps.reconstruction == "ace" & phen.reconstruction == "ace"){

      #######################
      ## Handle phen probs ##
      #######################

      ## CAREFUL--NOTE: names of phen (and snp row.names) are NOT in order 1:100
      ## (SHOULD THEY BE RE-ORDERED???) !!!!!!!!!!!!!!!!!!!!!!!!!!?????
      ## CAREFUL 2--NOTE: code below ASSUMES phen is numeric and binary!
      ## (Should add check for this above if not already there!!!!)

      ## get phen "lik.anc" (ie. states) for terminal nodes:
      phen.temp <- phen
      levels.phen <- unique(phen.temp)
      if(length(levels.phen) != 2) stop("At present, phenotype must be binary with two levels.")
      if(!is.numeric(phen.temp) | !all(levels.phen %in% c(0,1))){
        phen.temp <- as.numeric(as.factor(phen.temp))
        names(phen.temp) <- names(phen)
        phen.temp <- phen.temp - 1
      }
      phen.lik.term <- phen.temp

      ## get $lik.anc for internal nodes:
      phen.lik.int <- phen.ACE$lik.anc[,2]

      ## get phen.lik total:
      phen.lik <- c(phen.lik.term, phen.lik.int)


      ##############
      ## FOR LOOP ##
      ##############
      ace.score <- list()

      for(i in 1:length(snps.ACE)){

        #######################
        ## Handle SNPS probs ##
        #######################
        ## get $lik.anc for internal nodes:
        snp.lik.int <- snps.ACE[[i]]$lik.anc[,2]

        ## get "lik.anc" (ie. states) for terminal nodes:
        snp.lik.term <- snps[,i]

        ## get snp.lik total:
        snp.lik <- c(snp.lik.term, snp.lik.int)

        score <- list()

        ##############
        ## FOR LOOP ##
        ##############
        for(e in 1:nrow(edges)){

          pa <- phen.lik[edges[e,1]]
          pd <- phen.lik[edges[e,2]]
          sa <- snp.lik[edges[e,1]]
          sd <- snp.lik[edges[e,2]]
          bl <- tree$edge.length[e]

          score[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)

        } # end (e) for loop

        ace.score[[i]] <- abs(sum(as.vector(unlist(score))))

      } # end (i) for loop

      ace.score <- as.vector(unlist(ace.score))
      names(ace.score) <- colnames(snps)
      ace.score3.3 <- ace.score

      # head(sort(ace.score, decreasing=TRUE))
      # which(names(sort(ace.score, decreasing=TRUE)) %in% c("2", "4")) # 1132 2306

    }

    ############################################
    ## get values for duplicate snps columns: ##
    ############################################
    if(all.unique == TRUE){
      ace.score.complete <- ace.score
    }else{
      ace.score.complete <- rep(NA, ncol(snps.ori))
      for(i in 1:ncol(snps.unique)){
        ace.score.complete[which(index == i)] <- ace.score[i]
      }
    }

  } # end "subsequent" %in% test

  ## END SUBSEQUENT TEST ##

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  ## END ASSOCIATION TEST(S) ##

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  ###################
  ## return output ##
  ###################
  ## Return vector containing measure of
  ## correlation btw SNPi and PHENi diffs across nodes in tree
  ## for all snps in original snps matrix.
  return(ace.score.complete)

} # end ace.test




################
## get.score3 ##
################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param Pa A numeric value containing either the state,
#' or the probability of the state, of the phenotype at a given \emph{ancestral} node.
#' @param Pd A numeric value containing either the state,
#' or the probability of the state, of the phenotype at a given \emph{descendant} node.
#' @param Sa A numeric value containing either the state,
#' or the probability of the state, of SNPi at a given \emph{ancestral} node.
#' @paran Sd A numeric value containing either the state,
#' or the probability of the state, of SNPi at a given \emph{descendant} node.
#' @param l A numeric value specifying the length of the branch in the phylogenetic tree
#' that joins the ancestral and descendant node.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'

########################################################################

get.score3 <- function(Pa, Pd, Sa, Sd, l){

  score3 <- NULL

  ## CHECKS:
  if(length(Pa) > 1) stop("Pa must be of length one
                                    (i.e., (the probability of)
                          the phenotypic state of ONE ancestral node.")
  if(length(Pd) > 1) stop("Pd must be of length one
                                    (i.e., (the probability of)
                          the phenotypic state of ONE descendant node.")
  if(length(Sa) > 1) stop("Sa must be of length one
                                    (i.e., (the probability of)
                          the SNPi state of ONE ancestral node.")
  if(length(Sd) > 1) stop("Sd must be of length one
                                    (i.e., (the probability of)
                          the SNPi state of ONE descendant node.")

  score3 <- l*(((4/3)*Pa*Sa) +
                 ((2/3)*Pa*Sd) +
                 ((2/3)*Pd*Sa) +
                 ((4/3)*Pd*Sd) -
                 Pa -
                 Pd -
                 Sa -
                 Sd +
                 1)

  return(score3)

} # end get.score3





# ####################
# ## cor.tree.score ##
# ####################
#
# ########################################################################
#
# ###################
# ## DOCUMENTATION ##
# ###################
#
# #' Short one-phrase description.
# #'
# #' Longer proper discription of function...
# #'
# #' @param Pa A numeric value containing either the state,
# #' or the probability of the state, of the phenotype at a given \emph{ancestral} node.
# #' @param Pd A numeric value containing either the state,
# #' or the probability of the state, of the phenotype at a given \emph{descendant} node.
# #' @param Sa A numeric value containing either the state,
# #' or the probability of the state, of SNPi at a given \emph{ancestral} node.
# #' @paran Sd A numeric value containing either the state,
# #' or the probability of the state, of SNPi at a given \emph{descendant} node.
# #' @param l A numeric value specifying the length of the branch in the phylogenetic tree
# #' that joins the ancestral and descendant node.
# #'
# #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
# #' @export
# #' @examples
# #'
# #' ## basic use of fn
# #' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
# #'
#
# ########################################################################
#
# cor.tree.score <- function(Pa, Pd, Sa, Sd, l){
#   out <- NULL
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
##   BLACK & WHITE VERSION:
#
#   (((P1 | S1_branch) / sum(S1_branch_lengths)) –
#    ((P1 | S0_branch) / sum(S0_branch_lengths)))
#   +
#     (((S1 | P1_branch) / sum(P1_branch_lengths)) –
#      ((S1 | P0_branch) / sum(P0_branch_lengths)))
#
#
#   out <- ((Pa*Sa + (1-Pa)*(1-Sa) - Pa*(1-Sa) - (1-Pa)*Sa) +
#             (Pd*Sd + (1-Pd)*(1-Sd) - Pd*(1-Sd) - (1-Pd)*Sd)) * l/2
#
#   return(out)
# } # end corr.tree.score




######################
## get.branch.diffs ##
######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param var A vector containing a variable whose change across edges we want to examine.
#' @param edges A 2-column matrix containing the upstream and downstream nodes
#'  in columns 1 and 2 of a tree's edge matrix, as found in a phylo object's tree$edge slot.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'

########################################################################

## get diffs btw ace prob/likelihood at upstream vs. downstream node
## for a single variable for which you have a value for all terminal and internal nodes.
get.branch.diffs <- function(var, edges){
  ## CHECKS: ##
  ## var should be a vector or have 2 columns summing to 1 or 100:
  if(!is.null(dim(var))){
    if(ncol(var) > 2) warning("var contains more than one discrete variable;
                              selecting first variable.")
    var <- var[,2]
  }
  ## var and tree$edge should contain the same number of inds:
  if(length(var) != (nrow(edges)+1)) stop("var contains more
                                          individuals than tree$edge does.")

  ## ~ FOR LOOP ##
  diffs <- var[edges[,1]] - var[edges[,2]]

  return(as.vector(diffs))
} # end get.branch.diffs



###################
## get.ace.diffs ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param var A vector containing a variable whose change across edges we want to examine.
#' @param edges A 2-column matrix containing the upstream and downstream nodes
#'  in columns 1 and 2 of a tree's edge matrix, as found in a phylo object's tree$edge slot.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'

########################################################################

## get diffs btw ace prob/likelihood at upstream vs. downstream node
## for a single variable for which you have a value for all terminal and internal nodes.
get.ace.diffs <- function(var, edges){
  ## CHECKS: ##
  ## var should be a vector or have 2 columns summing to 1 or 100:
  if(!is.null(dim(var))){
    if(ncol(var) > 2) warning("var contains more than one discrete variable;
                              selecting first variable.")
    var <- var[,2]
  }
  ## var and tree$edge should contain the same number of inds:
  if(length(var) != (nrow(edges)+1)) stop("var contains more individuals than tree$edge does.")

  ## ~ FOR LOOP ##
  diffs <- var[edges[,1]] - var[edges[,2]]

  return(as.vector(diffs))
} # end get.ace.diffs




########################
## get.ancestral.pars ##
########################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Ancestral sequence reconstruction via parsimony
#'
#' A wrapper for the \code{ancestral.pars} function from \emph{ape}. Can perform
#' parsimonious ASR for variables in matrix or vector form.
#'
#' @param var A matrix or vector containing a variable whose state at ancestral nodes we want to infer.
#' @param tree A phylo object containing a phylogenetic tree whose tips contain the same individuals as are
#' in the elements of \code{var}, if \code{var} is a vector,
#' or in the rows of \code{var}, if \code{var} is a matrix.
#'
#' @details Note that the (row)names of \code{var} should match the tip.labels of \code{tree}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @import phangorn ape
#'
#' @export
#'

########################################################################

get.ancestral.pars <- function(var, tree){

  edges <- tree$edge

  #############################
  ## RUN PARSIMONY on MATRIX ##
  #############################

  if(is.matrix(var)){

    snps <- var

    ###############################
    ## get unique SNPs patterns: ##
    ###############################
    temp <- get.unique.matrix(snps, MARGIN=2)
    snps.unique <- temp$unique.data
    index <- temp$index

    if(ncol(snps.unique) == ncol(snps)){
      all.unique <- TRUE
    }else{
      all.unique <- FALSE
    }

    ## work w only unique snps:
    snps.ori <- snps
    snps <- snps.unique

    ##########################################
    ## get SNP states of all internal nodes ##
    ##########################################

    #############
    ## ACCTRAN ##
    #############

    ## pace == ancestral.pars ## (ape)
    ## parsimony version of ace ## based on the fitch algorithm.
    ## ACCTRAN = "ACCelerated TRANsformation"

    ## as.phyDat requires row/col names:
    if(is.null(row.names(snps))){
      if(!is.null(tree$tip.label)) row.names(snps) <- tree$tip.label
    }
    if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))
    ## get levels (ie. 0, 1)
    snps.levels <- unique(as.vector(snps))
    ## returns only unique patterns...
    snps.phyDat <- as.phyDat(as.matrix(snps),
                             type="USER", levels=snps.levels)
    ## get index of all original snps columns to map to unique pattern
    index <- attr(snps.phyDat, "index")

    ## pace == ancestral.pars
    pa.ACCTRAN <- pace(tree, snps.phyDat, type="ACCTRAN")

    ## NOTE: pace  --> diff resuls w MPR vs. ACCTRAN
    # pa.MPR <- pace(tree, snps.phyDat, type="MPR")
    #diffs <- sapply(c(1:length(pa.ACCTRAN)), function(e) identical(pa.MPR[[e]], pa.ACCTRAN[[e]]))

    ###########################################
    ## convert reconstruction back to snps.. ##
    ###########################################
    ## each of the n.ind elements of pa is a matrix w n.snps rows and either:
    ## 2 columns, for the 2 binary SNP states, or
    ## 4 columns, each for the 4 nts possible (acgt)

    # rec <- pa.MPR
    rec <- pa.ACCTRAN

    snps.rec <- vector("list", length(rec))
    ## Handle terminal nodes
    for(i in 1:nrow(snps)){
      snps.rec[[which(row.names(snps) == tree$tip.label[i])]] <- rec[[i]][,2]
    }
    ## Handle internal nodes
    for(i in (nrow(snps)+1):length(rec)){
      snps.rec[[i]] <- rec[[i]][,2]
    }
    ## bind rows together
    snps.rec <- do.call("rbind", snps.rec)
    ## assign rownames for all terminal and internal nodes
    rownames(snps.rec) <- c(rownames(snps), c((nrow(snps)+1):((nrow(snps)*2)-1)))
    colnames(snps.rec) <- c(1:length(snps.phyDat[[1]]))


    ###########################################
    ## get LOCATIONS (branches) of snps subs ##
    ###########################################
    subs.edges <- rep(list(NULL), ncol(snps.rec))
    for(i in 1:ncol(snps.rec)){
      snp <- snps.rec[, i]
      subs.logical <- sapply(c(1:nrow(edges)),
                             function(e)
                               snp[edges[e,1]]
                             ==
                               snp[edges[e,2]])
      ## get indices of all edges containing a substitution
      subs.total <- which(subs.logical == FALSE)
      ## get df of states of ancestor and descendants nodes on these edges
      df <- data.frame(snp[edges[subs.total,1]], snp[edges[subs.total,2]])
      names(df) <- c("anc", "dec")
      ## get indices of all edges w a positive sub (0 --> 1)
      subs.pos <- subs.total[which(df$anc==0)]
      ## get indices of all edges w a negative sub (1 --> 0)
      subs.neg <- subs.total[which(df$anc==1)]

      ## get output list
      subs.edges[[i]] <- rep(list(NULL), 3)
      names(subs.edges[[i]]) <- c("total", "pos", "neg")
      if(length(subs.total) > 0) subs.edges[[i]][["total"]] <- subs.total
      if(length(subs.pos) > 0) subs.edges[[i]][["pos"]] <- subs.pos
      if(length(subs.neg) > 0) subs.edges[[i]][["neg"]] <- subs.neg
    }

    ####################
    ## PLOT to CHECK? ##
    ####################
    ## for SNP1, does it identify the correct/reasonable branches?
    #     edgeCol <- rep("black", nrow(edges))
    #     edgeCol <- replace(edgeCol, subs.edges[[1]][["total"]], "green")
    #
    #     ## plot the i'th character's reconstruction on the tree:
    #     #require(adegenet)
    #     plotAnc(tree, pa.ACCTRAN, i=1,
    #             col=transp(c("red", "royalblue"), 0.75),
    #             cex.pie=0.1, pos=NULL,
    #             edge.color=edgeCol, edge.width=2, use.edge.length=FALSE, type="c")

    ################
    ## Get output ##
    ################

    ## get reconstruction for all original sites
    if(ncol(snps.ori) == ncol(snps.rec)){
      snps.rec.complete <- snps.rec
    }else{
      snps.rec.complete <- matrix(NA, nrow=nrow(snps.rec), ncol=ncol(snps.ori))
      for(i in 1:ncol(snps.rec)){
        snps.rec.complete[, which(index == i)] <- snps.rec[, i]
      }
      rownames(snps.rec.complete) <- rownames(snps.rec)
      colnames(snps.rec.complete) <- colnames(snps.ori)
    }

    ## get sub locations on branches for all original sites
    snps.subs.edges <- subs.edges
    if(ncol(snps.ori) == ncol(snps.rec)){
      snps.subs.edges.complete <- snps.subs.edges
    }else{
      snps.subs.edges.complete <- vector("list", ncol(snps.ori))
      for(i in 1:length(snps.subs.edges)){
        loc <- which(index == i)
        if(length(loc) == 1){
          snps.subs.edges.complete[[loc]] <- snps.subs.edges[[i]]
        }else{
          for(j in 1:length(loc)){
            snps.subs.edges.complete[[loc[j]]] <- snps.subs.edges[[i]]
          }
        }
      }
    }


    ## CHECK-- compare cost from fitch and pace: ##
    ## get n.subs per site by fitch:
    cost <- get.fitch.n.mts(snps, tree)
    ## get n.subs per site by pace:
    cost2 <- sapply(c(1:length(snps.subs.edges.complete)),
                    function(e) length(snps.subs.edges.complete[[e]][["total"]]))
    ## NOTE: cost2 differs somewhat noticeably from original fitch cost
    ## (ie. parsimony shifts distribution toward 1/reduces the weight of the upper tail...)
    ## WHY? Which should we use to get n.subs????????????????????????????????????????????????????????????

    ## Get final output list:
    var.rec <- snps.rec.complete
    subs.edges <- snps.subs.edges.complete

    out <- list("var.rec" = var.rec,
                "subs.edges" = subs.edges)

  }else{ # end matrix (snps) parsimony

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


    #############################
    ## RUN PARSIMONY on VECTOR ##
    #############################

    ## Eg. get PHEN states of all internal nodes

    phen <- var

    #############
    ## ACCTRAN ##
    #############

    phen.ori <- phen
    phen <- as.numeric(as.character(phen))
    ## as.phyDat requires names...
    if(is.null(names(phen))){
      if(!is.null(tree$tip.label)) names(phen) <- tree$tip.label
    }

    ## get levels (ie. 0, 1)
    phen.levels <- unique(phen)
    phen.phyDat <- as.phyDat(as.matrix(phen),
                             type="USER", levels=phen.levels)
    ## pace == ancestral.pars
    rec <- phen.pa.ACCTRAN <- pace(tree, phen.phyDat, type="ACCTRAN")

    ## get reconstruction:
    phen.rec <- list()
    for(i in 1:length(rec)){
      phen.rec[[i]] <- rec[[i]][,2]
    }
    phen.rec <- as.vector(unlist(phen.rec))
    names(phen.rec) <- c(names(phen), c((length(phen)+1):((length(phen)*2)-1)))

    ###########################################
    ## get LOCATIONS (branches) of phen subs ##
    ###########################################

    ## make empty output list
    phen.subs.edges <- rep(list(NULL), 3)
    names(phen.subs.edges) <- c("total", "pos", "neg")

    ## identify if subs occur on each branch:
    phen.subs.logical <- sapply(c(1:nrow(edges)),
                                function(e)
                                  phen.rec[edges[e,1]]
                                ==
                                  phen.rec[edges[e,2]])
    ## get indices of all edges containing a substitution
    phen.subs.total <- which(phen.subs.logical == FALSE)
    ## get df of states of ancestor and descendants nodes on these edges
    df <- data.frame(phen.rec[edges[phen.subs.total,1]], phen.rec[edges[phen.subs.total,2]])
    names(df) <- c("anc", "dec")
    ## get indices of all edges w a positive sub (0 --> 1)
    phen.subs.pos <- phen.subs.total[which(df$anc==0)]
    ## get indices of all edges w a negative sub (1 --> 0)
    phen.subs.neg <- phen.subs.total[which(df$anc==1)]

    ## get output list
    if(length(phen.subs.total) > 0) phen.subs.edges[["total"]] <- phen.subs.total
    if(length(phen.subs.pos) > 0) phen.subs.edges[["pos"]] <- phen.subs.pos
    if(length(phen.subs.neg) > 0) phen.subs.edges[["neg"]] <- phen.subs.neg

    ####################
    ## PLOT to CHECK? ##
    ####################
    ## for SNP1, does it identify the correct/reasonable branches?
    #     edgeCol <- rep("black", nrow(edges))
    #     edgeCol <- replace(edgeCol, phen.subs.edges[["total"]], "green")
    #
    #     ## plot the i'th character's reconstruction on the tree:
    #     #require(adegenet)
    #     plotAnc(tree, phen.pa.ACCTRAN, i=1,
    #             col=transp(c("red", "royalblue"), 0.75),
    #             cex.pie=0.1, pos=NULL,
    #             edge.color=edgeCol, edge.width=2, use.edge.length=FALSE, type="c")

    ################
    ## Get output ##
    ################
    var.rec <- phen.rec
    subs.edges <- phen.subs.edges
    out <- list("var.rec" = var.rec,
                "subs.edges" = subs.edges)

  } # end vector (phen) parsimony

  return(out)

} # end get.ancestral.pars





