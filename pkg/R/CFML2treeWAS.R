

#####################
## cluster.treeWAS ##
#####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' (**Internal fn for running on cluster**)
#'
#' Run treeWAS on processed CFML output. To run on the CLIMB cluster.
#'
#' @param prefix A character string specifying the filename prefix of the dataset to be analysed.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################

#####################
## cluster.treeWAS ##
#####################

# prefix <- "/home/caitiecollins/ClonalFrameML/src/pubMLST/Neisseria/B/Czech/WG/phylip/B_Czech_WG.fas.out"
# prefix <- "B_Czech_WG.fas.out"

CFML2treeWAS <- function(prefix){

  ##########################################################################################################
  ###############
  ## GET DATA: ##
  ###############
  phen <- snps <- tree <- out <- NULL

  phen <- get(load(sprintf('%s.phen_clean.Rdata', prefix)))
  snps <- get(load(sprintf('%s.snps_clean.Rdata', prefix)))
  tree <- get(load(sprintf('%s.tree_clean.Rdata', prefix)))


  ##########################################################################################################
  ## run treeWAS ##
  #################
  out <- treeWAS(snps = snps,
                 phen = phen,
                 tree =  tree,
                 n.snps.sim = 10*ncol(snps),
                 plot.tree = TRUE,
                 filename.plot = sprintf('%s.treeWAS_plots.pdf', prefix))

  print("treeWAS done")

  ##########################################################################################################
  ## Save output ##
  #################
  save(out, file=sprintf('%s.treeWAS_out.Rdata', prefix))


  output <- list("treeWAS.out" = out)

  return(output)

} # end cluster.treeWAS

#####################################################################################################################################
#####################################################################################################################################

#####################################################################################################################################
#####################################################################################################################################





#'
#'
#' ##################
#' ## CFML2treeWAS ##
#' ##################
#'
#' ########################################################################
#'
#' ###################
#' ## DOCUMENTATION ##
#' ###################
#'
#' #' (**Internal fn for running on cluster**)
#' #'
#' #' Convert CFML output to treeWAS output with one function. To run on the CLIMB cluster.
#' #'
#' #' @param prefix A character string specifying the filename prefix of the dataset to be analysed.
#' #'
#' #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' #' @export
#'
#' ########################################################################
#'
#'
#' ###############
#' ## PIPELINE: ##
#' ###############
#'
#' ## PC: ##
#' # (*) Move phenotype (phen.ori.Rdata) to "WG/phylip" folder
#' #     & rename "$prefix"phen_ori.Rdata, eg., "B_Czech_WG.fas.out.phen_ori.Rdata"
#'
#' ## COMMAND LINE (outside cluster): ##
#' # (*) Move to relevant PC folder w cd PATH/TO/FOLDER
#' # (*) Upload phenotype ("$prefix"phen_ori.Rdata) to cluster w/ scp "$prefix"phen_ori.Rdata caitlin@131.251.130.191:.
#'
#' ## CLUSTER: ##
#' # (*) Sign onto cluster w/ ssh caitlin@131.251.130.191
#' # (*) R
#' # (*) library(treeWAS)
#' # (*) run CFML2treeWAS(prefix = "B_Czech_WG.fas.out") # eg...
#' # (*) q()
#' # (*) htop OR ls # check that process is running or that files have been created...
#' # (*) exit
#'
#' ## COMMAND LINE (outside cluster): ##
#' # (*) Save files to PC w/ scp caitlin@131.251.130.191:"$prefix"phen_ori.Rdata ./
#'
#' ##################
#' ## CFML2treeWAS ##
#' ##################
#' ## Run these steps in a single fn... ##
#' ## (1) Run read.CFML
#' ## (2) Data cleaning, subsetting
#' ## (3) Run treeWAS (+ filename.plot arg??)
#' ## (4) Save output + plots
#'
#' # prefix <- "/home/caitiecollins/ClonalFrameML/src/pubMLST/Neisseria/B/Czech/WG/phylip/B_Czech_WG.fas.out"
#' # prefix <- "B_Czech_WG.fas.out"
#' CFML2treeWAS <- function(prefix){
#'
#'   ##########################################################################################################
#'   ## Run read.CFML ##
#'   ###################
#'   dat <- read.CFML(prefix = prefix, plot = FALSE)
#'
#'   ## Isolate (relevant) elements of output:
#'   tree <- dat$tree
#'   snps <- dat$snps
#'
#'   #########################
#'   ## SAVE original data: ##
#'   #########################
#'   save(snps, file=sprintf('%s.snps_ori.Rdata', prefix))
#'   save(tree, file=sprintf('%s.tree_ori.Rdata', prefix))
#'   save(dat, file=sprintf('%s.read.CFML_dat.Rdata', prefix))
#'
#'   rm(dat)
#'
#'
#'   ##########################################################################################################
#'   ## Data cleaning, subsetting ##
#'   ###############################
#'
#'   ####################
#'   ## get phenotype: ##
#'   ####################
#'   phen <- get(load(sprintf('%s.phen_ori.Rdata', prefix)))
#'
#'   ##################
#'   ## Subset data: ##
#'   ##################
#'
#'   ###################
#'   ## Missing phen? ##
#'   ###################
#'   phen.ori <- phen
#'   phen <- as.character(phen)
#'
#'   ## CHECK if phen is disease...
#'   if(any(phen %in% c("invasive (unspecified/other)", "meningitis and septicaemia", "meningitis", "septicaemia"))){
#'     ## unify invasive disease under one heading:
#'     phen <- replace(phen, which(phen == "invasive (unspecified/other)"), "invasive")
#'     phen <- replace(phen, which(phen == "meningitis and septicaemia"), "invasive")
#'     phen <- replace(phen, which(phen == "meningitis"), "invasive")
#'     phen <- replace(phen, which(phen == "septicaemia"), "invasive")
#'     # table(phen)
#'   }
#'
#'   ## re-add names before subsetting:
#'   if(is.null(names(phen))) names(phen) <- names(phen.ori)
#'
#'   ## remove missing/other/NA... from phen and snps:
#'   missing <- c("", "NA", NA, "other")
#'   toRemove <- names(phen)[which(phen %in% missing)]
#'
#'   ##########################
#'   ## Remove missing inds: ##
#'   ##########################
#'   if(length(toRemove) > 0){
#'     ## subset PHEN:
#'     phen <- phen[-which(names(phen) %in% toRemove)]
#'     ## subset SNPS:
#'     snps <- snps[-which(rownames(snps) %in% toRemove), ]
#'     ## subset TREE:
#'     tree <- drop.tip(tree, tip = toRemove)
#'   }
#'
#'
#'   ########################
#'   ## Missing SNPs rows? ##
#'   ########################
#'   toRemove <- NULL
#'
#'   ## Remove any ENTIRELY missing rows...
#'   NA.tab <- sapply(c(1:nrow(snps)), function(e) length(which(is.na(snps[e,]))))
#'   toRemove <- rownames(snps)[which(NA.tab == ncol(snps))]
#'
#'   ##########################
#'   ## Remove missing inds: ##
#'   ##########################
#'   if(length(toRemove) > 0){
#'     ## subset PHEN:
#'     phen <- phen[-which(names(phen) %in% toRemove)]
#'     ## subset SNPS:
#'     snps <- snps[-which(rownames(snps) %in% toRemove), ]
#'     ## subset TREE:
#'     tree <- drop.tip(tree, tip = toRemove)
#'   }
#'
#'   ########################
#'   ## SAVE cleaned data: ##
#'   ########################
#'   save(snps, file=sprintf('%s.snps_clean.Rdata', prefix))
#'   save(phen, file=sprintf('%s.phen_clean.Rdata', prefix))
#'   save(tree, file=sprintf('%s.tree_clean.Rdata', prefix))
#'
#'
#'
#'   ##########################################################################################################
#'   ## run treeWAS ##
#'   #################
#'   # out <- treeWAS(snps = snps,
#'   #                phen = phen,
#'   #                tree =  tree,
#'   #                n.snps.sim = 10*ncol(snps),
#'   #                plot.tree = TRUE,
#'   #                filename.plot = sprintf('%s.treeWAS_plots.pdf', prefix))
#'
#'
#'
#'   ##########################################################################################################
#'   ## Save output ##
#'   #################
#'   # save(out, file=sprintf('%s.treeWAS_out.Rdata', prefix))
#'
#'
#'   output <- list("snps" = snps,
#'                  "phen" = phen,
#'                  "tree" = tree) # "treeWAS.out" = out,
#'
#'   # "read.CFML.dat" = dat
#'
#'   return(output)
#'
#' } # end CFML2treeWAS
#'
#' #####################################################################################################################################
#' #####################################################################################################################################
#'
#' #####################################################################################################################################
#' #####################################################################################################################################
#'




###################
## scriptTreeWAS ##
###################

## Make a new scriptTree-like fn to run read.CFML & treeWAS on the cluster.

#####################################################################################################################################
#####################################################################################################################################

#####################################################################################################################################
#####################################################################################################################################

###################
## BASH COMMANDS ##
###################

################
## QUESTIONS: ##
################
## How to load R packages?
## What does the (sh -c ') part of the nohup line do?
## How to save product of R output into cluster enviro?
## How to write R commands (in regular R format? with each new line being a new command?)

############
## NOTES: ##
############
# "$1" # call the output arg $1
# >log"$1" # write the output to a file called log$1
# 2>err"$1" # write the standard error output to a file called err$1


################################################
# nohup sh -c '
# R

#####################################################################################################################################
#####################################################################################################################################

#####################################################################################################################################
#####################################################################################################################################

################
## scriptTree ##
################
# nohup sh -c '
# seqret -osformat phylip3 "$0" "$0".phylip
# phyml --quiet --no_memory_check -b 0 -v 0 -c 1 -s BEST --no_memory_check -q -i "$0".phylip
# ClonalFrameML "$0".phylip_phyml_tree.txt "$0" "$0".out -kappa 5' "$1" >log"$1" 2>err"$1" &

#####################################################################################################################################
#####################################################################################################################################

#####################################################################################################################################
#####################################################################################################################################

################
## scriptTemp ##
################

## Make a temporary practice script to learn BASH command coding...

# function e { echo $1 }
## e hello --> (prints hello)

# scriptTemp
#
# R
# phen_ori <- get(load("phen_ori.Rdata"))
# print(str(phen_ori))
# phen_new <- phen_ori[1:10]
# save(phen_new, file="./phen_new.Rdata")
# q()

#####################################################################################################################################
#####################################################################################################################################

###################
## scriptTreeWAS ##
###################


## Make ONE R fn that just uses prefix as only argument, runs everything, and saves output.
## Sign onto cluster, run R, run this fn, close R, sign off of cluster, save output to PC.



###########
## args: ##
###########

# $prefix # prefix for read.CFML





#




#



#
