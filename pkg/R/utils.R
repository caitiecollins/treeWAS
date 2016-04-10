
#############
## utils.R ##
#############

## useful little functions that get used within other functions

## NOTE: for package release should change all fns to .fns
## here and within all other fns s.t. no documentation required
## (unless we want to release these for public use?)

#####################
## get.unique.snps ##
#####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Get the order of the tip labels of a phylogenetic tree as plotted.
#'
#' Longer proper discription of function...
#'
#' @param snps A matrix of SNPs, potentially containing non-unique patterns.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' @import phangorn


########################################################################
###############################
## get unique SNPs patterns: ##
###############################

get.unique.snps <- function(snps){

  ## Identify unique SNP patterns:
  tab.out <- table.matrix(t(snps))
  snps.unique <- t(as.matrix(tab.out$unique.data))
  index <- tab.out$index
  colnames(snps.unique) <- c(1:ncol(snps.unique))

  out <- list(snps.unique=snps.unique,
              index=index)
  return(out)
} # end get.unique.snps



###################
## get.tip.order ##
###################
## fn getting order of tips as plotted

########################################################################

###################
## DOCUMENTATION ##
###################

#' Get the order of the tip labels of a phylogenetic tree as plotted.
#'
#' Longer proper discription of function...
#'
#' @param tree An object of class phylo containing a tree whose tip order is desired to be known.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
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
#' @import ape Hmisc


########################################################################


get.tip.order <- function(tree){
  require(ape)
  require(Hmisc)
  tree2 <- read.tree(text=write.tree(tree))
  if(all.is.numeric(tree2$tip.label)){
    out <- as.numeric(tree2$tip.label)
    out <- rev(out)
  }else{
    out <- sort(tree2$tip.label)
    out <- sapply(c(1:length(tree$tip.label)), function(e)
      which(out == tree$tip.label[e]))
  }
  return(out)
} # end get.tip.order

## OLD VERSION (pre-2016) --> NAs if tip labs not numeric!
# get.tip.order <- function(tree){
#   require(ape)
#   tree2 <- read.tree(text=write.tree(tree))
#   out <- as.numeric(tree2$tip.label)
#   out <- rev(out)
#   return(out)
# } # end get.tip.order


#########   ###   ###   ###   ###   ###   ###   ###   ###   ###   #########
#########   ###   ###   ###   ###   ###   ###   ###   ###   ###   #########

## NOTE: THESE ARE BEING RENAMED --> SET OF 4 FNS!
## NEED TO KEEP OLD 2 FOR NOW UNTIL YOU CAN SEARCH THROUGH
## ALL YOUR OTHER FNS FOR INSTANCES OF THEIR USE
## AND REPLACE THE OLD FN NAMES W THE NEW!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##################
## .substrRight ##
##################
## truncate character string from right
## ie. keep the LAST n characters of
## (each element of) x

#' @export

.substrRight <- function(x, n){
  sapply(x, function(e)
    substr(e, (nchar(e)-n+1), nchar(e))
  )
} # end .substrRight


#################
## .substrLeft ##
#################
## truncate character string from left
## ie. keep the FIRST n characters of
## (each element of) x

#' @export

.substrLeft <- function(x, n){
  sapply(x, function(e)
    substr(e, 0, n)
  )
} # end .substrLeft

#########   ###   ###   ###   ###   ###   ###   ###   ###   ###   #########

###############
## keepLastN ##
###############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Truncate to keep only the \emph{last} N characters.
#'
#' Truncate an element, or each element of a vector, by
#' removing all but the last N characters of each element.
#'
#' @param x A vector whose element(s) will be truncated.
#' @param n An integer specifying the number of characters to \emph{keep}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################

keepLastN <- function(x, n){
  sapply(x, function(e)
    substr(e, (nchar(e)-n+1), nchar(e))
  )
} # end keepLastN

################
## keepFirstN ##
################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Truncate to keep only the \emph{first} N characters.
#'
#' Truncate an element, or each element of a vector, by
#' removing all but the first N characters of each element.
#'
#' @param x A vector whose element(s) will be truncated.
#' @param n An integer specifying the number of characters to \emph{keep}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################

keepFirstN <- function(x, n){
  sapply(x, function(e)
    substr(e, (nchar(e)-n+1), nchar(e))
  )
} # end keepFirstN

#################
## removeLastN ##
#################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Truncate to remove all of the \emph{last} N characters.
#'
#' Truncate an element, or each element of a vector, by
#' removing the last N characters of each element.
#'
#' @param x A vector whose element(s) will be truncated.
#' @param n An integer specifying the number of characters to \emph{remove}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################

removeLastN <- function(x, n){
  sapply(x, function(e)
    substr(e, 0, (nchar(e)-n))
  )
} # end removeLastN


##################
## removeFirstN ##
##################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Truncate to remove all of the \emph{first} N characters.
#'
#' Truncate an element, or each element of a vector, by
#' removing the first N characters of each element.
#'
#' @param x A vector whose element(s) will be truncated.
#' @param n An integer specifying the number of characters to \emph{remove}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################

removeFirstN <- function(x, n){
  sapply(x, function(e)
    substr(e, n+1, nchar(e))
  )
} # end removeFirstN


#########   ###   ###   ###   ###   ###   ###   ###   ###   ###   #########
#########   ###   ###   ###   ###   ###   ###   ###   ###   ###   #########

##################
## .is.integer0 ##
##################
## mini fn testing for output==integer(0)

#' @export

.is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
} # end .is.integer0()


########################
## Parity-testing fns ##
########################

#' @export

.is.even <- function(x) x %% 2 == 0

#' @export

.is.odd <- function(x) x %% 2 != 0


#########################
## selectBiallelicSNP: ##
#########################
## fn that returns the alternative nt| the nt input

## NOTE-- while this is not inherently the
## definition of a biallelic SNP,
## it is currently suiting my purposes
#### by fulfilling the function of
## ensuring that sites always revert
## back and forth between one state and ONE other,
#### hence never creating any
## triallelic sites or tetralellic sites
## --> binary encoding guaranteed to work fine.

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param x A character vector of length 1 containing a nucleotide to be converted.
#' @param DNA logical; if TRUE (default), uses DNA bases (ACGT), if FALSE, uses RNA bases (ACGU).
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################


selectBiallelicSNP <- function(x, DNA=TRUE){

  out <- NULL

  if(!is.null(x)){

    x <- as.character(x)

    ## for DNA encoding:
    if(DNA==TRUE){
      if(x=="a") out <- "t"
      if(x=="t") out <- "a"
      if(x=="c") out <- "g"
      if(x=="g") out <- "c"

      if(x=="A") out <- "T"
      if(x=="T") out <- "A"
      if(x=="C") out <- "G"
      if(x=="G") out <- "C"
    } # end DNA

    ## for RNA encoding:
    if(DNA==FALSE){
      if(x=="a") out <- "u"
      if(x=="u") out <- "a"
      if(x=="c") out <- "g"
      if(x=="g") out <- "c"

      if(x=="A") out <- "U"
      if(x=="U") out <- "A"
      if(x=="C") out <- "G"
      if(x=="G") out <- "C"
    } # end RNA
  } # end !is.null(x)

  return(out)

} # end selectBiallelicSNP
#######################################################


#####################
## .switch.phen fn ##
#####################

#' @export

.switch.phen <- function(x){
  out <- NULL
  ## A/B coding
  if(x == "A") out <- "B"
  if(x == "B") out <- "A"
  ## R/S coding
  if(x == "R") out <- "S"
  if(x == "S") out <- "R"
  return(out)
} # end .switch.phen

###############
## .getFixed ##
###############

#' @export

.getFixed <- function(locus, posi,
                     exp.char=c("a","t","g","c"),
                     polyThres=0.01){
  vec <- as.character(locus)
  vec[!vec %in% exp.char] <- NA
  N <- sum(!is.na(vec)) # N: number of sequences
  if(N==0 || sum(table(vec)/N >= polyThres )<2){
    return(TRUE) # escape if untyped locus or no SNPs
  }else{
    return(FALSE)
  }
} # end getFixed
