
#############
## utils.R ##
#############

## useful little functions that get used within other functions

## NOTE: for package release should change all fns to .fns
## here and within all other fns s.t. no documentation required
## (unless we want to release these for public use?)

#######################
## get.unique.matrix ##
#######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Get unique rows/columns of a matrix with an index vector.
#'
#' A wrapper for the \code{table.matrix} function that assigns consecutive
#' row or column names to the output matrix's unique rows or columns.
#'
#' @param data A matrix or data.frame, potentially containing
#' non-unique patterns in its rows or columns.
#' @param MARGIN A single integer specifying the array margin to be held fixed.
#' (To get unique \emph{rows}, select \code{MARGIN} = 1;
#' for unique \emph{columns}, select \code{MARGIN} = 2.)
#'
#' @details An extension of the base \code{unique.matrix} function that returns
#' a unique matrix (by removing duplicate rows or columns) and also
#' an index vector containing the indices (row or column numbers),
#' in the matrix composed only of unique rows or columns,
#' to which each row or column in the original matrix corresponds.
#'
#' @return A list with the following elements:
#' \itemize{
#'    \item{\code{index} \item{An index vector containing the indices (row numbers),
#'          in a matrix composed only of unique rows,
#'          to which each row in the original matrix maps.}}
#'    \item{\code{unique.data} \item{A new matrix
#'          containing only the unique rows of the input matrix.}}
#' }
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export


########################################################################
## get unique SNPs column patterns: ##
get.unique.matrix <- function(data, MARGIN=2, silent=TRUE){

  ## Identify unique SNP row/column patterns:
  tab.out <- table.matrix(data, MARGIN=MARGIN)
  unique.data <- as.matrix(tab.out$unique.data)
  index <- tab.out$index

  if(MARGIN == 1){
    row.names(unique.data) <- c(1:nrow(unique.data))
    if(length(unique(index)) == nrow(data)){
      if(silent == FALSE){
      warning("Data inputted was already unique along the selected MARGIN.")
      }
    }
  }
  if(MARGIN == 2){
    colnames(unique.data) <- c(1:ncol(unique.data))
    if(length(unique(index)) == ncol(data)){
      if(silent == FALSE){
      warning("Data inputted was already unique along the selected MARGIN.")
      }
    }
  }

  ## Get output:
  out <- list(unique.data=unique.data,
              index=index)
  return(out)
} # end get.unique.matrix


## assign new name to old name (NOT SURE THIS WORKS... CHECK!)
get.unique.snps <- get.unique.matrix


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
#' @param tree An object of class phylo containing a tree
#' whose tip order is desired to be known.
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
#' @import ape
#' @importFrom Hmisc all.is.numeric

########################################################################

get.tip.order <- function(tree){

  # require(ape)
  # require(Hmisc)

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



##################
## table.matrix ##
##################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Cross-tabulate the rows or columns of a matrix.
#'
#' A version of the base \code{table} function designed for matrices.
#' Taking a matrix as input, \code{table.matrix} returns a contingency table,
#' index vector, and unique matrix.
#'
#' @param data A matrix or data.frame, potentially containing
#' non-unique patterns in its rows or columns.
#' @param MARGIN A single integer specifying the array margin to be held fixed.
#' (To get unique \emph{rows}, select \code{MARGIN} = 1;
#' for unique \emph{columns}, select \code{MARGIN} = 2.)
#'
#' @details To apply this function to the \emph{columns} of a matrix, simply
#' transpose the matrix before executing the command, as in:
#' \code{table.matrix(t(data))}.
#'
#' @return A list with the following elements:
#' \itemize{
#'    \item{\code{table} \item{A contingency table of the counts of the
#'          number of occurrences of each unique row in the matrix.}}
#'    \item{\code{index} \item{An index vector containing the indices (row numbers),
#'          in a matrix composed only of unique rows,
#'          to which each row in the original matrix maps.}}
#'    \item{\code{unique.data} \item{A new matrix
#'          containing only the unique rows of the input matrix.}}
#' }
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#'
#' ## load example data:
#' data("snps.ace")
#' x <- snps.ace
#'
#' ## basic use of fn on rows of x:
#' tab.out <- table.matrix(x)
#'
#' ## apply fn to columns of x:
#' tab.out <- table.matrix(t(x))
#'
#' @export

########################################################################

table.matrix <- function(data, MARGIN=1){

  ## handle MARGIN argument:
  if(is.character(MARGIN)){
    MARGIN <- tolower(MARGIN)
    if(MARGIN %in% c("row", "rows", "r")){
      MARGIN <- 1
    }else{
      if(MARGIN %in%
         c("column", "columns", "col", "cols", "c")){
        MARGIN <- 2
      }else{
        stop("MARGIN argument should be either 1 or 2;
             or, 'rows' or 'columns'.")
      }
    }
  }
  ## for columns, transpose matrix at beginning and end:
  if(MARGIN == 2) data <- t(data)


  ## get df
  if(!is.data.frame(data)){
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  ## concatenate all rows into single elements
  dat <- do.call("paste", c(data, sep = "\r"))
  ## get unique rows
  unique.inds <- !duplicated(dat)
  ## keep only unique rows
  levels <- dat[unique.inds]
  ## make a factor in which each level
  ## is the smushed single-element vector
  ## of each unique row
  cat <- factor(dat, levels = levels)
  ## get n. unique levels
  n.levels <- length(levels(cat))
  ## get the unique index that
  ## each original ind/row should map to:
  map.to <- (as.integer(cat) - 1)
  map.to <- map.to[!is.na(map.to)]
  if(length(map.to)) map.to <- map.to + 1
  ## get the number of inds at each unique level:
  tab <- tabulate(map.to, n.levels)

  ## get output
  if(MARGIN == 1){
    out <- list(table = tab,
                index = map.to,
                unique.data = data[unique.inds, ])
  }
  if(MARGIN == 2){
    out <- list(table = tab,
                index = map.to,
                unique.data = as.data.frame(t(data)[, unique.inds]))
  }

  return(out)

} # end table.matrix


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
    substr(e, 1, n)
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

      if(x == "0") out <- "1"
      if(x == "1") out <- "0"
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

      if(x == "0") out <- "1"
      if(x == "1") out <- "0"
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
  ## 0/1 coding:
  if(x == "0") out <- "1"
  if(x == "1") out <- "0"
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
