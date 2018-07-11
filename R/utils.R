
#############
## utils.R ##
#############

## useful little functions that get used within other functions

## NOTE: for package release should change all fns to .fns
## here and within all other fns s.t. no documentation required
## (unless we want to release these for public use?)


################################################################################


#############
## memfree ##
#############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Get the current amount of available memory.
#'
#' Function to determine how much memory (in GB) is currently available for use
#' on your PC.
#'
#'
#' @param OS A character string indicating the operating system of the machine in question.
#' Can be one of "Windows", "Mac" (or "Darwin"), or "Linux". If OS is NULL (the default),
#' OS will be set to Sys.info()["sysname"].
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export


########################################################################

#############
## memfree ##
#############
memfree <- function(OS = NULL){

  mem <- NULL

  if(is.null(OS)) OS <- Sys.info()["sysname"]
  OS <- tolower(OS)

  if(OS == "darwin") OS <- "mac"

  if(OS == "Windows"){
    ## WORKS ON WINDOWS: ##
    mem <- as.numeric(shell("wmic OS get FreePhysicalMemory",
                                intern=TRUE)[2])/1000000 # (KB -> GB)
  }

  if(OS == "linux"){
    ## WORKS ON LINUX: ##
    mem <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",
                                 intern=TRUE))/1000000 # (KB -> GB)
  }

  if(OS == "mac"){
    ## WORKS ON MAC: ##
    mem <- system("vm_stat", intern=TRUE)[2] # [2] = "Pages free:    #."
    ## Remove text before and (.) after
    mem <- removeLastN(mem, 1)
    mem <- removeFirstN(mem, nchar("Pages free:"))
    ## Convert to numeric:
    mem <- as.numeric(mem)
    ## Must convert: 1 Page = 4096 bytes --> GB:
    mem <- mem*4096/1000000000
  }

  # mem
  ## Units = GB

  return(mem)

} # end memfree


################################################################################









#####################
## get.binary.snps ##
#####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Reduce a genetic data matrix to only necessary columns.
#'
#' Function to reduce a genetic data matrix containing multiple columns per locus
#' to one column for each binary locus and N columns for each N-allelic non-binary locus.
#'
#'
#' @param snps A genetic data matrix.
#'
#' @details This funtion identifies the number of alleles at each locus by assuming that
#' the allele of each column is contained in the last two characters of each column name.
#' We recommend that the columns of \code{snps} be labelled using the following four suffixes:
#' ".a", ".c", ".g", ".t" (e.g., "Locus_123243.a", "Locus_123243.g").
#' If you are using an alternative naming convention,
#' but the allele is also always being denoted using the last two characters
#' (e.g., "Locus_123243_1", "Locus_123243_2"),
#' the function will still work if you set the argument \code{force = TRUE}.
#' Please also be careful not to accidentally remove any purposeful duplications with repeated names;
#' for example, if you have deliberately duplicated unique columns
#' (e.g., by expanding according to an index returned by ClonalFrameML).
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export


########################################################################

get.binary.snps <- function(snps, force=FALSE){

  suffixes <- keepLastN(colnames(snps), 2)
  suffixes <- unique(suffixes)

  if(all(suffixes %in% c(".a", ".c", ".g", ".t")) | force == TRUE){

    noms <- removeLastN(colnames(snps), 2)
    tab <- table(noms) ## note: tab in character order of colnames!

    ## if any non-binary loci, set these columns aside and keep all N of them:
    if(length(which(tab != 2))){

      ## For NON-BINARY loci, KEEP all columns, PLUS EVERY OTHER BINARY column:
      toKeep <- names(which(tab != 2))
      cols.toKeep <- which(noms %in% toKeep)
      ## make a cols.toKeep remaining vector:
      cols.nonBin <- cols.toKeep

      ## Identify number of columns to be in reduced matrix:
      cols <- c(1:ncol(snps))
      cols <- cols[-cols.toKeep]
      ncol.red <- (length(cols)/2) + length(cols.toKeep)

      ## make reduced snps matrix...
      ## (with one column for each binary SNP and N columns for the non-binary loci):

      COLS <- list()
      COLS[[1]] <- c(1:(cols.toKeep[1] - 1))
      COLS[[1]] <- COLS[[1]][seq(1, length(COLS[[1]]), 2)]
      ## Isolate non-binary columns to add here:
      cols.toAdd <- which(noms %in% noms[cols.nonBin[1]])
      ## Append these columns to current sequence:
      COLS[[1]] <- c(COLS[[1]], cols.toAdd)
      ## Remove cols.toAdd from cols.nonBin:
      cols.nonBin <- cols.nonBin[-which(cols.nonBin %in% cols.toAdd)]

      ## FOR LOOP: ##
      if(length(toKeep) > 1) for(i in 2:length(toKeep)){
        ## Isolate non-binary columns to add here:
        cols.toAdd <- which(noms %in% noms[cols.nonBin[1]])
        ## Get sequence from end of last COLS element to start of cols.toAdd
        from <- COLS[[(i-1)]]
        from <- (from[(length(from))] + 1)
        to <- (cols.toAdd[1] - 1)
        ## Append non-binary columns to binary sequence:
        ## (Unless there are two non-binary columns in a row!)
        if(to > from){
          COLS[[i]] <- seq(from, to, 2)
          COLS[[i]] <- c(COLS[[i]], cols.toAdd)
        }else{
          COLS[[i]] <- cols.toAdd
        }
        ## Remove cols.toAdd from cols.nonBin:
        cols.nonBin <- cols.nonBin[-which(cols.nonBin %in% cols.toAdd)]
      } # end for loop

      COLS.ori <- COLS
      COLS <- as.vector(unlist(COLS))

      ## Finally, add seq to max new ncol(snps)...
      if(length(COLS) < ncol.red){
        from <- (COLS[length(COLS)] + 1)
        COLS <- c(COLS, seq(from=from, by=2, length.out=(ncol.red - length(COLS))))
      }

      ## Get the appropriately reduced set of loci:
      snps <- snps[, COLS]

    }else{
      ## For only BINARY loci, REMOVE 2nd column:
      ## Keep every other snps column:
      toKeep <- seq(1, ncol(snps), 2)
      snps <- snps[, toKeep]
    }
  }else{
    warning("This function requires column names using these suffixes:
            '.a', '.c', '.g', '.t' (e.g., 'Locus_123243.a').
            If there are redundant columns, please remove these by hand,
            or see ?get.binary.snps for more.")
  }
  return(snps)
} # end get.binary.snps

################################################################################


##############
## set.args ##
##############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Set a list of arguments.
#'
#' Function to set a list of arguments without having to remove commas.
#' Useful for troubleshooting. For example, if attempting to run a function
#' (particualrly one with many arguments) line by line,
#' \code{set.args} can be used to set a list of arguments in one go, by copying a
#' comma-separated set of arguments from an existing function call or a new call to \code{args(fn)}.
#'
#'
#' @param args A named list of arguments.
#' @param envir The environment in which these arguments will set.
#'
#' @details Please note that unless the \code{envir} argument is changed from its default (\code{sys.frame}),
#' any arguments set with \code{set.args} will \emph{over-ride} any values currently assigned to those names.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export


########################################################################

set.args <- function(args, envir=sys.frame(which=0L)){

  if(!is.list(args)) args <- as.list(args)
  if(length(names(args)) != length(args)) stop("All elements of the args list must be named.")

  sapply(c(1:length(args)), function(e) assign(names(args)[e], args[[e]], envir=envir))

  return(NULL)
} # end set.args




###############
## ggplot.bg ##
###############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Mimic ggplot2 Background
#'
#' Get an imitation ggplot2-style background for plots made outside ggplot2
#'
#'
#' @details This function must be sandwiched between two instances
#' of the function used to generate the (foreground) plot
#' to which you are hoping to add this background.
#' \emph{Before} running the \code{ggplot.bg} function, you need to run your plot function
#' so that \code{ggplot.bg} knnows how to set the axes.
#' \emph{After} running the \code{ggplot.bg} function, you need to run your plot function
#' again \emph{with the added argument} \code{add=TRUE}
#' so that your plot can be overlayed on top of the background.
#'
#' @param bg The background colour, by default ``lightgray'' with 50\% transparency.
#' @param x.ax A logical specifying whether to re-draw the x-axis.
#' @param y.ax A logical specifying whether to re-draw the y-axis.
#' @param box A logical specifying whether to draw a box around the plotting area.
#' @param grid A logical specifying whether to draw a grid across the background within the plotting area.
#' @param grid.col The color of the gridlines, ``white'' by default. Only used if grid is set to TRUE.
#' @param grid.nx An optional integer to specify the number of gridlines to be drawn along the x-axis.
#' @param grid.ny An optional integer to specify the number of gridlines to be drawn along the y-axis.
#' @param grid.lwd An integer specifying the lwd (line weight) of the gridlines; by default, set to 1.
#' @param grid.lty An integer specifying the line type to be used for the gridlines; by default, set to 1 (i.e., solid lines).
#'
#' @importFrom adegenet transp
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export


########################################################################

## importFrom adegenet transp

ggplotbg <- function(bg=transp("lightgray", 0.5),
                      x.ax=FALSE, y.ax=FALSE, box=TRUE,
                      grid=TRUE, grid.col="white",
                      grid.nx=NULL, grid.ny=NULL, grid.lwd=1, grid.lty=1){

  ## get user plotting parameters:
  lim <- par("usr")
  rect(lim[1],  lim[3], lim[2], lim[4], col="white")

  par(bg=bg)
  # rect(lim[1],  lim[3], lim[2], lim[4], col=bg)

  ## add axes back
  if(!is.null(x.ax)) x.ax
  if(!is.null(y.ax)) y.ax

  ## add grid:
  if(grid == TRUE){
    grid(nx=grid.nx, ny=grid.ny, col=grid.col, lwd=1, lty=1)
  }

  if(box == TRUE) box()   ## and the plot frame

  return(NULL)
} # end ggplotbg

ggplot.bg <- ggplotbg

################################################################################



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

    if(!is.null(row.names(data))){
      row.names(unique.data) <- row.names(data)[!duplicated(index)] ## TO DO: CHECK (MAKE SURE THIS IS OK)
    }else{
      row.names(unique.data) <- c(1:nrow(unique.data))
    }
    if(length(unique(index)) == nrow(data)){
      if(silent == FALSE){
      warning("Data inputted was already unique along the selected MARGIN.")
      }
    }
  }
  if(MARGIN == 2){
    if(!is.null(colnames(data))){
      colnames(unique.data) <- colnames(data)[!duplicated(index)]
    }else{
      colnames(unique.data) <- c(1:ncol(unique.data)) ## TO DO: CHECK! (NOT SURE WHY I DID THIS INSTEAD OF JUST USING THE FIRST INSTANCE.. (PROBLEMATIC ANYWHERE? MISLEADING?))
      ## PLUS -- WOULD IT BE BETTER HERE TO USE THE NON-DUPLICATED INDEX RATHER THAN JUST NUMBERING 1:NCOL?
    }
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
# @examples
#
# ## load example data:
# data("snps.ace")
# x <- snps.ace
#
# ## basic use of fn on rows of x:
# tab.out <- table.matrix(x)
#
# ## apply fn to columns of x:
# tab.out <- table.matrix(t(x))
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

      if(x == FALSE) out <- TRUE
      if(x == TRUE) out <- FALSE
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

      if(x == FALSE) out <- TRUE
      if(x == TRUE) out <- FALSE
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
  ## T/F coding:
  if(x == FALSE) out <- TRUE
  if(x == TRUE) out <- FALSE
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
