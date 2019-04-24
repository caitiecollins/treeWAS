



############################################################################################################################################




###############
## read.CFML ##
###############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Convert ClonalFrameML output.
#'
#' Convert the output of ClonalFrameML into a form usable within \code{treeWAS}.
#'
#' @param prefix A character string containing the prefix of all file names to be read in.
#'
#' @details The \code{prefix} must be the prefix to three files ending in:
#' (i) "labelled_tree.newick", (ii) "ML_sequence.fasta", (iii) "position_cross_reference.txt".
#'
#' @return read.CFML returns a list containing:
#' (i) \code{tree}: The phylogenetic tree.
#' (ii) \code{snps}: The binary genetic data matrix of polymorphic loci.
#' (iii) \code{snps.rec}: The genetic data reconstruction matrix.
#' (iv) \code{seqs}: The genetic data sequences (polymorphic loci only), a \code{DNAbin} object.
#' (v) \code{index}: The index vector, indicating for each column in \code{seqs}
#' the unique polymorphic column pattern to which it corresponds (0 = non-polymorphic).
#' (vi) \code{n.subs}: The distribution of the number of substitutions per site.
#' Note that all genetic data elements (ii - iv) are returned in expanded form; that is,
#' they contain both unique and duplicate column patterns for all polymorphic loci as indicated in the \code{index} vector.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)

########################################################################

## This function reads the output of CFML and returns a LIST containing:
## the phylogenetic tree, converted to be binary and post-ordered (tree),
## the snps matrix
## the snps.rec matrix
## the set of sequences containing duplicate column patterns (seqs),
## the mapping index of all the original sequences in the set of unique sequence columns (index).
## the distribution of substitutions per site (n.subs),


read.CFML <- function(prefix, tree=NULL, plot=TRUE, suff.length = 2) {

  # require(ape)
  # require(adegenet)

  if(is.null(tree)){
    tree <- read.tree(sprintf('%s.labelled_tree.newick', prefix))
  }else{
    ## if tree is filename, read in:
    if(class(tree) == "character"){
      tree <- read.tree(tree)
    }
  }
  seqs <- read.dna(sprintf('%s.ML_sequence.fasta', prefix), format='fasta')
  mapping <- scan(sprintf('%s.position_cross_reference.txt', prefix), sep=',', quiet=T)
  l <- length(seqs[1,])

  ## check tree (violations cause problems, eg. in phen.sim)
  if(!is.binary.tree(tree)) tree <- multi2di(tree)
  # if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder") ## not sure this makes any difference..

  ## Modify the edge matrix so that it uses the same indices as the fasta file
  labs <- labels(seqs)
  edges <- tree$edge
  treelabs <- c(tree$tip.label, tree$node.label)

  ## If treelabs not of length labs (eg tree$node.label is NULL), fill in remainder w labs from seqs??
  # if(length(treelabs) < length(labs))  treelabs <- c(treelabs, labs[c((length(treelabs)+1):length(labs))])
  if(length(treelabs) > length(labs) | length(treelabs) < length(labs)){
    stopHere <- TRUE
    ## (*** ONLY FOR NON-CFML/non-labelled.tree.newick trees: ***)
    ## If treelabs are a subset of labs (i.e., if seqs are all labelled but one of node or tip labs are missing from tree):
    if(all(treelabs %in% labs) & (length(treelabs) > 0) & (length(labs) == dim(seqs)[1] & (is.null(tree$tip.label) | is.null(tree$node.label)))){
      ## replace empty nodelabs (probably OK)
      if(is.null(tree$node.label)){
        nodelabs <- labs[which(!labs %in% treelabs)]
        treelabs <- c(tree$tip.label, nodelabs)
        if(all(treelabs %in% labs)  & all(labs %in% treelabs)) stopHere <- FALSE
      }
      ## replace empty tiplabs (may be DANGEROUS...)
      # if(is.null(tree$tip.label)){
      #   tiplabs <- labs[which(!labs %in% treelabs)]
      #   treelabs <- c(tiplabs, tree$node.label)
      #   if(all(treelabs %in% labs)  & all(labs %in% treelabs)) stopHere <- FALSE
      # }
    }
    if(stopHere == TRUE){
      stop("The number of individuals (terminal and internal nodes) labelled in the tree
            exceeds the number of individuals (rows) labelled in the sequences.")
    }
  } # end treelabs check

  ## Double check that treelabs and labs match up (order does NOT matter yet):
  if(!all(labs %in% treelabs)){
    stop("Sequence labels do not match tree labels.")
  }
  ## Realign order of labs and treelabs within edge mat:
  for (i in 1:nrow(edges)) for (j in 1:ncol(edges)) edges[i,j] <- which(labs==treelabs[edges[i,j]])

  #Count substitutions for each site
  subs <- rep(0,l) # Number of substitutitions for patterns
  num <- rep(0,l) # Number of times a given pattern is used

  ## SLOW STEP!!!
  # system.time(
  for(i in 1:l) {
    if(length(unique(seqs[,i])) == 2){
      num[i] <- sum(mapping == i) # Only count biallelic sites
      n <- length(which(seqs[edges[, 1], i] != seqs[edges[, 2], i])) ## faster
      subs[i] <- subs[i]+n
    } ## NOTE--only counting biallelic sites for num, and for subs dist (below)...
  } # end for loop
  # ) # end system.time

  ## Build distribution
  dist <- rep(0, max(subs)) # there should be no trailing zeros!
  ## trailing zeros still happening!
  ## bc.
  for (i in 1:max(subs)) dist[i] <- sum(num[which(subs==i)])
  ## assign names:
  names(dist) <- c(1:length(dist))

  ## Plot distribution
  if(plot==TRUE){
    barplot(dist,
            main="Number of substitutions per site",
            names.arg=names(dist),
            ylab="Frequency",
            col=transp("royalblue", alpha=0.5))
  }



  ## Expand seqs with mapping/index:
  seqs <- seqs[,mapping]

  ## Convert DNAbin object:
  snps.rec <- DNAbin2genind(seqs, polyThres=0)
  snps.rec <- snps.rec@tab

  ## Get binary loci:
  snps.rec.bin <- get.binary.snps(snps.rec)

  ## Expand snps:
  snps.rec <- snps.rec.bin

  ## Get snps only:
  N <- nrow(snps.rec)-tree$Nnode
  toKeep <- 1:N
  snps <- snps.rec[toKeep, ]


  ## NEW ##

  ## keep original names...
  ########################
  ## get ORIGINAL positions of snps.loci: ##
  ## name index according to ORIGINAL indices:
  indx <- mapping
  names(indx) <- c(1:length(indx))
  ## get POLYMORPHIC indices:
  ix <- indx[indx != 0]

  ####    ####    ####    ####    ####
  ## ENTER SNP.locus POSITION:
  if(is.null(suff.length)) suff.length <- 0
  l.ori <- removeLastN(colnames(snps), suff.length)
  l.ori <- as.numeric(l.ori)
  ## and SEQS positions:
  l.ori.sq <- 1:ncol(seqs)

  ## Identify locus original name among polymorphic indices:
  l <- as.numeric(names(ix[l.ori])) # snps
  l.sq <- as.numeric(names(ix[l.ori.sq])) # seqs
  # head(l, 20) # original positions
  ## Check (should be no NAs):
  if(length(which(is.na(l))) > 0){
    warning("Oops, there should be no NAs in the set of original positions, but NAs have been generated.")
  }

  ## Set originalPositions.allele as snps colnames:
  colnames(snps) <- paste(l, keepLastN(colnames(snps), suff.length), sep="")
  colnames(snps.rec) <- colnames(snps)
  colnames(seqs) <- l.sq
  ########################


  ## NOTE: BOTH seqs and snps are being returned in EXPANDED form (POLYMORPHIC only).
  ## This means that the returned index element is not needed for expansion, but instead
  ## serves as a record of the expansion/duplication of polymorphic columns.

  out <- list(tree = tree,
              snps = snps,
              snps.rec = snps.rec,
              seqs = seqs,
              index = mapping,
              n.subs = dist)

  return(out)

} # end read.CFML





#

##

####

########

#################################################################################################################################







############################################################################################################################################




#######################
## get.original.loci ##
#######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' \code{(read.CFML+)} Get original sequence positions of polymorphic loci.
#'
#' If you ran \code{read.CFML} on ClonalFrameML output before running \code{treeWAS},
#' this function can be used to identify the original sequence positions of your polymorphic loci.
#' E.g., If \code{treeWAS} identified loci "1417.a" and "2017.g" as significant, \code{get.original.loci}
#' can identify corresponding sequence positions "1165743" and "1741392" and return
#' flanking sequence segments.
#'
#' @param seqs A \code{DNAbin} object containing the original sequences
#' input into ClonalFrameML (see details).
#' @param dat An object containing the output of the \code{read.CFML} function.
#' @param sig.snps.names A character vector containing the names of
#' polymorphic loci whose original sequence positions you desire (see details).
#' @param n.bp An integer specifying the desired length of the flanking
#' sequence to be returned; by default, 50 (see details).
#' @param suff.length An integer specifying the suffix length
#' of \code{snps} elements; by default, 2 (see details).
#' @param csv A logical indicating whether to save the results as a CSV file.
#' @param csv.prefix An optional character vector specifying a directory and
#' filename prefix for the CSV file (if \code{csv=TRUE}); default name/suffix, "sig_loci.csv".
#' \emph{Please be careful: Any existing file of that name will be overwritten!}
#' @param NA.thresh A number between 0 and 1 indicating the max allowable
#' proportion of NAs that the output sequence fragments can contain.
#' (if a sequence fragment from row 1 exceeds this threshold,
#' a sufficiently complete sequence fragment will be sought in subsequent rows); by default, 0.2.
#'
#' @details \strong{seqs} must contain ClonalFrameML \emph{input*},
#' which can be read in from fasta with \code{read.dna("FILENAME.fasta", format="fasta")}
#' (*not the ClonalFrameML output file "ML_sequence.fasta" or the \code{seqs} element of \code{read.CFML} output).\cr\cr
#' \strong{sig.snps.names} can contain any set of \code{colnames(snps)}, for example,
#' the set of significant loci identified by \code{treeWAS} (\code{out$treeWAS.combined$treeWAS.combined}).\cr\cr
#' \strong{n.bp} specifies the total length of flanking sequence
#' (drawn from the first row of \code{seqs} only),
#' half of which will be on either side of each locus in \code{sig.snps.names}.
#' Each such sequence will be of total length \code{n.bp+1}, arranged (e.g., with \code{n.bp = 50}) as:\cr
#' <---25bp---><locus.i><---25bp--->.\cr\cr
#' \strong{suff.length} tells the \code{removeLastN} function how many characters are used to specify
#' the allele in \code{sig.snps.names} and \code{colnames(snps)}. For names of the form:
#' "1234.a", \code{suff.length = 2} (note that the decimal counts as a character).
#' If \code{snps} names are purely numeric with no alleles indicated
#' (i.e., they already match names in \code{seqs}), then set \code{suff.length = 0}.
#'
#' @return \code{get.original.loci} returns a list containing:
#' \enumerate{
#' \item \code{loci}: The original sequence positions for all polymorphic loci in \code{seqs}.
#' \item \code{loci.sig}: The original sequence positions for all polymorphic loci in \code{sig.snps.names}.
#' \item \code{seq.sig}: A list of length \code{sig.snps.names} containing sequence fragments of length \code{n.bp}.
#' }
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' @importFrom ape read.dna
#' @importFrom utils write.csv

########################################################################

# ################################################
# fasta <- "/home/caitiecollins/ClonalFrameML/src/pubMLST/Gono/Grad2014_WG.fas"
# prefix <- "/home/caitiecollins/ClonalFrameML/src/pubMLST/Gono/Grad2014_WG.fas.out"
#
# ## read in original fasta sequence:
# seqs <- read.dna(fasta, format="fasta")
#
# ## load read.CFML_dat.Rdata
# dat <- get(load(sprintf('%s.read.CFML_dat.Rdata', prefix)))
#
# ## get sig.snps.names:
# out <- get(load("/home/caitiecollins/ClonalFrameML/src/pubMLST/Gono/CRO/Grad2014_WG.fas.out.treeWAS_out.Rdata"))
# # out <- get(load("/home/caitiecollins/ClonalFrameML/src/pubMLST/Gono/CRO/Grad2014_WG.fas.out.treeWAS_PA_out.Rdata"))
# sig.snps.names <- out$treeWAS.combined$treeWAS.combined
#
# foo <- get.original.loci(seqs, dat, sig.snps.names, n.bp=40, csv=T, csv.prefix="/home/caitiecollins/ClonalFrameML/src/pubMLST/Gono/CRO/Grad2014_CRO")
################################################
################################################
# fasta <- "/media/caitiecollins/Seagate Backup Plus Drive/ClonalFrameML/src/pubMLST/Neisseria/B/penicillin_range/WG/phylip/B_penicillin_range_WG.fas"
# prefix <- "/media/caitiecollins/Seagate Backup Plus Drive/ClonalFrameML/src/pubMLST/Neisseria/B/penicillin_range/WG/phylip/B_penicillin_range_WG.fas.out"
#
# ## read in original fasta sequence:
# seqs <- read.dna(fasta, format="fasta")
#
# ## load read.CFML_dat.Rdata
# dat <- get(load(sprintf('%s.read.CFML_dat.Rdata', prefix)))
#
# ## get sig.snps.names from treeWAS output:
# ## core SNPs results:
# out <- get(load(sprintf('%s.treeWAS_out.Rdata', prefix)))
# ## accessory gene presence/absence results:
# # out <- get(load("/media/caitiecollins/Seagate Backup Plus Drive/ClonalFrameML/src/pubMLST/Neisseria/B/penicillin_range/WG/phylip/B_penicillin_range_WG.fas.out.treeWAS_PA_out.Rdata"))
# sig.snps.names <- out$treeWAS.combined$treeWAS.combined
#
# foo <- get.original.loci(seqs, dat, sig.snps.names, n.bp=40, csv=T, csv.prefix=prefix)
################################################

get.original.loci <- function(seqs, dat, sig.snps.names, n.bp = 50, suff.length = 2, csv = TRUE, csv.prefix = NULL, NA.thresh = 0.2){

  ## Extract elements from read.CFML dat:
  index <- dat$index
  snps <- dat$snps
  # seqs <- dat$seqs

  ## (if necessary) get ORIGINAL positions of snps.loci: ##
  ## name index according to ORIGINAL indices:
  indx <- index
  names(indx) <- c(1:length(indx))
  ## get POLYMORPHIC indices:
  ix <- indx[which(indx != 0)]

  ####    ####    ####    ####    ####
  ## CHECK---old snps.names (!= seqs) or new (matching seqs)?
  if(max(as.numeric(names(ix))) == max(as.numeric(unique(removeLastN(colnames(snps), suff.length))))){
  # which.max(mapping)
  # identical(as.numeric(unique(removeLastN(colnames(snps), suff.length))), as.numeric(unique(colnames(seqs))))

    ## Get set of ORIGINAL POSITIONS (total or for sig.loci only if provided): ##
    sig.snps <- sapply(c(1:length(sig.snps.names)), function(e) which(colnames(snps) == sig.snps.names[e]))
    l.sig <- l[sig.snps]
    names(l.sig) <- paste(l[sig.snps], keepLastN(sig.snps.names, suff.length), sep="")

  }else{
    ####    ####    ####    ####    ####
    ## ENTER SNP.locus POSITION:
    if(is.null(suff.length)) suff.length <- 0
    # l.ori <- 1:ncol(snps)
    l.ori <- removeLastN(colnames(snps), suff.length)
    l.ori <- as.numeric(l.ori)

    ## Identify sig.locus original name among polymorphic indices:
    l <- as.numeric(names(ix[l.ori]))
    # head(l, 20) # original positions
    ## Check (should be no NAs):
    if(length(which(is.na(l))) > 0){
      warning("Oops, there should be no NAs in the set of original positions, but NAs have been generated.")
    }

    ## Get set of ORIGINAL POSITIONS (total or for sig.loci only if provided): ##
    sig.snps <- sapply(c(1:length(sig.snps.names)), function(e) which(colnames(snps) == sig.snps.names[e]))
    # sig.snps <- which(colnames(snps) %in% sig.snps.names)
    l.sig <- l[sig.snps]
    names(l.sig) <- paste(l[sig.snps], keepLastN(sig.snps.names, suff.length), sep="")
  }


  ###############################################################
  ## get ALLELES of ORIGINAL position of sig.locus +/- n.bp/2: ##
  ###############################################################
  if(is.null(n.bp)) n.bp <- 50
  bp.back <- floor(n.bp/2)
  bp.fwd <- ceiling(n.bp/2)

  ####################################
  SEQ.l <- list()
  miss <- c("-", " ", "", NA, "NA", "na", ".", "99")
  for(i in 1:length(l.sig)){
    ## get fragment start & end:
    start <- l.sig[i]-bp.back
    end <- l.sig[i]+bp.fwd

    ## get row 1 sequence:
    rowN <- 1
    seq <- as.character(seqs[rowN, start:end])
    ## check that proportion NAs is lower than threshold:
    while((length(which(seq %in% miss)) > length(seq)*NA.thresh) & rowN < nrow(seqs)){
      rowN <- rowN+1
      seq.ori <- seq
      seq <- as.character(seqs[rowN, start:end])
      ## keep the more complete sequence:
      if(length(which(seq %in% miss)) > length(which(seq.ori %in% miss))){
        seq <- seq.ori
      }
    } # end while loop
    if((length(which(seq %in% miss)) > length(seq)*NA.thresh)){
      cat("Sequence fragment ", i, ": ", length(which(seq %in% miss))/length(seq), "% missing.", sep="")
    }

    ## Paste seq together:
    SEQ <- toupper(seq)
    SEQ <- paste0(SEQ, collapse = "")
    SEQ.l[[i]] <- SEQ

  } # end for loop
  names(SEQ.l) <- l.sig
  # SEQ.l[[1]]
  ####################################

  ## save loci as csv? ##
  if(csv==TRUE){
    df <- data.frame("SNP.locus" = sig.snps.names,
                     "Original.position" = as.vector(unlist(l.sig)),
                     "Original.locus" = as.vector(unlist(names(l.sig))),
                     "Original.seq" = as.vector(unlist(SEQ.l)))

    if(is.null(csv.prefix)){
      write.csv(df, file='sig_loci.csv')
    }else{
      write.csv(df, file=sprintf('%s_sig_loci.csv', csv.prefix))
    }
  } # end csv


  ## get OUTPUT (real positions (total & for sig.loci) & seqs (for sig.loci)): ##
  foo <- list("loci" = l,
              "loci.sig" = l.sig,
              "seq.sig" = SEQ.l)

  ## return output:
  return(foo)

} # end get.original.loci

#################################################################################################################################



