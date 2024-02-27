


################
## pair.tests ##
################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Pairwise tests for categorical phenotypes
#'
#' Internal function to calculate treeWAS 
#' terminal, simultaneous, subsequent tests, 
#' and chi-squared p-values for a given snp across pairs of
#' phenotype levels.  
#'
#' @param x A contingency table (snps[,i] x phen) for score 1 (\code{terminal.test} 
#'          with \code{correct.prop = TRUE}, \code{categorical = TRUE}).
#' @param y A vector of values containing pairwise score 2 (\code{simultaneous.test} 
#'          with \code{categorical = TRUE}) results for snps[,i].
#' @param z A contingency table (snps.rec[,i] x phen.rec) for score 3 (\code{subsequent.test} 
#'          with \code{correct.prop = TRUE}, \code{categorical = TRUE}).
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#' ## Example ##
#' \dontrun{
#' ## basic use of fn
#' out <- pair.tests(x, y, z)
#' }
#' 
#' @importFrom stats chisq.test
#'

########################################################################


pair.tests <- function (x, y, z, 
                        method = "bonf", digits = 3){
  n <- nrow(x)
  N <- n * (n - 1)/2
  df <- data.frame(phen.pair = rep("A", N), stringsAsFactors = FALSE)
  p.chisq <- rep(NA, N)
  phi <- rep(NA, N)
  phi.rec <- rep(NA, N)
  k <- 0
  
  for (a in 1:(n - 1)) {
    for (b in (a + 1):n) {
      k <- k + 1
      ## Get phen pair:
      nom.a <- as.character(rownames(x)[a])
      nom.b <- as.character(rownames(x)[b])
      mat <- matrix(c(x[a, ], x[b, ]), nrow = 2, byrow = TRUE)
      mat.rec <- matrix(c(z[a, ], z[b, ]), nrow = 2, byrow = TRUE)
      df$phen.pair[k] <- paste0(nom.a, " : ", nom.b)
      ## Calculate scores 1, 3, chisq.p values:
      x2 <- suppressWarnings(chisq.test(mat, correct=FALSE))
      x2.rec <- suppressWarnings(chisq.test(mat.rec, correct=FALSE))
      p.chisq[k] <- signif(x2$p.value, digits = digits)
      phi[k] <- signif(sqrt(x2$statistic/sum(mat)), digits = digits)
      phi.rec[k] <- signif(sqrt(x2.rec$statistic/sum(mat.rec)), digits = digits)
    } # end for (b) loop
  } # end for (a) loop
  
  ## Reorder pairwise score 2:
  ox <- rep(NA, length(y))
  for(pp in 1:length(y)){
    noms.pp <- strsplit(names(y)[pp], " : ")[[1]]
    ox[pp] <- which(sapply(c(1:nrow(df)), 
                           function(e) 
                             all(strsplit(df$phen.pair[e], " : ")[[1]] %in% noms.pp)))
  } # end for (pp) loop
  
  df$terminal <- phi
  df$simultaneous <- y[ox]
  df$subsequent <- phi.rec
  df$p.chisq <- p.chisq
  df$p.adj.chisq <- signif(p.adjust(df$p.chisq, method = method), 
                           digits = digits)
  return(df)
  
} # end pair.tests


## eg. output:
# PT[[snps.sig[j]]]
#         phen.pair terminal simultaneous subsequent  p.chisq p.adj.chisq
# 1   chicken : cow    0.876            9      0.852 1.62e-12    4.86e-12
# 2 chicken : human    0.435           -2      0.497 4.53e-04    1.36e-03
# 3     cow : human    0.530            3      0.417 9.25e-06    2.77e-05
