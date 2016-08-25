

############################
## NEW SCORE 3 (integral) ##
############################

# Pa <- 0.2
# Pd <- 0.8
# Sa <- 0.3
# Sd <- 0.7
# l <- 1
# x <- l

Pa <- 0
Pd <- 1
Sa <- 0
Sd <- 1
l <- 1
x <- l


Px <- (((Pa*(l - x)) + (Pd*x))/l)

Sx <- (((Sa*(l - x)) + (Sd*x))/l)

Cx <- Px*Sx + (1 - Px)*(1 - Sx) - Px*(1 - Sx) - (1 - Px)*Sx

# score3 <- integrate(f = Cx, lower = 0, upper = l)


## Long-form Cx function ##
get.Cx <- function(x=1, Pa=0, Pd=1, Sa=0, Sd=1, l=1){

  Cx <- NULL

  Cx <- (((Pa*(l - x)) + (Pd*x))/l) * (((Sa*(l - x)) + (Sd*x))/l) +
        (1 - (((Pa*(l - x)) + (Pd*x))/l)) * (1 - (((Sa*(l - x)) + (Sd*x))/l)) -
        (((Pa*(l - x)) + (Pd*x))/l) * (1 - (((Sa*(l - x)) + (Sd*x))/l)) -
        (1 - (((Pa*(l - x)) + (Pd*x))/l)) * (((Sa*(l - x)) + (Sd*x))/l)

  return(Cx)
} # end get.Cx

score3 <- integrate(f = get.Cx, lower = 0, upper = l)
score3


##################
## MATHEMATICA: ##
##################

## Need single-character variables/parameters:
a <- Pa
b <- Pd
c <- Sa
d <- Sd

## SUBMIT THIS to mathematica (without the "Cx <- "): ##
Cx <- (((a*(l - x)) + (b*x))/l) * (((c*(l - x)) + (d*x))/l) +
      (1 - (((a*(l - x)) + (b*x))/l)) * (1 - (((c*(l - x)) + (d*x))/l)) -
      (((a*(l - x)) + (b*x))/l) * (1 - (((c*(l - x)) + (d*x))/l)) -
      (1 - (((a*(l - x)) + (b*x))/l)) * (((c*(l - x)) + (d*x))/l)

## SIMPLIFIED INPUT: ##
integrate(((a*(l - x) + b*x)/l)*((c*(l - x) + d*x)/l) +
            (1 - (a*(l - x) + b*x)/l)* (1 - (c*(l - x) + d*x)/l) -
            ((a*(l - x) + b*x)/l)* (1 - (c*(l - x) + d*x)/l) -
            (1 - (a*(l - x) + b*x)/l)*((c*(l - x) + d*x)/l),
            lower = 0, upper = x)

## SOLN: ##
((-1 + 2*a)*(-1 + 2*c)*l^2*x - (-a + b - c + 4*a*c - 2*b*c + d - 2*a*d)*l*x^2 + (4*(a - b)*(c - d)*x^3)/3)/l^2
###########

## Re-writing soln in terms of Pa, Pd, Sa, Sd: ##
score3.i <- ((-1 + 2*Pa)*(-1 + 2*Sa)*l^2*x -
               (-Pa + Pd - Sa + 4*Pa*Sa - 2*Pd*Sa + Sd - 2*Pa*Sd)*l*x^2 +
               (4*(Pa - Pd)*(Sa - Sd)*x^3)/3) / l^2
score3.i

## Re-writing SIMPLIFIED soln, WHERE x = l... ##
score3 <- l*(((4/3)*Pa*Sa) + ((2/3)*Pa*Sd) + ((2/3)*Pd*Sa) + ((4/3)*Pd*Sd) - Pa - Pd - Sa - Sd + 1)


get.score3 <- function(Pa, Pd, Sa, Sd, l){
  score3 <- NULL
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

l <- 1
Pa <- 0.2
Pd <- 0.49
Sa <- 0.6
Sd <- 0.9

get.score3(Pa, Pd, Sa, Sd, l) # - 0.126

l <- 1
Pa <- 0.49
Pd <- 0.2
Sa <- 0.6
Sd <- 0.9

get.score3(Pa, Pd, Sa, Sd, l) # - 0.184

l <- 1
Pa <- 0
Pd <- 1
Sa <- 0
Sd <- 1

get.score3(Pa, Pd, Sa, Sd, l) # 0.333333


##################
## CONTOUR PLOT ##
##################
library(lattice)

################################################
## Run f2 function to get unstructured output ##
################################################

res <- f2(n=5000)

contourplot(score ~ (Pa - Pd)*(Sa - Sd), data=res)
     # xlab="(Pa - Pd) * (Sa - Sd)", ylab="Score")

scoreSpecial


# l <- 1
# dx <- deriv(score3.i, "x")
# x <- 0:1
# eval(dx)

#####################################################################################################################


## Run the following two functions (i.e., run all the code below) to examine the results of our
## third association testing score.

## The first function's output examines the score that results when the ancestral state of BOTH
## the phenotype and SNPi is a variable in [0,1] and the descendant state of BOTH is (1 - ancestor).

## The second function's output examines the score that results when all four states--
## Pa (phenotype of ancestor), Pd (phenotype of descendant), Sa (SNPi of ancestor), Sd (SNPi of descendant)--
## are random variables drawn from [0,1] with no relationship to each other.

## Note that when I say "the state" above, this may also be interpreted as "probability of the state",
## if we are working with a probabilistic inference of a binary variable.

## Note also that all scores below are calculated only for a single branch, always of length 1.

#######################
## Structured output ##
#######################
## (descendant = 1 - ancestor)

########
## f1 ##
########
f1 <- function(x,y){

  out <- NULL

  a <- b <- c <- d <- c(x,y)
  df <- expand.grid(a,b,c,d)
  names(df) <- c("Pa", "Pd", "Sa", "Sd")

  score <- Pa <- Pd <- Sa <- Sd <- list()

  for(i in 1:nrow(df)){
    pa <- Pa[[i]] <- df[i,1]
    pd <- Pd[[i]] <- df[i,2]
    sa <- Sa[[i]] <- df[i,3]
    sd <- Sd[[i]] <- df[i,4]
    l <- 1

    #     score[[i]] <- abs(((pa*sa + (1-pa)*(1-sa) - pa*(1-sa) - (1-pa)*sa) +
    #                          (pd*sd + (1-pd)*(1-sd) - pd*(1-sd) - (1-pd)*sd)) * l/2)

    score[[i]] <- get.score3(Pa=pa, Pd=pd, Sa=sa, Sd=sd, l=l)
  }

  score <- as.vector(unlist(score))

  out <- cbind(df, score)

  return(out)
} # end f1

##############################################
## Run f1 function to get structured output ##
##############################################

out <- list()
for(i in 1:100){
  x <- 1/i
  y <- 1-x
  out[[i]] <- f1(x,y)
}
res <- do.call("rbind", out)

plot((res$Pa-res$Pd)*(res$Sa-res$Sd), res$score,
     xlab="(Pa - Pd) * (Sa - Sd)", ylab="Score")


####### END of STRUCTURED OUTPUT ###################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###







#########################
## UNSTRUCTURED OUTPUT ##
#########################
## descendant & ancestor = random variables in [0,1]

########
## f2 ##
########

f2 <- function(n){

  df <- list()

  for(i in 1:n){

    pa <- sample(seq(0, 1, 0.01), 1)
    pd <- sample(seq(0, 1, 0.01), 1)
    sa <- sample(seq(0, 1, 0.01), 1)
    sd <- sample(seq(0, 1, 0.01), 1)
    l <- 1

    #     score <- abs(((pa*sa + (1-pa)*(1-sa) - pa*(1-sa) - (1-pa)*sa) +
    #                     (pd*sd + (1-pd)*(1-sd) - pd*(1-sd) - (1-pd)*sd)) * l/2)

    score <- get.score3(Pa=pa, Pd=pd, Sa=sa, Sd=sd, l=l)

    df[[i]] <- data.frame(pa, pd, sa, sd, score)
    names(df[[i]]) <- c("Pa", "Pd", "Sa", "Sd", "score")

  }

  out <- do.call("rbind", df)
  return(out)

} # end f2

################################################
## Run f2 function to get unstructured output ##
################################################

res <- f2(n=5000)

plot((res$Pa-res$Pd)*(res$Sa-res$Sd), res$score,
     xlab="(Pa - Pd) * (Sa - Sd)", ylab="Score")

#######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###









#######

# score <- as.vector(unlist(score))
# Pa <- as.vector(unlist(Pa))
# Pd <- as.vector(unlist(Pd))
# Sa <- as.vector(unlist(Sa))
# Sd <- as.vector(unlist(Sd))
# PaPd <- score ~ Pa*Pd
# df <- data.frame(score, Pa, Pd, Sa, Sd)
#
# plot(Sa-Pd, score)
#
# fit <- lm(formula=score~Pa*Pd, data=df, x=TRUE)
# plot(fit)
#
# library(interplot)
# PaPd <- lm(score ~ Pa * Pd, data=df)
# # summary(PaPd)
# interplot(m=PaPd, var1="Pa", var2="Pd", xmin=0, xmax=1)
#
# SdPd <- lm(score ~ Sd * Pd, data=df)
# # summary(PaPd)
# interplot(m=SdPd, var1="Sd", var2="Pd", xmin=0, xmax=1)
#
# SdPd <- lm(score ~ Sd * Pd, data=df)
# # summary(PaPd)
# interplot(m=SdPd, var1="Sd", var2="Pd", xmin=0, xmax=1)
#
# df2 <- data.frame(score, PaPd=(Pa-Pd), SaSd=(Sa-Sd))
#
# SaSdPaPd <- lm(score ~ (SaSd) * (PaPd), data=df2)
# # summary(PaPd)
# interplot(m=SaSdPaPd, var1="SaSd", var2="PaPd", xmin=0, xmax=1)
# interplot(m=SaSdPaPd, var1="PaPd", var2="SaSd", xmin=0, xmax=1)
#
# plot((Sa-Sd)*(Pa-Pd), score)
#
# library(cluster)
# library(fpc)
#
# data(iris)
# dat <- iris[, -5] # without known classification
# # Kmeans clustre analysis
# clus <- kmeans(dat, centers=3)
#
# with(iris, pairs(dat, col=c(1:3)[clus$cluster]))







#########

## TESTING SCORE 3 ALTERNATIVES:

## Original SCORE 3:

## Should have MAX score:
get.score3(Pa=1, Sa=1, Pd=0, Sd=0, l=1) # 0.33

get.score3(Pa=0, Sa=0, Pd=1, Sd=1, l=1) # 0.33

get.score3(Pa=1, Sa=0, Pd=0, Sd=1, l=1) # -0.33

## Should have <= Max score..
get.score3(Pa=1, Sa=1, Pd=1, Sd=1, l=1) # 1

get.score3(Pa=0, Sa=0, Pd=0, Sd=0, l=1) # 1

get.score3(Pa=1, Sa=0, Pd=1, Sd=0, l=1) # -1

## Should have ZERO score:
get.score3(Pa=0.5, Sa=0.5, Pd=0.5, Sd=0.5, l=1) # 0

get.score3(Pa=0.5, Sa=0, Pd=0.5, Sd=1, l=1) # 0

get.score3(Pa=0, Sa=1, Pd=1, Sd=1, l=1) # 0 ## Doesn't this make it a NON-subsequent score???!

######

## SCORE 1 - like score 3:

## ie:
# score3 <- (l/2)*(
#   Pa*Sa + (1 - Pa)*(1 - Sa) -
#   Pa*(1 - Sa) - (1 - Pa)*Sa +
#   Pd*Sd + (1 - Pd)*(1 - Sd) -
#   Pd*(1 - Sd) - (1 - Pd)*Sd)

## Should have MAX score:
get.score3(Pa=1, Sa=1, Pd=0, Sd=0, l=1) # 1

get.score3(Pa=0, Sa=0, Pd=1, Sd=1, l=1) # 1

get.score3(Pa=1, Sa=0, Pd=0, Sd=1, l=1) # -1

## Should have <= Max score..
get.score3(Pa=1, Sa=1, Pd=1, Sd=1, l=1) # 1

get.score3(Pa=0, Sa=0, Pd=0, Sd=0, l=1) # 1

get.score3(Pa=1, Sa=0, Pd=1, Sd=0, l=1) # -1

## Should have ZERO score:
get.score3(Pa=0.5, Sa=0.5, Pd=0.5, Sd=0.5, l=1) # 0

get.score3(Pa=0.5, Sa=0, Pd=0.5, Sd=1, l=1) # 0

get.score3(Pa=1, Sa=1, Pd=0, Sd=1, l=1) # 0 (Same state at top (+1), but opposite states at bottom (-1))

## in between scores:
get.score3(Pa=1, Sa=1, Pd=0.25, Sd=0.25, l=1) # 0.625

get.score3(Pa=1, Sa=0.25, Pd=0.25, Sd=1, l=1) # -0.5

get.score3(Pa=0.75, Sa=0.25, Pd=0.25, Sd=0.75, l=1) # -0.25

get.score3(Pa=0.25, Sa=0.25, Pd=0.75, Sd=0.75, l=1) # 0.25


## SCORE 1 - like score 3 (recentred around 0):
get.score3(Pa=1, Sa=1, Pd=-1, Sd=-1, l=1) # 5

get.score3(Pa=1, Sa=-1, Pd=-1, Sd=1, l=1) # -3

get.score3(Pa=1, Sa=1, Pd=1, Sd=-1, l=1) # -1

##


## Phylo Independent Contrasts?
# pic.X <- pic(x, tree)
# pic.Y <- pic(y, tree)
# pic.out <- cor.test(pic.X, pic.Y)
## --> try cor.test fn for paired samples (ie. Pa, Sa and Pd, Sd (ie P and S for all nodes)):

COR <- list()
PIC <- PIC.X <- list()
PIC.Y <- pic(phen, tree)
for(i in 1:ncol(snps.rec)){
  # COR[[i]] <- cor.test(x = snps.rec[,i], y = phen.rec,
  #                      alternative = "two.sided") # same as reg cor + test..

  PIC.X[[i]] <- pic(snps[,i], tree)
  PIC[[i]] <- cor.test(PIC.X[[i]], PIC.Y,
                       alternative = "two.sided")


}

PIC.p <- sapply(c(1:length(PIC)), function(e) PIC[[e]]$p.value)
PIC.cor <- sapply(c(1:length(PIC)), function(e) PIC[[e]]$estimate)
