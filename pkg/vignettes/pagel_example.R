

################
## Pagel 1994 ##
################


#######################
## Example w my data ##
#######################

library(devtools)
install_github("liamrevell/phytools")
require(phytools)

## Old phytools (??)
## Version: 0.5-25

## New phytools
## Github dl April 7, 2016

data("snps.ace")
data("phen.ace")
data("tree.ace")
snps <- snps.ace
phen <- phen.ace
tree <- tree.ace

## Attempt w SNP1 ##
var <- snps[,5]

x <- var
y <- phen

names(x) <- names(phen)

fit.ape.eq <- fitPagel(tree,x,y, equal=TRUE, method="ace") # suppressWarnings(fitPagel(tree,x,y, method="ace")) # warnings OK (?? if perfectly correlated only??)
Q <- fit.ape.eq$independent.Q
Q
get.n.subs(Q, tree)

## MODIFY fitPagel to return ace output,
## at least for independent model:
lik.anc <- fit.iQ$lik.anc
lik.anc <- lik.anc*100
lik.anc <- data.frame((lik.anc[,1] + lik.anc[,2]), (lik.anc[,3] + lik.anc[,4]))
names(lik.anc) <- c("0", "1")
#lik.anc <- replace(lik.anc, which(lik.anc < 20), 0)
head(round(lik.anc,4))
round(lik.anc,4)


## RUN ace ##

## on x (snp)
ace.snp <- ace(x, tree, type="discrete")
## NB: need to remove "call" from list to test id...
identical(ace.snp[-(6)], ace.eq.snp[-(6)]) # TRUE
identical(ace.snp[-(6)], ace.diff.snp[-(6)]) # FALSE

ace.eq.snp <- ace(x, tree, type="discrete", model= matrix(c(0, 1, 1, 0), 2))
str(ace.eq.snp)
ace.eq.snp$rates # 1 rate

ace.diff.snp <- ace(x, tree, type="discrete", model= matrix(c(0, 1, 2, 0), 2))
str(ace.diff.snp)
ace.diff.snp$rates # 2 rates

## on y (phen) ##
ace.eq.phen <- ace(phen, tree, type="discrete", model= matrix(c(0, 1, 1, 0), 2))
str(ace.eq.phen)
ace.eq.phen$rates

ace.diff.phen <- ace(phen, tree, type="discrete", model= matrix(c(0, 1, 2, 0), 2))
str(ace.diff.phen)
ace.diff.phen$rates

## get probs from ACE rates (????)
exp(ace.eq.snp$rates)
exp(ace.diff.snp$rates)
exp(ace.eq.phen$rates)
exp(ace.diff.phen$rates)

## get probs from fitPagel rate matrices:
Q.fitPagel.i <- fit.ape1$independent.Q
Q.fitPagel.d <- fit.ape1$dependent.Q

P.fitPagel.i <- Q2P(Q=Q.fitPagel.i, type="independent")
P.fitPagel.d <- Q2P(Q=Q.fitPagel.d, type="dependent")

P.fitPagel.i$probs
P.fitPagel.d$probs


## RUN fitPagel ##
#fit.ape1 <- suppressWarnings(fitPagel(tree,x,y, method="ace")) # warnings OK (?? if perfectly correlated only??)
## warnings = In nlminb... imaginary parts discarded in coercion
#save(fit.ape1, file="~/treeWAS/misc/fit.ape1_old.phytools.Rdata")
fit.ape1 <- get(load("~/treeWAS/misc/fit.ape1_old.phytools.Rdata"))
str(fit.ape1)

fit.ape2 <- fitPagel(tree,x,y, equal=FALSE, method="ace") # suppressWarnings(fitPagel(tree,x,y, method="ace")) # warnings OK (?? if perfectly correlated only??)
## warnings = In nlminb... imaginary parts discarded in coercion
str(fit.ape2)
Q <- fit.ape2$independent.Q
get.n.subs(Q, tree)

fit.ape.eq <- fitPagel(tree,x,y, equal=TRUE, method="ace") # suppressWarnings(fitPagel(tree,x,y, method="ace")) # warnings OK (?? if perfectly correlated only??)
## warnings = In nlminb... imaginary parts discarded in coercion
# str(fit.ape.eq)
Q <- fit.ape.eq$independent.Q
get.n.subs(Q, tree)

fit.ape.eq2 <- fitPagel(tree,x,y, equal=TRUE, method="ace") # suppressWarnings(fitPagel(tree,x,y, method="ace")) # warnings OK (?? if perfectly correlated only??)
## warnings = In nlminb... imaginary parts discarded in coercion
str(fit.ape.eq2)
Q <- fit.ape.eq2$independent.Q
get.n.subs(Q, tree)
# identical(fit.ape1, fit.ape2) # TRUE


## RUN fitMk ? ##
#fit.Mk1 <- fitPagel(tree,x,y, method="fitMk")
# Error in fitMk(tree, xy, model = iQ, ...) :
#   model does not have the right number of columns
fit.Mk2 <- fitPagel(tree,x,y, method="fitMk")
## No error!
str(fit.Mk2)
## p-val = 7.49e-05
## BUT p-val = 1 for snps[,2] :(


## RUN fitDiscrete ?
## MUCH SLOWER THAN APE/ACE
#fit.disc1 <- fitPagel(tree,x,y, method="fitDiscrete")
#save(fit.disc1, file="~/treeWAS/misc/fit.disc1_old.phytools.Rdata")
fit.disc1 <- get(load("~/treeWAS/misc/fit.disc1_old.phytools.Rdata"))
str(fit.disc1)
# p-val == 1

fit.disc2 <- fitPagel(tree,x,y, method="fitDiscrete")
str(fit.disc2)

identical(fit.disc1, fit.disc2)
## FALSE!
## BUT p-val still = 1 :(
## Everything looks identical, except lik.ratio
## which was 0 in fit.disc1 (old version of phytools)
## and is -2.02e-10 in fit.disc2 (new version)...
###############

mat <- seq(0,1,1/10000)[-10001]
str(mat)
res=0
for (i in mat)
  res=res+A(i)
res=res/10000

fn <- function(Q){
  return(as.vector(Q))
}

integrate(fn, 0, 1)

## From Q2P fn:
P <- mexp(0.05 * Q)
row.names(P) <- row.names(Q)
colnames(P) <- colnames(Q)
P

## CHECK:
## Make sure pairs of cells add up to same probabilities.
if(any(c(sum(P[1,1:2]) != sum(P[2,1:2]),
         sum(P[1,3:4]) != sum(P[2,3:4]),
         sum(P[3,1:2]) != sum(P[4,1:2]),
         sum(P[3,3:4]) != sum(P[4,3:4])))) stop("Oops--something went wrong!
                                                Pairs of cells in the probability matrix
                                                do not sum to the same thing.")

P00 <- sum(P[1, 1:2])
P01 <- sum(P[1, 3:4])
P11 <- sum(P[3, 1:2])
P10 <- sum(P[3, 3:4])

#####

p.val <- fit.ape$P
p.val

## Get INDEPENDENT rates
rate.mat <- fit.ape1$independent.Q
rate.mat

Q <- rate.mat

probs.fitPagel <- list(Q, P, P00, P01, P10, P11)

probs.ace <- list(Q, P, P00, P01, P10, P11)


################################################################################################################################################################

##########
## PLOT ##
##########

probs.d <- lik.anc[,2]

## convert probs to show on tree pies:
probs.ori <- probs.d
df <- data.frame(100-probs.ori, probs.ori)
names(df) <- c("prob.A", "prob.B")
head(df)
# save(df, file="./MRC_probs_ACE.Rdata")

lab <- as.character(tree$tip.label)
temp <- sapply(c(1:length(lab)),
               function(e)
                 paste(
                   paste(rep(c("-"), 3-nchar(lab[e])), collapse=""),
                   "-", lab[e], sep=""))
tiptext <- as.vector(unlist(temp))

par(mar=c(1,1,1,2)-0.5)
edgeCol <- "black"
plot(tree, show.tip=FALSE, edge.width=2,
     edge.color=edgeCol, type="c", use.edge.length=FALSE)
## PHEN-coloured tip labels:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(phen)))
tiplabels(text=tiptext, cex=0.6, adj=c(-.4, 0), col=myCol[tipCol.pattern], frame="none")


## try plotting ~ example: ##
#### Showing the probs and inverse probs on each node:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(var[1:(tree$Nnode+1)])))
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
nodelabels(thermo = df/100, piecol = transp(myCol, 0.8), cex = 0.25)
nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp("blue", 0.3), adj=2)


##########################################################################################################################################




rates <- as.vector(unlist(t(rate.mat)))[c(3,9)]
rate01 <- rates[1]
rate10 <- rates[2]

rate01
rate10

exp(rate01)
exp(rate10)

sum(rate01, rate10)

sum(exp(rate01), exp(rate10))

exp(sum(rate01, rate10))

## Get DEPENDENT rates
rate.mat.d <- fit.ape$dependent.Q
rate.mat.d

## unbind the matrix by rows:
rates <- as.vector(t(rate.mat))
## get (4) unique values (for independent rate matrices)
rates <- unique(rates[rates>0])

rates[2]
rates[4]
rates[2]+rates[4]

## Note (from Pagel 94): ##
## If the 2 variables change INDEPENDENTLY of each other,
## then their JOINT Prob's of change on any given branch
## is equal to the PRODUCT of their SEPARATE prob's of change.
ACE <- ace(x, tree, type="discrete")
str(ACE)

ACE.eq <- ace(snps[,1], tree, type="discrete", model= matrix(c(0, 1, 1, 0), 2))
# str(ACE.eq)
ACE.eq$rates

ACE.diff <- ace(snps[,1], tree, type="discrete", model= matrix(c(0, 1, 2, 0), 2))
# str(ACE.diff)
ACE.diff$rates


## P(t) = exp[Qt] ##

r <- rates

exp.mat <- expm(rate.mat)

## to get the

r[3]/r[4]

r[2]/r[1]
#


#######################

fit.ape <- fitPagel(tree,x,y, method="ace") # warnings OK (?? if perfectly correlated only??)
fit.ape

fit.Mk<-fitPagel(tree,x,y, method="fitMk")
fit.Mk

## Not sure why this one is giving me p-values of 1 and disagreeing w ape's ace...??
fit.geiger<-fitPagel(tree,x,y,method="fitDiscrete") # warnings OK
fit.geiger

###############

n.ind <- 100

set.seed(1)
tree <- pbtree(n=n.ind, scale=1)


###############   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

#################
## DIVERSITREE ##
#################
require(diversitree)
require(adegenet) # transp

## run make.mk2() or make.mkn()
## then run asr.marginal()

## make.mkn.multitrait ##

var <- snps[,5]

x <- var
y <- phen

#if(method=="fitMkn")

## fit independent model (4 rates)
iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
rownames(iQ)<-colnames(iQ)<-levels(xy)

tr <- mkn.multitrait.translate(2)
#tr

pars <- c(.05, .05, .0,       # q00.10, q00.01, q00.11
          .05, .0,  .05,      # q10.00, q10.01, q10.11
          .05, .0,  .05,      # q01.00, q01.10, q01.11
          .0,  .05, .05)      # q11.00, q11.10, q11.01)

states <- expand.grid(A=0:1, B=0:1)[phy$tip.state,]
rownames(states) <- phy$tip.label

states <- data.frame(x, y)
rownames(states) <- tree$tip.label

lik <- make.mkn.multitrait(tree, states, strict=FALSE)
# pars <- Q[c(1:length(Q))] ## pars should contain all 8 rates for independent & 4 rates for dependent...
pars <- c(1:8)
fit <- find.mle(lik, pars, method="subplex") ## CHECK method ## CHECK if this should be anova/amova instead??
st <- asr.marginal(lik, coef(fit))

## fit dependent model (8 states)
dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
rownames(dQ)<-colnames(dQ)<-levels(xy)

#######################################################
## EXAMPLE ##
#############

## The translation between these two bases is fairly straightforward; if
## we have a vector of parameters in our new basis 'p' we can convert it
## into the original MuSSE basis ('q') through this matrix:
tr <- musse.multitrait.translate(2)
tr

## Notice that the rows that correspond to transitions in multiple
## traits are all zero by default; this means that these q values will
## be zero regardless of the parameter vector used.
tr["q00.11",]

## And here is the section of the transition matrix corresponding to the
## lambda values; every rate gets a contribution from the intercept term
## (lambda0), lambda10 and lambda11 get a contribution from lambdaA, etc.
tr[1:4,1:4]

## There is currently no nice simulation support for this, so bear with
## an ugly script to generate the tree and traits.
pars <- c(.10, .15, .20, .25, # lambda 00, 10, 01, 11
          .03, .03, .03, .03, # mu 00, 10, 01, 11
          .05, .05, .0,       # q00.10, q00.01, q00.11
          .05, .0,  .05,      # q10.00, q10.01, q10.11
          .05, .0,  .05,      # q01.00, q01.10, q01.11
          .0,  .05, .05)      # q11.00, q11.10, q11.01
set.seed(2)
phy <- tree.musse(pars, 60, x0=1)

states <- expand.grid(A=0:1, B=0:1)[phy$tip.state,]
rownames(states) <- phy$tip.label


####

####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####

############################
## make.mkn with my data! ##
############################

states <- as.character(xy)
states <- replace(states, which(states=="0|0"), 1)
states <- replace(states, which(states=="0|1"), 2)
states <- replace(states, which(states=="1|0"), 3)
states <- replace(states, which(states=="1|1"), 4)
states <- as.numeric(states)
names(states) <- names(xy)

# k <- 4

## fit independent model (4 rates)
# iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
# rownames(iQ)<-colnames(iQ)<-  1:k # levels(xy)

## fit independent model (4 rates)
# iQ<-matrix(c(0,1,0,2,0,3,4,0,1,0,4,0),4,3,byrow=TRUE)
# rownames(iQ)<-colnames(iQ)<- levels(as.factor(states))

## dealing w empty levels:
k <- length(levels(as.factor(states)))
if(max(as.numeric(states)) > k){
  states <- replace(states, which(states==max(as.numeric(states))), k)
}
iQ <- matrix(rep(0, k*k), k, k, byrow=TRUE)
iQ[1,2] <- iQ[2,1] <- 1
iQ[2,3] <- iQ[3,2] <- 2

pars <- c(iQ[1,c(2,3)],
          iQ[2,c(1,3)],
          iQ[3,c(1,2)])

lik.mkn <- make.mkn(tree, states, k=k, strict=FALSE)
fit.mkn.i <- find.mle(lik.mkn, pars, method="subplex")
st.i <- asr.marginal(lik.mkn, coef(fit.mkn.i))


## fit dependent model (8 states)
# dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
# rownames(dQ)<-colnames(dQ)<-levels(xy)
dQ <- matrix(rep(0, k*k), k, k, byrow=TRUE)
dQ[1,2] <- 1
dQ[2,1] <- 2
dQ[2,3] <- 3
dQ[3,2] <- 4

pars <- c(dQ[1,c(2,3)],
          dQ[2,c(1,3)],
          dQ[3,c(1,2)])

#lik.mkn <- make.mkn(tree, states, k=k, strict=FALSE)
fit.mkn.d <- find.mle(lik.mkn, pars, method="subplex")
st.d <- asr.marginal(lik.mkn, coef(fit.mkn.d))

all.equal(st.i, st.d)
# "Mean relative difference: 3.46061e-06"

## SEE EXAMPLE IN ?asr.mkn for where to go from here
## wrt testing sig diff btw independent and dependent outputs of asr.marginal generated above

##

####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####




# tr <- mkn.multitrait.translate(2)
# tr
#
# pars <- c(.05, .05, .0,       # q00.10, q00.01, q00.11
#           .05, .0,  .05,      # q10.00, q10.01, q10.11
#           .05, .0,  .05,      # q01.00, q01.10, q01.11
#           .0,  .05, .05)      # q11.00, q11.10, q11.01)
#
# states <- expand.grid(A=0:1, B=0:1)[phy$tip.state,]
# rownames(states) <- phy$tip.label
#
# states <- data.frame(x, y)
# rownames(states) <- tree$tip.label
#
# lik <- make.mkn.multitrait(tree, states, strict=FALSE)
# # pars <- Q[c(1:length(Q))] ## pars should contain all 8 rates for independent & 4 rates for dependent...
# pars <- c(1:8)
# fit <- find.mle(lik, pars, method="subplex") ## CHECK method ## CHECK if this should be anova/amova instead??
# st <- asr.marginal(lik, coef(fit))
#######################
## make.mkn example: ##
#######################
## Simulate a tree and character distribution.  This is on a birth-death
## tree, with high rates of character evolution and an asymmetry in the
## character transition rates.
pars <- c(.1, .1, .03, .03, .1, .2)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa=25, max.t=Inf, x0=0)[[1]]

## Here is the 25 species tree with the true character history coded.
## Red is state '1', which has twice the character transition rate of
## black (state '0').
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy)

## Maximum likelihood parameter estimation:
p <- c(.1, .1) # initial parameter guess

lik <- make.mk2(phy, phy$tip.state)
fit.mk2 <- find.mle(lik, p)
coef(fit.mk2)   # q10 >> q01
logLik(fit.mk2) # -10.9057

## This can also be done using the more general Mk-n.
## This uses an approximation for the likelihood calculations.  make.mkn
## assumes that states are numbered 1, 2, ..., k, so 1 needs to be added
## to the states returned by trees.
lik.mkn <- make.mkn(phy, phy$tip.state + 1, 2)
fit.mkn <- find.mle(lik.mkn, p)
fit.mkn[1:2]

## These are the same (except for the naming of arguments)
all.equal(fit.mkn[-7], fit.mk2[-7], check.attr=FALSE, tolerance=1e-7)

## Equivalence to ape's \link{ace} function:
model <- matrix(c(0, 2, 1, 0), 2)
fit.ape <- ace(phy$tip.state, phy, "discrete", model=model, ip=p)

## To do the comparison, we need to rerun the diversitree version with
## the same root conditions as ape.
fit.mk2 <- find.mle(lik, p, root=ROOT.GIVEN, root.p=c(1,1))

## These are the same to a reasonable degree of accuracy, too (the
## matrix exponentiation is slightly less accurate than the ODE
## solving approach.  The make.mk2 version is exact)
all.equal(fit.ape[c("rates", "loglik")], fit.mk2[1:2],
          check.attributes=FALSE, tolerance=1e-4)

## The ODE calculation method may be useful when there are a large
## number of possible states (say, over 20).
lik.ode <- make.mkn(phy, phy$tip.state + 1, 2,
                    control=list(method="ode"))
fit.ode <- find.mle(lik.ode, p)
fit.ode[1:2]

all.equal(fit.ode[-7], fit.mkn[-7], tolerance=1e-7)

####

foo <- make.mk2(tree) ####

lik.m <- make.mk2(tree, var)
fit.m <- find.mle(lik.m, pars[5:6], method="subplex")
st.m <- asr.marginal(lik.m, coef(fit.m))

#######################################
## comparing bisse and mk2 w my data ##
#######################################
## BiSSE ancestral state reconstructions under the ML model
# lik <- make.bisse(tree, var)
# fit <- find.mle(lik, pars, method="subplex")
# st <- asr.marginal(lik, coef(fit))
system.time(st <- asr.marginal(make.bisse(tree, var), coef(find.mle(lik, pars, method="subplex"))))
myCol <- c("blue", "red")
nodelabels(thermo=t(st), piecol=transp(myCol, 0.8), cex=.25, adj=-0.5)

## Mk2 ancestral state reconstructions, ignoring the shifts in
## diversification rates:
# lik.m <- make.mk2(tree, var)
# fit.m <- find.mle(lik.m, pars[5:6], method="subplex")
# st.m <- asr.marginal(lik.m, coef(fit.m))
system.time(st.m <- asr.marginal(make.mk2(tree, var), coef(find.mle(lik.m, pars[5:6], method="subplex"))))
## The Mk2 results have more uncertainty at the root, but both are
## similar.
nodelabels(thermo=t(st.m), piecol=transp(myCol, 0.8), cex=.25, adj=-1)

#############################
## ALTERNATIVES TO ace FN: ##
#############################
## EG: from diversitree: bisse, mk2...

## Start with a simple tree evolved under a BiSSE with all rates
## asymmetric:
pars <- c(.1, .2, .03, .06, .01, .02)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa=50, max.t=Inf, x0=0)[[1]]

## Here is the true history
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy, main="True history")

## BiSSE ancestral state reconstructions under the ML model
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars, method="subplex")
st <- asr.marginal(lik, coef(fit))
nodelabels(thermo=t(st), piecol=1:2, cex=.5)

## Mk2 ancestral state reconstructions, ignoring the shifts in
## diversification rates:
lik.m <- make.mk2(phy, phy$tip.state)
fit.m <- find.mle(lik.m, pars[5:6], method="subplex")
st.m <- asr.marginal(lik.m, coef(fit.m))
## The Mk2 results have more uncertainty at the root, but both are
## similar.
nodelabels(thermo=t(st.m), piecol=1:2, cex=.5, adj=-.5)

#####
## section Not run .. ##
#####

## Equivalency of Mk2 and BiSSE where diversification is state
## independent. For any values of lambda/mu (here .1 and .03) where
## these do not vary across character states, these two methods will
## give essentially identical marginal ancestral state reconstructions.
st.id <- asr.marginal(lik, c(.1, .1, .03, .03, coef(fit.m)))
st.id.m <- asr.marginal(lik.m, coef(fit.m))
## Reconstructions are identical to a relative tolerance of 1e-7
## (0.0000001), which is similar to the expected tolerance of the BiSSE
## calculations.
all.equal(st.id, st.id.m, tolerance=1e-7)
## Equivalency of BiSSE and MuSSE reconstructions for two states:
lik.b <- make.bisse(phy, phy$tip.state)
lik.m <- make.musse(phy, phy$tip.state + 1, 2)
st.b <- asr.marginal(lik.b, coef(fit))
st.m <- asr.marginal(lik.m, coef(fit))
all.equal(st.b, st.m)

###############   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
################
## get.n.subs ##
################

## identify the n.subs generated by sim.history
## if the transition rate matrix Q is -1, 1, 1, -1 divided by n.

## HYPOTHESIS:
## Setting n to sum(tree$edge.length) should make
## the expected value of n.subs to be 1 sub per tree...

## EG:
# foo <- replicate(100, get.n.subs(sum(tree$edge.length)))
# hist(foo)
# summary(foo)




## plot
#library(adegenet)
#par(mfrow=c(1,1))
plotSimmap(tt1,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1)
edgelabels(c(1:nrow(tree$edge)), cex=0.5, bg=transp("yellow", 0.3))


tt2<-sim.history(tree,Q)
## Done simulation(s).

## these are uncorrelated, see:
par(mfrow=c(1,2))
plotSimmap(tt1,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1)
plotSimmap(tt2,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1,direction="leftwards")



## print method for objects of class "fitMk"
print.fitMk<-function(x,digits=6,...){
  cat("Object of class \"fitMk\".\n\n")
  cat("Fitted (or set) value of Q:\n")
  Q<-matrix(NA,length(x$states),length(x$states))
  Q[]<-c(0,x$rates)[x$index.matrix+1]
  diag(Q)<-0
  diag(Q)<--rowSums(Q)
  colnames(Q)<-rownames(Q)<-x$states
  print(round(Q,digits))
  cat("\nFitted (or set) value of pi:\n")
  print(x$pi)
  cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n\n"))
}

## summary method for objects of class "fitMk"
summary.fitMk<-function(object,...){
  if(hasArg(digits)) digits<-list(...)$digits
  else digits<-6
  if(hasArg(quiet)) quiet<-list(...)$quiet
  else quiet<-FALSE
  if(!quiet) cat("Fitted (or set) value of Q:\n")
  Q<-matrix(NA,length(object$states),length(object$states))
  Q[]<-c(0,object$rates)[object$index.matrix+1]
  diag(Q)<-0
  diag(Q)<--rowSums(Q)
  colnames(Q)<-rownames(Q)<-object$states
  if(!quiet) print(round(Q,digits))
  if(!quiet) cat(paste("\nLog-likelihood:",round(object$logLik,digits),"\n\n"))
  invisible(list(Q=Q,logLik=object$logLik))
}

## logLik method for objects of class "fitMk"
logLik.fitMk<-function(object,...) object$logLik

## AIC method
AIC.fitMk<-function(object,...,k=2){
  np<-length(object$rates)
  -2*logLik(object)+np*k
}


fitPagel<-function(tree,x,y,method="fitMk",...){
  if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
  if(method=="fitDiscrete"){
    chk<-.check.pkg("geiger")
    if(!chk){
      cat("  method = \"fitDiscrete\" requires the package \"geiger\"\n")
      cat("  Defaulting to method = \"fitMk\"\n\n")
      method<-"fitMk"
      fitDiscrete<-function(...) NULL
    }
  }
  if(!is.factor(x)) x<-as.factor(x)
  levels.x<-levels(x)
  if(!is.factor(y)) y<-as.factor(y)
  levels.y<-levels(y)
  y<-y[names(x)]
  if(length(levels.x)!=2||length(levels.y)!=2)
    stop("Only binary characters for x & y currently permitted.")
  xy<-setNames(factor(paste(x,y,sep="|"),
                      levels=sapply(levels.x,paste,levels.y,sep="|")),
               names(x))
  ## fit independent model
  iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
  rownames(iQ)<-colnames(iQ)<-levels(xy)
  fit.iQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=iQ,...)
  else if(method=="ace") ace(xy,tree,type="discrete",model=iQ,...)
  else if(method=="fitMk") fitMk(tree,xy,model=iQ,...)
  #else if(method=="fitMkn")
  ## fit dependendent model
  dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
  rownames(dQ)<-colnames(dQ)<-levels(xy)
  fit.dQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=dQ,...)
  else if(method=="ace") ace(xy,tree,type="discrete",model=dQ,...)
  else if(method=="fitMk") fitMk(tree,xy,model=dQ,...)
  #else if(method=="fitMkn")
  ## back translate independent model
  if(method=="fitDiscrete") iQ<-.Qmatrix.from.gfit(fit.iQ)
  else {
    I<-fit.iQ$index.matrix
    I[I==0]<-NA
    iQ<-apply(I,2,function(i,x) x[i],x=fit.iQ$rates)
    iQ[is.na(iQ)]<-0
    diag(iQ)<--rowSums(iQ)
    rownames(iQ)<-colnames(iQ)
  }
  ## dependent model
  if(method=="fitDiscrete") dQ<-.Qmatrix.from.gfit(fit.dQ)
  else {
    I<-fit.dQ$index.matrix
    I[I==0]<-NA
    dQ<-apply(I,2,function(i,x) x[i],x=fit.dQ$rates)
    dQ[is.na(dQ)]<-0
    diag(dQ)<--rowSums(dQ)
    rownames(dQ)<-colnames(dQ)
  }
  ## assemble object to return
  obj<-list(independent.Q=iQ,
            dependent.Q=dQ,
            independent.logL=logLik(fit.iQ),
            dependent.logL=logLik(fit.dQ),
            lik.ratio=2*(logLik(fit.dQ)-logLik(fit.iQ)),
            P=pchisq(2*(logLik(fit.dQ)-logLik(fit.iQ)),
                     df=length(levels(x))+length(levels(y)),
                     lower.tail=FALSE),
            method=method)
  class(obj)<-"fitPagel"
  obj
}

## print method for objects of class "fitPagel"
## written by Liam J. Revell 2014

print.fitPagel<-function(x,...){
  cat("\n  Pagel's binary character correlation test:\n")
  cat("\nIndependent model rate matrix:\n")
  print(x$independent.Q)
  cat("\nDependent model rate matrix:\n")
  print(x$dependent.Q)
  cat("\nModel fit:\n")
  obj<-matrix(c(x$independent.logL,x$dependent.logL),2,1)
  rownames(obj)<-c("independent","dependent")
  colnames(obj)<-"log-likelihood"
  print(obj)
  cat("\nHypothesis test result:\n")
  cat(paste("  likelihood-ratio: ",signif(x$lik.ratio,7),"\n"))
  cat(paste("  p-value: ",signif(x$P,7),"\n"))
  cat(paste("\nModel fitting method used was",x$method,"\n\n"))
}

## function borrowed from geiger to pull the Q-matrix from a fit returned by fitDiscrete

.Qmatrix.from.gfit<-function(x){
  if(!.check.pkg("geiger")) argn<-function(...) NULL
  lik=x$lik
  numberize=function(x){
    y=gsub("q","",x)
    sp=(nn<-nchar(y))/2
    as.numeric(c(substring(y,1,sp),substring(y,sp+1,
                                             nn)))
  }
  att=attributes(lik)
  att$k=length(att$levels)
  Qmat=matrix(0,att$k,att$k)
  nms=att$argn[att$trns]
  other=att$argn[!att$trns]
  if("constrained"%in%class(lik)){
    cpars=x$opt[argn(lik)]
    apars=names(lik(unlist(cpars),pars.only=TRUE))
    nms=apars[!apars%in%other]
  }
  trns=x$opt[nms]
  for(i in 1:length(trns)){
    nm=names(trns)[i]
    idx=numberize(nm)
    Qmat[idx[1],idx[2]]=trns[[i]]
  }
  diag(Qmat)=-rowSums(Qmat)
  rownames(Qmat)<-colnames(Qmat)<-levels(lik)
  Qmat
}

## BELOW NOW IN PHYTOOLS (?!) ######################################################################



##################################################################################
## EXAMPLE FROM: #################################################################
##################################################################################
## http://blog.phytools.org/2014/12/r-function-for-pagels-1994-correlation.html ##
##################################################################################

## first load packages & source code
library(phytools)
library(geiger)
# source("fitPagel.R")
.check.pkg<-phytools:::.check.pkg

##################
## UNcorrelated ##
##################

set.seed(1)
## now let's simulate some uncorrelated data
tree<-pbtree(n=300,scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-letters[1:2]
tt1<-sim.history(tree,Q)
## Done simulation(s).

tt2<-sim.history(tree,Q)
## Done simulation(s).

## these are uncorrelated, see:
par(mfrow=c(1,2))
plotSimmap(tt1,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1)
plotSimmap(tt2,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1,direction="leftwards")


## run Pagel's binary character correlation test:
x<-tt1$states
y<-tt2$states

fit.ape<-fitPagel(tree,x,y, method="ace")
fit.ape

fit.Mk <- fitPagel(tree, x, y, method="fitMk")
fit.Mk

fit.geiger<-fitPagel(tree,x,y,method="fitDiscrete") # warnings OK
fit.geiger


################
## CORRELATED ##
################

Q<-matrix(c(0,0.5,0.5,0,2,0,0,2,2,0,0,2,0,0.5,0.5,0),4,4,byrow=TRUE)
rownames(Q)<-colnames(Q)<-c("aa","ab","ba","bb")
diag(Q)<--rowSums(Q)
tt<-sim.history(tree,t(Q))

tt1<-mergeMappedStates(tt,c("aa","ab"),"a")
tt1<-mergeMappedStates(tt1,c("ba","bb"),"b")
tt2<-mergeMappedStates(tt,c("aa","ba"),"a")
tt2<-mergeMappedStates(tt2,c("ab","bb"),"b")

## these data are correlated, see:
par(mfrow=c(1,2))
plotSimmap(tt1,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1)
plotSimmap(tt2,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1,direction="leftwards")

x<-getStates(tt1,"tips")
y<-getStates(tt2,"tips")

fit.ape<-fitPagel(tree,x,y, method="ace")
fit.ape

fit.Mk <- fitPagel(tree, x, y, method="fitMk")
fit.Mk

fit.geiger<-fitPagel(tree,x,y,method="fitDiscrete") # warnings OK
fit.geiger












# set.seed(4)
# tree <- rcoal(10)
# plot(tree)
# snp <- c(0,0,1,1,0,0,0,1,1,0)
# x <- y <- var <- snp
# names(x) <- names(y) <- tree$tip.label


# ######
# library(ape)
#
# x <- paste("AJ5345", 26:49, sep = "")
# x <- c("Z73494", x)
# sylvia.seq <- read.GenBank(x)
#
# sylvia.clus <- clustal(sylvia.seq)
#
# library(devtools)
# install_github()
# install_github("fmichonneau/phyloch")
# #https://github.com/fmichonneau/phyloch.git
#
# library(phyloch)
# sylvia.maff <- mafft(sylvia.seq)
# identical(sylvia.clus[x, ], sylvia.maff[x, ])
#
# taxa.sylvia <- attr(sylvia.seq, "species")
# names(taxa.sylvia) <- names(sylvia.seq)
# rm(sylvia.seq)
# taxa.sylvia[1] <- "Sylvia_atricapilla"
# taxa.sylvia[24] <- "Sylvia_abyssinica"
#
# sylvia.eco <- read.table("sylvia_data.txt")
# str(sylvia.eco)
# rownames(sylvia.eco)
# save(sylvia.clu, taxa.sylvia, sylvia.eco,
#      file = "sylvia.RData")
#
#
# sylvia <- read.table(file="~/treeWAS/misc/sylvia_data.txt", header=TRUE)
# head(sylvia)
# table(DF$geo.range, DF$mig.behav)
# tr <- rcoal(26, tip.label=row.names(DF)) ## can't find real tree...
# plot(tr)
#
# syl.er <- ace(DF$geo.range, tr, type = "d")
# str(syl.er)
#
# syl.sym <- ace(DF$geo.range, tr, type = "d", model = "SYM")#
# str(syl.sim)
#
# 1 - pchisq(2*(syl.sym$loglik - syl.er$loglik), 2)
#
# ######
