
#############################
## Running PLINK from R ?? ##
#############################

## EXAMPLE: PLINK TUTORIAL (http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml)

?shell

## Set working directory to plink program 
setwd("C:/Program Files/plink-1.07-dos")

## Confirm you are in that directory
shell("cd") 

## At this point, the command "plink" should be recognized (though you may get a different error)
shell("plink")

## open the plink example file ##
shell("plink --file C:\\PLINK\\hapmap1")

## make a binary PED file ##
## (provide the full path, not just the file name)
shell("plink --file C:\\PLINK\\hapmap1 --make-bed --out C:\\PLINK\\hapmap1")

## use --bfile to work with the BINARY file 
# (same as --file, but loads the binary one and prints summary stats)
shell("plink --bfile C:\\PLINK\\hapmap1")



################################
## BASIC ASSOCIATION ANALYSIS ##
################################

## --assoc --> 1df chi-square allelic test
## --fisher --> fisher's exact test


## perform a basic association analysis on the disease trait for all single SNPs
shell("plink --bfile C:\\PLINK\\hapmap1 --assoc --out C:\\PLINK\\as1")

## to view the file you created, just read it in with R:
as1 <- read.table("C:/PLINK/as1.assoc", header=TRUE)
head(as1)

## to get a SORTED list of association results
## that also includes significance values that are 
## ADJUSTED FOR MULTIPLE TESTING, use --adjust:
## --> TWO output files (as2.assoc & as2.assoc.adjusted)
shell("plink --bfile C:\\PLINK\\hapmap1 --assoc --adjust --out C:\\PLINK\\as2")


## again, to view the file you created, just read it in with R:
as2 <- read.table("C:/PLINK/as2.assoc", header=TRUE)
head(as2)
## and the adjusted version:
as2.adjust <- read.table("C:/PLINK/as2.assoc.adjusted", header=TRUE)
head(as2.adjust)

## Check the LOG FILE, which is where --adjust records the:
## (1) Genomic inflation factor (for Genomic Control analysis)
## (2) Mean Chi-squared statistic (should be 1, under the null)

## use "more" to look at the most significant associations
## NOTE: This prints a fucktardedly LONG table!! 
#shell("more C:\\PLINK\\as2.assoc.adjusted")


## Given that we already know about the Chinese/Japanese sub-populations, 
## we can look at the INFLATION FACTOR that results from using POPULATION MEMBERSHIP
## as the "PHENOTYPE" in a CASE/CONTROL analysis. 
## --> Replace the disease phenotype with pop.phe...
shell("plink --bfile C:\\PLINK\\hapmap1 --pheno C:\\PLINK\\pop.phe --assoc --adjust --out C:\\PLINK\\as3")

## --> OBS: Some departure from null dist (based on genomic inflation factor, chi-squared stat)
## ie. GI = 1.7
## THEREFORE: 1.7 is the max possible inflation factor that could occur IF the disease were perfectly correlated
## with the Chinese/Japanese pop split. 
## NOTE: Additional (unaccounted for) within-pop structure may also increase SNP-disease FPR. 


############################################
## GENOTYPIC and OTHER ASSOCIATION MODELS ##
############################################

## ALTERNATIVE OPTIONS:
## --> using --model (instead of --assoc)
## 1) Calculate association stats based on the 2x3 genotype table (instead of standard allelic test)
## 2) Calculate tests that assuem DOMINANT or RECESSIVE action of the minor allele
## 3) Perform the Cochran-Armitage trend test (instead of the baseic allelic test)

## Eg. To run genotypic tests for a particular SNP:
shell("plink --bfile C:\\PLINK\\hapmap1 --model --cell 0 --snp rs2222162 --out C:\\PLINK\\mod2")
## NOTE: Adding "--cell 0" allows you to over-ride the default requirement for every cell in the 2X3 table
## to have 5 observations, which does not hold here... 
## "--cell 0" changes that requirement to set the min number in each cell to 0. 

######################################################################################################

#############################
## STRATIFICATION ANALYSIS ##
#############################

## DEALING WITH POPULATION STRUCTURE...
## Use whole-genome data to cluster indivudals into homogeneous groups.
## Many options/ ways of doing this in PLINK

## Eg. IBS CLUSTERING
## (ie. CLUSTER ANALYSIS that pairs up individuals on the basis of GENETIC IDENTITY.)
## NOTE: May take a COUPLE OF MINUTES TO RUN!!

shell("plink --bfile C:\\PLINK\\hapmap1 --cluster --mc 2 --ppc 0.05 --out C:\\PLINK\\str1")
shell("more C:\\PLINK\\str1.cluster1")


###################################################
## ASSOCIATION ANALYSIS, ACCOUNTING FOR CLUSTERS ##
###################################################

## Having performed the IBS clustering/matching, we can now perform the association test, conditional on the matching. 
## Here, we use the file str1.cluster2, which contains the same info as str1.cluster1 but
## in the format of a cluster variable file that can be used with the --within option. 

##################
## CLUSTERING 1 ##
##################

## The Cochran-Mantel-Haenszel (CMH) association stat is used here
## The CMH stat tests for SNP-disease association CONDITIONAL on the clustering supplied by the cluster file. 
## --> effectively reduces the GC inflation factor to ~ 1.00. 
## The --adjust option allows us to get a SORTED list of CMH association results. 

shell("plink --bfile C:\\PLINK\\hapmap1 --mh --within C:\\PLINK\\str1.cluster2 --adjust --out C:\\PLINK\\aac1")

## --> Log file shows that 45 clusters were found (1 single individual, 44 pairs of individuals)

## To look at the adjusted results file: (stupidly long)
#shell("more C:\\PLINK\\aac1.cmh.adjusted")
## OR, BETTER:
aac1.cmh.adj <- read.table("C:\\PLINK\\aac1.cmh.adjusted", header=TRUE)
head(aac1.cmh.adj)

## --> OBS: the "disease" SNP (rs2222162) has moved from being 2nd to 1st ranked, 
## BUT remains insignificant after genome-wide correction... 

##################
## CLUSTERING 2 ##
##################

## NOTE: When CLUSTERING, above, we specified that PLINK should pair up the most similar individuals. 
## We can also perform clustering with DIFFERENT CONSTRAINTS:
## Eg. 
## Instead of imposing a MAXIMUM cluster size, 
## we can rather requiest that EACH cluster contains AT LEAST 1 case AND 1 control (s.t. it is informative for association)
## --cc option
## Also, we can specify a threshold of 0.01 for --ppc. 

shell("plink --bfile C:\\PLINK\\hapmap1 --cluster --cc --ppc 0.01 --out C:\\PLINK\\version2")
## --> new final solution "version2.cluster1"

## We can now re-run our association analysis with the new clustering scheme:
shell("plink --bfile C:\\PLINK\\hapmap1 --mh --within C:\\PLINK\\version2.cluster2 --adjust --out C:\\PLINK\\aac2")

## --> Log file now shows that FIVE clusters were found with a LOW inflation factor. 

## To look at the new adjusted results file:
aac2.cmh.adj <- read.table("C:\\PLINK\\aac2.cmh.adjusted", header=TRUE)
head(aac2.cmh.adj)

## --> OBS: Now our "disease" SNP is 1st-ranked AND has achieved genome-wide signficance (even w Bonf corr)!! 


##################
## CLUSTERING 3 ##
##################

## We can also perform the stratification analysis by SPECIFYING THE NUMBER OF CLUSTERS DESIRED. 
## Eg. Here we specify TWO final clusters (--K option), 
## and remove the significance test constraint ( ie. set "--ppc 0" by simply omitting that option)

shell("plink --bfile C:\\PLINK\\hapmap1 --cluster --K 2 --out C:\\PLINK\\version3")

## --> 2-class solution

## Using this clustering in our association analysis:
shell("plink --bfile C:\\PLINK\\hapmap1 --mh --within C:\\PLINK\\version3.cluster2 --adjust --out C:\\PLINK\\aac3")

## --> again, g-w significance :)

##################
## CLUSTERING 4 ##
##################

## Finally, given that in this case the ACTUAL ANCESTRY IS KNOWN, 
## we can use this EXTERNAL CLUSTERING in the analysis:

shell("plink --bfile C:\\PLINK\\hapmap1 --mh --within C:\\PLINK\\pop.phe --adjust --out C:\\PLINK\\aac3")

## --> Gives very similar results to 2-cluster soln (no surprise there)


##########################
## SUMMARY : CLUSTERING ##
##########################

## - Simple IBS-clustering works well wrt differentiating btw Chinese and Japanese individuals (with this n.SNPs)
## - Accounting for this pop structure can lower false positives rates AND increase power :) ...
##### ... the disease variant in only g-w sig AFTER performing a stratified analysis. 
## - Different approaches to clustering were applied. 

## - In general: When a SMALL NUMBER of DISCRETE SUB-POPULATIONS exist in the sample, ...
#### ... then a CLUSTER SOLUTION that most closely resembles this structure might be expected to work well. 

## - By contrast: IF, instead of a small number of discrete homogeneous clusters, the  sample acually contains ...
#### ... a COMPLEX MIXTURE OF INDIVUDALS from ACROSS a range of CLINES OF ANCESTRY, ...
#### ... then we might expect the approaches that form a LARGE NUMBER OF SMALLER CLUSTERS (eg. MATCHING PAIRS) to perform better.


###############################
## VISUALISATION: CLUSTERING ##
###############################

## SKIPPED ##


#################################
## QUANTITATIVE TRAIT ANALYSIS ##
#################################

## SKIPPED ##



##################################
## EXTRACTING A SNP OF INTEREST ##
##################################

## Having identified a SNP/ set of SNPs/ region of interest. you may want
## to extract those SNPs as a separate, smaller, more manageable file. 

## Will need to CONVERT from BINARY PED file --> STANDARD PED file format. 
## Use: --recode options
## Eg. --recodeAD --> GENERATES GENOTYPES CONVENIENT FOR SUBSEQUENT ANALYSIS IN R! 

## Extract a single SNP:
shell("plink --bfile C:\\PLINK\\hapmap1 --snp rs2222162 --recodeAD --out C:\\PLINK\\rec_snp1")

## READ INTO R:
#dat <- read.table("C:\\PLINK\\rec_snp1.recode.raw", header=TRUE) ## recode??
dat <- read.table("C:\\PLINK\\rec_snp1.raw", header=TRUE)

head(dat)

## Run a glm, get summary:
#summary(glm(PHENOTYPE-1 ~ rs2222162_A, data=dat, family="binomial")) ## failed to recode 1/2 into A/X??
summary(glm(PHENOTYPE-1 ~ rs2222162_1, data=dat, family="binomial"))

###






#####################     #####################     #####################     #####################     #####################







#############################################################################################
######## PLINK - R SUPPORT ##################################################################
#############################################################################################

## from: http://pngu.mgh.harvard.edu/~purcell/plink/rfunc.shtml

## PLINK expects a function in the EXACT form,
## to be defined in the supplied file:

Rplink <- function(PHENO, GENO, CLUSTER, COVAR){}

## PHENO # vector of phenotypes(n)
## GENO # matrix of genotypes (n x l)
## CLUSTER # vector of cluster membership coes (n)
## COVAR # matrix of covariates (n x c)

## Where: 
## n = the number of individuals in the analysis
## l = the number of SNPs*
## c = the number of covariates (if any)

## * In practice, "l", the number of SNPs, will probably be SMALLER than the ACTUAL number of SNPs in the file. 
## PLINK passes genotype data into R in batches, rather than all in one go. 

## GENOTYPES are CODED 0,1, or 2 copies of the MINOR allele, and NA, as per the --recodeA option. 

## For EACH SNP, the fn must return a numeric vector of values of the form: c(length(r), r)
## where "r" is the desired return vector. 
## ie. For each SNP, PLINK expects to get from the fn, c("how many values to read for SNPi", "values for SNPi")


## Function expected to RETURN a NUMERIC VECTOR with as many elements as there are SNPs. 

## Internally, PLINK will call the Rplink fn. 




## .. STOPPED .. ##








#####################     #####################     #####################     #####################     #####################







#############################################################################################
######## PLINK DATA FORMAT ##################################################################
#############################################################################################

## data formats # http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
## binary ped (bed) # http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

################
## CONVERTING ##
################

## ATTEMPT: CONVERT A COALESCENT.SIM + PHENO OUTPUT FILE --> .PED + .MAP FILES (to be input into PLINK) ##


## get sample geno data (from file...):

#plink.wd <- getwd()
setwd("C:/Cait 2012/Work/Xavier/Sims/set3/")

snps <- get(load("./set3_1_snps.Rdata"))
phen <- get(load("./set3_1_phen.Rdata"))

## inspect data
str(snps)
str(phen)


##################################
## convert snps --> ped format: ##
##################################

## PED FORMAT: ##

## 6 "MANDATORY" COLUMNS first:
## FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype
#### ... EXCLUDE OTHER "mandatory" columns w PLINK commands
#### --> (IndividualID, Phenotype) 
## 7-to-p: GENOTYPE COLUMNS
## NO HEADER ROW (belongs in .MAP file)


## SNPs can NOT be 0 --> convert from 0/1 to 1/2:
snps <- snps+1
row.names(snps) <- paste("ind", row.names(snps), sep=".")
## save row and column names for elsewhere:
individualID <- row.names(snps)
loci.names <- colnames(snps)

## convert to matrix:
snps <- matrix(snps, byrow=FALSE, ncol=ncol(snps))
## Replace every other column with a copy of the column before it:
snps[,seq(2, ncol(snps), 2)] <- snps[,seq(1, ncol(snps), 2)]


## add first "mandatory" columns (ie. individualID, phenotype)
#head(phen)
#head(as.numeric(phen))

## recode phen: S as 0, R as 1:
phen <- as.character(phen)
phen <- replace(phen, which(phen=="A"), 1)
phen <- replace(phen, which(phen=="B"), 2)
## check:
#head(phen)
#table(phen)

## bind individualID, phen, snps
dat <- cbind(individualID, phen, snps)
colnames(dat) <- NULL
#dat[1:10,1:10]

ped <- dat

## save dat.ped as Rdata
save(ped, file="C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_ped.Rdata")
#ped <- get(load("C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_ped.Rdata"))

## save as text!
#ped <- dat
write.table(ped, file="C:/PLINK/set3.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

## Do NOT save as .PED
## convert from text to PED!!
shell("mv C:\\PLINK\\set3.txt C:\\PLINK\\set3.ped")


##

####################################
## convert snp meta-data --> .map ##
####################################

## MAP format: ##

## EACH LINE of a .map file describes a SINGLE marker.
## Contains 4 "MANDATORY" COLUMNS:
## Chromosome, rs#/SNP identifier, (Genetic distance), Base-pair position.


## loci.names:
## SHOULD BE ONE NAME PER SITE (ie. PER TWO LOCI): ie. L001, NOT L001.1, L001.2 !!!!
head(loci.names) 
length(loci.names)

## remove decimal
#loci.names <- gsub("[.]", "", loci.names)
## remove last character(ie digit following decimal)
#loci.names <- substr(loci.names, 1, nchar(loci.names)-2)

## remove last TWO characters (ie. decimal and trailing digit):
loci.names <- substr(loci.names, 1, nchar(loci.names)-2)
## keep only every other:
loci.names <- loci.names[seq(1, length(loci.names), 2)]
## OR: # loci.names <- unique(loci.names)
head(loci.names)


chromosome <- rep(0, length(loci.names))
gen.dist <- rep(0, length(loci.names))
## get base-pair posi (ie. loci name - L):
bp <- loci.names
bp <- as.numeric(gsub("L", "", bp))
tail(bp)
dat <- data.frame(chromosome, loci.names, gen.dist, bp)

## as matrix, no header:
dat <- as.matrix(dat, byrow=FALSE, ncol=ncol(dat))
colnames(dat) <- NULL
head(dat)

map <- dat

## save as Rdata
save(map, file="C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_map.Rdata") 
#dat <- get(load("./set1_1_map.Rdata"))

## save properly as text:
write.table(map, file="C:/Cait 2012/Work/Xavier/Sims/set3/map.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(map, file="C:/PLINK/set3.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

## Do NOT save as .MAP
## convert from text to MAP!!
shell("mv C:\\PLINK\\set3.txt C:\\PLINK\\set3.map")



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


###################
## RUNNING PLINK ##
###################

## Set working directory to plink program 
setwd("C:/Program Files/plink-1.07-dos")

## Confirm you are in that directory
shell("cd") 

## open your ped and map files##
## ped options needed: --no-fid --no-parents --no-sex --missing-genotype NA --1
## map options needed: (--map3 (?) ) # ONLY if NO genetic_distance column in map file
shell("plink --file C:\\PLINK\\set3 --no-fid --no-parents --no-sex --allow-no-sex")

## make a binary PED file ##
## (provide the full path, not just the file name)
shell("plink --file C:\\PLINK\\set3 --make-bed --no-fid --no-parents --no-sex --allow-no-sex --out C:\\PLINK\\set3")

## use --bfile to work with the BINARY file 
# (same as --file, but loads the binary one and prints summary stats)
shell("plink --bfile C:\\PLINK\\set3 --no-fid --no-parents --no-sex --allow-no-sex")

# ## check freq of SNPs...?
# shell("plink --file C:\\PLINK\\set3  --no-fid --no-parents --no-sex --allow-no-sex --freq --out C:\\PLINK\\set3_freq")
# ## yay! 
# set3_freq <- read.table("C:/PLINK/set3_freq.frq", header=TRUE)
# head(set3_freq)


################################
## BASIC ASSOCIATION ANALYSIS ##
################################

## --assoc runs "basic association test based on comparing allele freqs in cases and controls"
## != --fisher

## inspect snps assoc from file:
# snps.assoc <- get(load("C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_performance.Rdata"))
# snps.assoc <- snps.assoc[[1]]
# snps.assoc
##

## perform a basic association analysis on the disease trait for all single SNPs
shell("plink --bfile C:\\PLINK\\set3 --assoc --counts --allow-no-sex --out C:\\PLINK\\aset3")

## to view the file you created, just read it in with R:
aset3 <- read.table("C:/PLINK/aset3.assoc", header=TRUE)
head(aset3)

## check
which.min(aset3$P) # yay! 

## get p.vals
p.vals <- aset3$P

## get sig ##

## Bonferonni ##
p.vals.bonf <- p.adjust(p.vals, "bonferroni")
p.thresh <- 0.05 # 0.01 # 0.001
p.sig <- which(p.vals.bonf < p.thresh)
length(p.sig)

## FDR ##
p.vals.fdr <- p.adjust(p.vals, "fdr")
p.thresh <- 0.001 # 0.05 # 0.01 # 0.001
p.sig <- which(p.vals.fdr < p.thresh)
length(p.sig)

#########
## Fisher test:
shell("plink --file C:\\PLINK\\set3 --no-fid --no-parents --no-sex  --allow-no-sex --fisher --out  C:\\PLINK\\aset3_fisher")

aset3_fisher <- read.table("C:/PLINK/aset3_fisher.assoc.fisher", header=TRUE)
head(aset3_fisher)
which.min(aset3_fisher$P) # yay! 


## get p.vals
p.vals <- aset3_fisher$P

## get sig ##
p.thresh <- 0.001 # 0.05 # 0.01 # 0.001

## Bonferonni ##
p.vals.bonf <- p.adjust(p.vals, "bonferroni")
p.sig <- which(p.vals.bonf < p.thresh)
length(p.sig)

## FDR ##
p.vals.fdr <- p.adjust(p.vals, "fdr")
p.sig <- which(p.vals.fdr < p.thresh)
length(p.sig)

#######
## TROUBLESHOOTING: FISHER TEST (PLINK) FINDING TOO MANY SNPS SIG...
## PROBLEM: FISHER TEST PLINK != FISHER TEST R! WHY??
p.vals <- sapply(c(1:ncol(snps)), function(e) 
              fisher.test(snps[,e], y=phen, 
              alternative="two.sided")$p.value) 

## get sig ##
p.thresh <- 0.001 # 0.05 # 0.01 # 0.001

## Bonferonni ##
p.vals.bonf <- p.adjust(p.vals, "bonferroni")
p.sig <- which(p.vals.bonf < p.thresh)
length(p.sig)

## FDR ##
p.vals.fdr <- p.adjust(p.vals, "fdr")
p.sig <- which(p.vals.fdr < p.thresh)
length(p.sig)

#

  
















#####################     #####################     #####################     #####################     #####################



#############################################################################################
###### NEXT UP ... (see below for where to START) ###########################################
#############################################################################################

## JUST BEFORE THIS NEXT COMMAND:



##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

# Hi,
# 
# I came here looking for an *answer* to this question, 
# so I only have so much to offer, but I think I managed 
# to get PLINK's initial steps to work in R using the shell function...
# 
# This is what worked for me:
# 
# ***NOT* in R:**
# 
# - Install PLINK and add its location to your PATH.
# - Download the example files from PLINK's tutorial 
#   (http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml) 
#   and put them in a folder whose path contains NO spaces 
#   (unless you know something I don't, in which case, space it up). 
# 
# **Then, in R:**
# 
# ## Set your working directory as the path to the PLINK program files: ##
# setwd("C:/Program Files/plink-1.07-dos")
# 
# ##  Use shell to check that you are now in the right directory: ##
# shell("cd")
# 
# ## At this point, the command "plink" should be at least be recognized
# # (though you may get a different error)
# shell("plink")
# 
# ## Open the PLINK example files ##
# # FYI mine are in "C:/PLINK/", so replace that accordingly...
# shell("plink --file C:\\PLINK\\hapmap1") 
# 
# ## Make a binary PED file ##
# # (provide the full path, not just the file name)
# shell("plink --file C:\\PLINK\\hapmap1 --make-bed --out C:\\PLINK\\hapmap1")
# 
# ... and so on.
# 
# That's all I've done so far. But with any luck, mirroring the structure 
# and general format of those lines of code should allow you to do 
# what you like with PLINK from within R. 
# 
# *Hope that helps! <br>*
# 
# PS. The PLINK output should just print in your R console when you run the lines above.
# 
# *All the best, <br>
# - CC.* 



