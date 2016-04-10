

## Campy data ##

#require(ape)

ST21.dna <- fasta2DNAbin("~/treeWAS/misc/Campy/ST21.fasta")
str(ST21.dna)

ST21.gen <- DNAbin2genind(ST21.dna)
str(ST21.gen@tab) # ncol = 94,711

###
phen <- read.table("~/treeWAS/misc/Campy/ST21-hosts.csv", header=FALSE, sep=",")
head(phen)
noms.phen <- phen$V1
phen <- phen$V2
rownames(snps)

prefix <- "~/treeWAS/misc/Campy/ST21.fasta.out"
tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
l<-length(seqs[1,])
l # l = 12,522

table(mapping)[1:10]

plot(tree)

dist <- readCFML(prefix=prefix, plot=TRUE)
dist

str(seqs)
tips <- which(rownames(seqs) %in% tree$tip.label)
snps <- DNAbin2genind(seqs[tips,], polyThres = 0)
snps <- snps@tab
str(snps) # ncol = 25,712
snps[1:10,1:10]

coln <- colnames(snps)
coln <- removeLastN(coln, 2)
coln <- as.numeric(coln)
length(unique(coln)) # 12,552
length(which(table(coln) > 2)) # 657

## WHAT TO DO ABOUT NON-BINARY SNPS??


