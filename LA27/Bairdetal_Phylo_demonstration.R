#a demonstration of methods used in Baird et al

#add a path that includes these files
#"Bairdetal_Phylo_functions.R"
#"Grasses_27sp_traits.csv"
#"sack27_combined.tre"
setwd( )

rm(list = ls())

#phylo libraries & source code with functions
library(ape)
library(phytools)
library(PHYLOGR)
library(nlme)
library(nortest)
source("Bairdetal_Phylo_functions.R")

#load phylogeny and traits, and make sure the two have matching "tip.label"

#phylogeny
tree27 <- read.nexus("sack27_combined.tre")
x11()
plotTree(tree27)
#this shows that tip label formatting is uneven
#clean up the tip labels to make life a bit easier
tree27$tip.label <- gsub("\\.", "_", tree27$tip.label)
t27tl_long_ind <- grep("\\w+_\\w+_", tree27$tip.label)
t27tl_long <- tree27$tip.label[t27tl_long_ind]
el23 <- function(.){ paste(.[2], .[3], sep = "_") }
t27tl_long <- sapply(strsplit(t27tl_long, "_"), el23)
tree27$tip.label <- replace(tree27$tip.label, t27tl_long_ind, t27tl_long)
#fix a spelling mistake
tree27$tip.label <- gsub("Erharta", "Ehrharta", tree27$tip.label)
#fix a species assignment that differs from the tree/latest nomenclature
tree27$tip.label <- gsub("Pennisetum_setaceaum", "Cenchrus_setaceus", tree27$tip.label)

plotTree(tree27)

#dataset
veins<-read.csv("Grasses_27sp_traits.csv")
#needs a set of tip.labels that match the tree
veins$tip.label <- veins$Species
veins$tip.label <- gsub(" ", "_", veins$tip.label)
veins$tip.label[!veins$tip.label %in% tree27$tip.label]
veins$tip.label <- gsub("_ssp._\\w+", "", veins$tip.label)
veins$tip.label %in% tree27$tip.label
all(veins$tip.label == tree27$tip.label)
#order of taxa currently does not match

#match up the order to prevent any surprises
veins$tip.label <- factor(veins$tip.label, levels = tree27$tip.label)
veins <- veins[order(veins$tip.label),]
all(veins$tip.label == tree27$tip.label)

#rename rows in the data
row.names(veins) <- veins$tip.label

#select trait columns with no NA values for further analysis
not.na <- unlist(lapply(veins, function(.){ all(!is.na(.)) }))
is.num <- unlist(lapply(veins, is.numeric))
#all columns no NA
no.na <- names(veins)[not.na]
head(veins[ , no.na])
#just traits that are not.na
no.na.t <- names(veins)[not.na & is.num]
no.na.t
#check this
head(veins[ , no.na.t])

#all trait headers irrespective of whether there is missing data
all.t <- names(veins)[is.num]

#note, there is missing data for some traits - these will require separate analysis
length(all.t)
length(no.na.t)
all.t[!all.t %in% no.na.t]

#produce a log[10]-transformed version of all the data 
log.veins <- data.frame(veins[ , c("tip.label", "C3C4")],
												log(veins[ , all.t], base = 10)
												)

#example ANOVA syntax/output for different combinations of data

#print output for ANOVA on a single trait
#warnings will arise because of the way corPagel is initiated
C3C4.ANOVA(trt = "LW", data = veins, tree = tree27)$summary

#for log transformed data
C3C4.ANOVA(trt = "LW", data = log.veins, tree = tree27)$summary

#ANOVA on a subset of traits
C3C4.ANOVA(trt = no.na.t[c(1:3)], data = veins, tree = tree27 )$summary

#example REGRESSION syntax/output for different combinations of data
#these will produce warnings linked with phyl.RMA slope test that is not important here

#single pair of traits
phylo.REGRESSION(no.na.t[1], no.na.t[2], data = veins, tree = tree27)$summary

#one trait against a pair of traits
phylo.REGRESSION(no.na.t[1], no.na.t[c(2:3)], data = veins, tree = tree27)$summary

#log transformed data
phylo.REGRESSION(no.na.t[1], no.na.t[c(2:3)], data = log.veins, tree = tree27)$summary

#note that self-self are trimmed out so the below only returns two rows
phylo.REGRESSION(no.na.t[1], no.na.t[c(1:3)], data = log.veins, tree = tree27)$summary

#trying to evaluate a trait with missing values will fail 
C3C4.ANOVA(trt = "CD2", data = veins, tree = tree27 )$summary

#so generate versions of the data and tree that drop the missing value
veinsCD2 <- veins[!is.na(veins$CD2), ]
dropCD2 <- as.character(veins$tip.label[!veins$tip.label %in% veinsCD2$tip.label])
treeCD2 <- drop.tip(tree27, dropCD2)

#now this works
C3C4.ANOVA(trt = "CD2", data = veinsCD2, tree = treeCD2)$summary

#same for regression
#this fails
phylo.REGRESSION("CD2", "LW", data = veins, tree = tree27)$summary
#this works
phylo.REGRESSION("CD2", "LW", data = veinsCD2, tree = treeCD2)$summary
