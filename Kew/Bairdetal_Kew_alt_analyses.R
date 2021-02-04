##Kew leaf trait dataset, climate traits and phylogeny
##multiple regression and AIC comparisons with approximately scaled,  and centred data
##alternative models are fit to check impacts of left skew in MAT data
##needs:
###"Bairdetal_Kew_functions.R"
###"data.newclimate.csv"
###"C3_C4_coding.txt"
###"Poaceae_hyp1_allC4_maxcred.txt"
##writes csv outputs, and saves the workspace image

setwd( )

rm(list = ls())

library(ape)
library(phytools)
library(PHYLOGR)
library(nlme)
library(nortest)

#load functions
source("Bairdetal_Kew_functions.R")

#This provided by Alec, already has tip.labels to match Stips below
Kew <- read.csv("data.newclimate.csv", row.names = 1)
names(Kew)
#for convenience later, rename the growing season variables
names(Kew) <- gsub("AvgTempGS", "GST", names(Kew))
names(Kew) <- gsub("AvgPrecipGS", "GSP", names(Kew))
#and renames MAT.20 to be more informative
names(Kew) <- gsub("MAT.20", "MATplus20", names(Kew))

#The Kew trait dataset is a subset of the Spriggs phylogeny
#and CH is a subset of the other sets
#so we need to clip the phylogeny to match the trait data
#before doing this, make sure that tip label formatting matches

#phylogeny tip labels
Stips <- read.table("C3_C4_coding.txt", header = TRUE)
head(Stips)
table(Stips$tip.label %in% Kew$tip.label)
#this looks right...

#merge the trait data with the phylogeny tip labels
Kew.Stips <- merge(Kew, Stips, all.x = TRUE)
#match the formatting in the C3C4 columns from the two sources and check they correspond
Kew.Stips$stateC3C4 <- replace(Kew.Stips$state.C3.0.C4.1.,
																Kew.Stips$state.C3.0.C4.1. == 0,
																"C3"
																)
Kew.Stips$stateC3C4 <- replace(Kew.Stips$stateC3C4,
																Kew.Stips$stateC3C4 == 1,
																"C4"
																)
all(Kew.Stips$C3C4 == Kew.Stips$stateC3C4)

#create a new object
#retaining only species common to the data and tree tip labels
Kew.Spriggs <- Kew.Stips[!is.na(Kew.Stips$stateC3C4), ]

#Import the Spriggs tree
Stree <- read.tree("Poaceae_hyp1_allC4_maxcred.txt")
head(Stree$tip.label)

#check how many of the species in the Kew dataset are also in the tree
intree <- Stree$tip.label %in% Kew.Spriggs$tip.label
table(intree)
#1752 species included in phylogeny also have data
#1843 species included in phylogeny don't have data

#trim down Stree to produce a Kew tree
Ktree <- drop.tip(Stree, Stree$tip.label[!intree])
length(Ktree$tip.label)

#reorder tip labels in Kew.Spriggs to match those in Ktree
Kew.Spriggs$tip.label <- factor(Kew.Spriggs$tip.label, levels = Ktree$tip.label)
Kew.Spriggs <- Kew.Spriggs[order(Kew.Spriggs$tip.label), ]
row.names(Kew.Spriggs) <- Kew.Spriggs$tip.label
head(Kew.Spriggs)

#set up transformed variables for analysis

#MAP / 50 to get similar scaling to MATplus20
Kew.Spriggs$MAPover50 <- Kew.Spriggs$MAP / 50
#GSP / 100 to get similar scaling to GST
Kew.Spriggs$GSPover100 <- Kew.Spriggs$GSP / 100

#log and center all the independent variables
#so that analysis can be interpreted as multiplicative
#keep the naming as above for convenience in plotting functions later
#alt combinations for analysis are:
#MAPover50*MATplus20
#GSPover100*GST
Kew.Spriggs$MAPover50 <- log(Kew.Spriggs$MAPover50, base = 10)
Kew.Spriggs$MAPover50 <- Kew.Spriggs$MAPover50 - mean(Kew.Spriggs$MAPover50)

Kew.Spriggs$MATplus20 <- log(Kew.Spriggs$MATplus20, base = 10)
Kew.Spriggs$MATplus20 <- Kew.Spriggs$MATplus20 - mean(Kew.Spriggs$MATplus20)

Kew.Spriggs$GSPover100 <- log(Kew.Spriggs$GSPover100, base = 10)
Kew.Spriggs$GSPover100 <- Kew.Spriggs$GSPover100 - mean(Kew.Spriggs$GSPover100)

#so that original GST values are retained in the dataframe
Kew.Spriggs$GSToriginal <- Kew.Spriggs$GST

Kew.Spriggs$GST <- log(Kew.Spriggs$GST, base = 10)
Kew.Spriggs$GST <- Kew.Spriggs$GST - mean(Kew.Spriggs$GST)

#also log the dependent variables
#again, keeping the original name for purposes of plotting later
#and copying the original data in case it is needed later
Kew.Spriggs$LLoriginal <- Kew.Spriggs$LL
Kew.Spriggs$LL <- log(Kew.Spriggs$LL, base = 10)
Kew.Spriggs$LWoriginal <- Kew.Spriggs$LW
Kew.Spriggs$LW <- log(Kew.Spriggs$LW, base = 10)
Kew.Spriggs$LAoriginal <- Kew.Spriggs$LA
Kew.Spriggs$LA <- log(Kew.Spriggs$LA, base = 10)
Kew.Spriggs$CHoriginal <- Kew.Spriggs$CH
Kew.Spriggs$CH <- log(Kew.Spriggs$CH, base = 10)

#traits for scaling analysis with n = 1752 are 
pl.trts <- c("LL", "LW", "LA")
#because
table(!is.na(Kew.Spriggs$CH))
#climate traits are
MA <- c("MAPover50", "MATplus20")
GS <- c("GSPover100", "GST")

#looksee
plot(Kew.Spriggs[, c(pl.trts, MA)])
x11()
plot(Kew.Spriggs[, c(pl.trts, GS)])

#run analyses
MA.1752 <- allmods(pl.trts, MA, Ktree, Kew.Spriggs)
write.csv(MA.1752$stats_table, "Bairdetal_Kew_alt_MA_1752.csv")
GS.1752 <- allmods(pl.trts, GS, Ktree, Kew.Spriggs)
write.csv(GS.1752$stats_table, "Bairdetal_Kew_alt_GS_1752.csv")

#Additional models testing
#relationship between climate and CH
#and
#the effects of removing left skewed MAT values
#require a different phylogeny because N is different

#CH
Kew.Spriggs.CH <- Kew.Spriggs[!is.na(Kew.Spriggs$CH), ]
intree.CH <- Stree$tip.label %in% Kew.Spriggs.CH$tip.label
table(intree.CH)
#1729 species included in phylogeny also have data
#1866 species included in phylogeny don't have data

#trim down Stree to produce a Kew tree for CH
Ktree.CH <- drop.tip(Stree, Stree$tip.label[!intree.CH])
length(Ktree.CH$tip.label)

#reorder tip labels in Kew.Spriggs.CH to match those in Ktree.CH
Kew.Spriggs.CH$tip.label <- factor(Kew.Spriggs.CH$tip.label, levels = Ktree.CH$tip.label)
Kew.Spriggs.CH <- Kew.Spriggs.CH[order(Kew.Spriggs.CH$tip.label), ]
row.names(Kew.Spriggs.CH) <- Kew.Spriggs.CH$tip.label
head(Kew.Spriggs.CH)

MA.CH <- allmods("CH", MA, Ktree.CH, Kew.Spriggs.CH)
write.csv(MA.CH$stats_table, "Bairdetal_Kew_alt_MA_CH.csv")
GS.CH <- allmods("CH", GS, Ktree.CH, Kew.Spriggs.CH)
write.csv(GS.CH$stats_table, "Bairdetal_Kew_alt_GS_CH.csv")

#removing left skewed MAT as MAT>0
Kew.Spriggs.MAT0 <- Kew.Spriggs[Kew.Spriggs$MAT >= 0, ]
intree.MAT0 <- Stree$tip.label %in% Kew.Spriggs.MAT0$tip.label
table(intree.MAT0)
#1723 species included in phylogeny also have data (29 species dropped)
Ktree.MAT0 <- drop.tip(Stree, Stree$tip.label[!intree.MAT0])
length(Ktree.MAT0$tip.label)

MA.MAT0 <- allmods(pl.trts, MA, Ktree.MAT0, Kew.Spriggs.MAT0)
write.csv(MA.MAT0$stats_table, "Bairdetal_Kew_alt_MA_MAT0.csv")
GS.MAT0 <- allmods(pl.trts, GS, Ktree.MAT0, Kew.Spriggs.MAT0)
write.csv(GS.MAT0$stats_table, "Bairdetal_Kew_alt_GS_MAT0.csv")

#for CH with these exclusions a further subset is needed
table(!is.na(Kew.Spriggs.MAT0$CH))

Kew.Spriggs.MAT0.CH <- Kew.Spriggs.MAT0[!is.na(Kew.Spriggs.MAT0$CH), ]
intree.MAT0.CH <- Stree$tip.label %in% Kew.Spriggs.MAT0.CH$tip.label
table(intree.MAT0.CH)
#1700 species
Ktree.MAT0.CH <- drop.tip(Stree, Stree$tip.label[!intree.MAT0.CH])
length(Ktree.MAT0.CH$tip.label)

MA.MAT0.CH <- allmods("CH", MA, Ktree.MAT0.CH, Kew.Spriggs.MAT0.CH)
write.csv(MA.MAT0.CH$stats_table, "Bairdetal_Kew_alt_MA_MAT0_CH.csv")
GS.MAT0.CH <- allmods("CH", GS, Ktree.MAT0.CH, Kew.Spriggs.MAT0.CH)
write.csv(GS.MAT0.CH$stats_table, "Bairdetal_Kew_alt_GS_MAT0_CH.csv")

save.image("Bairdetal_Kew_alt_analyses.Rdata")
